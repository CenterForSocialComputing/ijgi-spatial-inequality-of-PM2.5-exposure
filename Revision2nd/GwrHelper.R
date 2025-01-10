library(spatialreg)
library(spdep)
library(sf)
library(spgwr)

TEMPLATE_GWR_SEL_TABLE <- data.frame(
  adapt = as.character(c(
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)),
  method = c(
    'cv', 'cv', 'cv', 'aic', 'aic', 'aic',
    'cv', 'cv', 'cv', 'aic', 'aic', 'aic'
  ),
  gweight = c(
    'gwr.Gauss', 'gwr.bisquare', 'gwr.tricube',
    'gwr.Gauss', 'gwr.bisquare', 'gwr.tricube',
    'gwr.Gauss', 'gwr.bisquare', 'gwr.tricube',
    'gwr.Gauss', 'gwr.bisquare', 'gwr.tricube'
  )
)
 
dnearneight.find <- function(
    spatial.polygon
) {
  crds <- coordinates(spatial.polygon)
  
  max.nn <- knn2nb(knearneigh(crds, k = 1))
  max.nn.dist <- max(unlist(nbdists(max.nn, crds)))
  
  return (max.nn.dist)
}

dnearneight.listw <- function(
    spatial.polygon,
    max.nn.dist
) {
  crds <- coordinates(spatial.polygon)
  nb2listw(dnearneigh(crds, 0, max.nn.dist))
}

gwr.sel.stat <- function(data, formula, listw) {
  results <- data.frame(
    adapt = logical(),
    method = character(),
    gweight = character(),
    gwr.moran = numeric(), 
    gwr.moran.p = numeric(),
    lm.moran = numeric(), 
    lm.moran.p = numeric(), 
    F1 = numeric(), 
    F1.p = numeric(), 
    F2 = numeric(), 
    F2.p = numeric(), 
    ssr.ols = numeric(), 
    ssr.gwr = numeric(), 
    improve.gwr = numeric()
  )
  for (adapt in c(FALSE, TRUE)) {
    for (method in c('cv', 'aic')) {
      for (gweight in c('gwr.Gauss', 'gwr.bisquare', 'gwr.tricube')) {
        cat.gwr.sel(adapt, method, gweight)
        tryCatch(
          expr = {
            result <- gwr.gof(
              formula = formula, 
              data = data,
              listw = listw,
              adapt = adapt,
              method = method,
              gweight = gweight
            )
            results <- rbind(results, result)
          },
          error = function(e) {
            cat('Error, pass.\n')
          }
        )
      }
    }
  }
  colnames(results) <- c(
    'adapt', 'method', 'gweight',
    'gwr.moran','gwr.moran.p', 'lm.moran', 'lm.moran.p', 
    'F1', 'F1.p', 'F2', 'F2.p', 
    'ssr.ols', 'ssr.gwr', 'improve.gwr'
  )
  return(results)
}

gwr.gof <- function(
    formula, 
    data,
    listw,
    adapt = FALSE, 
    method = 'cv', 
    gweight = 'gwr.Gauss',
    show.error.messages = TRUE
){
  if (gweight == 'gwr.Gauss') {
    gweight_method <- gwr.Gauss
  }
  else if (gweight == 'gwr.bisquare') {
    gweight_method <- gwr.bisquare
  }
  else {
    gweight_method <- gwr.tricube
  }
  sel <- gwr.sel(
    formula = formula,
    data = data,
    adapt = adapt,
    method = method,
    gweight = gweight_method,
    show.error.messages = show.error.messages
  )
  
  if (adapt) {
    fit <- gwr(
      formula = formula, 
      data = data,
      adapt = sel, 
      hatmatrix = T
    )
  }
  else {
    fit <- gwr(
      formula = formula, 
      data = data,
      bandwidth = sel, 
      hatmatrix = T
    )
  }
  
  gwr.moran <- gwr.morantest(fit, listw)
  lm.moran <- lm(formula, data = data) %>%
    lm.morantest(listw)
  F1 <- LMZ.F1GWR.test(fit)
  F2 <- LMZ.F2GWR.test(fit)
  result <- c(
    adapt = adapt,
    method = method,
    gweight = gweight,
    gwr.moran = gwr.moran[["estimate"]] %>% unname(), 
    gwr.moran.p = gwr.moran[["p.value"]] %>% unname(),
    lm.moran = lm.moran[["estimate"]][1] %>% unname(), 
    lm.moran.p = lm.moran[["p.value"]] %>% unname(), 
    F1 = F1[["statistic"]] %>% unname(), 
    F1.p = F1[["p.value"]] %>% unname(), 
    F2 = F2[["statistic"]] %>% unname(), 
    F2.p = F2[["p.value"]] %>% unname(), 
    ssr.ols = F1[["estimates"]][["SS OLS residuals"]] %>% unname(), 
    ssr.gwr = F1[["estimates"]][["SS GWR residuals"]] %>% unname(), 
    improve.gwr = F2[["estimates"]][["SS GWR improvement"]] %>% unname()
  )
  return(result)
}

gwr.sig <- function(
    gwr_model, 
    original_data, 
    positive_sig, 
    varname
) {
  gwr_model$SDF@data <- mutate(
    gwr_model$SDF@data,
    GridID = original_data@data$GridID,
    pseudo.t = .data[[varname]] / .data[[paste0(varname, '_se')]],
    p.value = 2 * pt(0 - abs(pseudo.t), gwr_model$results$edf),
    sig = ifelse(
      rep(positive_sig, dim(gwr_model$SDF@data)[1]),
      ifelse(p.value < 0.05 & pseudo.t > 0, 1, 0),
      ifelse(p.value < 0.05 & pseudo.t < 0, 1, 0)
    )
  )
  
  colnames(gwr_model$SDF@data)[dim(gwr_model$SDF@data)[2]] <-
    ifelse(positive_sig, 'positive_sig', 'negative_sig')
  
  gwr_model
}

extract_gwr_sig <- function(sp, positive) {
  sig_name <- ifelse(positive, 'positive_sig', 'negative_sig')
  st_as_sf(sp) %>%
    st_drop_geometry() %>%
    dplyr::select(GridID, any_of(sig_name))
}

format_gwr_sel <- function(
    linear, 
    poly
) {
  linear.ols.ssr <- paste(
    sprintf('%0.0f', as.numeric(linear[1, 'ssr.ols'])),
    'for OLS (Linear)'
  )
  linear.ols.moran <- as.numeric(linear[1, 'lm.moran'])
  linear.ols.moran.str <- paste(
    sprintf('%0.3f', linear.ols.moran),
    'for OLS (Linear)'
  )
  
  poly.ols.ssr <- paste(
    sprintf('%0.0f', as.numeric(poly[1, 'ssr.ols'])),
    'for OLS (Polynomial)'
  )
  poly.ols.moran <- as.numeric(poly[1, 'lm.moran'])
  poly.ols.moran.str <- paste(
    sprintf('%0.3f', poly.ols.moran),
    'for OLS (Polynomial)'
  )
  
  linear <- linear %>%
    dplyr::select(
      adapt, 
      method, 
      gweight,
      gwr.moran,
      gwr.moran.p,
      ssr.gwr
    ) %>%
    mutate(
      gwr.moran = paste0(
        sprintf('%0.3f', as.numeric(gwr.moran)),
        sig_to_star(as.numeric(gwr.moran.p))
      ),
      ssr.gwr = sprintf('%0.3f', as.numeric(ssr.gwr))
    )
  
  poly <- poly %>%
    dplyr::select(
      adapt, 
      method, 
      gweight,
      gwr.moran,
      gwr.moran.p,
      ssr.gwr
    ) %>%
    mutate(
      gwr.moran = paste0(
        sprintf('%0.3f', as.numeric(gwr.moran)),
        sig_to_star(as.numeric(gwr.moran.p))
      ),
      ssr.gwr = sprintf('%0.3f', as.numeric(ssr.gwr))
    )
  
  format_result <- TEMPLATE_GWR_SEL_TABLE %>%
    full_join(
      linear, 
      by = c('adapt', 'method', 'gweight')
    ) %>%
    full_join(
      poly, 
      by = c('adapt', 'method', 'gweight'),
      suffix = c('.linear', '.poly')
    ) %>%
    dplyr::select(
      adapt, 
      method, 
      gweight,
      gwr.moran.linear,
      gwr.moran.poly,
      ssr.gwr.linear,
      ssr.gwr.poly
    )
  
  colnames(format_result) <- c(
    'adapt', 
    'method', 
    'gweight',
    linear.ols.moran.str,
    poly.ols.moran.str,
    linear.ols.ssr,
    poly.ols.ssr
  )
  return(format_result)
}

sig_to_star <- function(
    sig
) {
  ifelse(
    sig > 0.05,
    '',
    ifelse(
      sig > 0.01,
      '*',
      ifelse(
        sig > 0.001,
        '**',
        '***'
      )
    )
  )
}

cat.gwr.sel <- function(
    adapt, 
    method, 
    gweight
) {
  cat('------------------------------------------------\n')
  cat(paste('adapt: ', adapt, '; method: ', method, '; gweight: ', gweight, '\n'))
}

cat.data_and_formula <- function(
    quarter, 
    formula_label, 
    formula_type_label
) {
  cat('================================================\n')
  cat(paste('quarter: ', quarter, '; formula: ', formula_label, 
            '; formula type: ', formula_type_label, '\n'))
  cat('================================================\n')
}
