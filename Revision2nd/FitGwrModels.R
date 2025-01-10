library(tidyverse)
library(openxlsx)
source('./GwrHelper.R')

option(warn = 1)

# Define constants
IN_DSN <- './Gwr.gdb'
LOG_FILE <- './Results/GwrSel.log'
GWR_SEL_IMAGE <- './Results/GwrSelStat.RData'
GWR_SEL_EXCEL <- './Results/GwrSel.xlsx'
GWR_MODEL_IMAGE <- './Results/GwrModel.RData'

# Define ids
QUARTERS <- paste0('Q', 1:4)
DATA_LABELS <- c('NonLocal', 'HousePrice')
FORMULA_LABELS <- c('NonLocal.Day', 'NonLocal.Night', 'HousePrice')
FORMULA_TYPE_LABELS <- c('Linear', 'Polynomial')

# Define templates
TEMPLATE_IN_LAYER <- 'Gwr_2021'
TEMPLATE_OUT_LAYER <- './Results/GwrModel_2021'

# Define the container
processed_datasets <- list()
list_of_nn_dist <- c()
list_of_listw <- list()
list_of_gwr.sel.stat <- list()

# Define formulas
formulas <- list(
  NonLocal.Day = list(
    Linear = PolDay ~ NonLocal,
    Polynomial = PolDay ~ I(NonLocal ^ 3) + I(NonLocal ^ 2) + NonLocal
  ),
  NonLocal.Night = list(
    Linear = PolNight ~ NonLocal,
    Polynomial = PolNight ~ I(NonLocal ^ 3) + I(NonLocal ^ 2) + NonLocal
  ),
  HousePrice = list(
    Linear = PolNight ~ HousePrice,
    Polynomial = PolNight ~ I(HousePrice ^ 3) + I(HousePrice ^ 2) + HousePrice
  )
)

# Start
layer <- paste0(TEMPLATE_IN_LAYER, QUARTERS[1])
grid_shape <- read_sf(dsn = IN_DSN, layer = layer) %>%
  dplyr::select(GridID)

for (quarter in QUARTERS) {
  # Read dataset
  layer <- paste0(TEMPLATE_IN_LAYER, quarter)
  grid <- read_sf(dsn = IN_DSN, layer = layer)
  
  # Process dataset
  for (data_label in DATA_LABELS) {
    # Drop missing value
    processed_dataset <- grid %>%
      dplyr::filter(!is.na(.data[[data_label]])) %>%
      as_Spatial()
    processed_datasets[[quarter]][[data_label]] <- processed_dataset
    
    # Find the distance between the nearest neighbors
    nn_dist <- dnearneight.find(processed_datasets[[quarter]][[data_label]])
    list_of_nn_dist <- c(list_of_nn_dist, nn_dist)
  }
}

# Get the maximum of the distance between the nearest neighbors
nn_dist.max <- max(list_of_nn_dist)

# Create lists of weights
for (quarter in QUARTERS) {
  for (data_label in DATA_LABELS) {
    processed_dataset <- processed_datasets[[quarter]][[data_label]]
    
    list_of_listw[[quarter]][[data_label]] <- 
      dnearneight.listw(processed_dataset, nn_dist.max)
  }
}

# Gwr Selector
log_file <- file(LOG_FILE)
sink(log_file, append = TRUE)
sink(log_file, append = TRUE, type = 'message')

for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    data_label <- str_split(formula_label, '\\.')[[1]][1]
      
    for (formula_type_label in FORMULA_TYPE_LABELS) {
      data <- processed_datasets[[quarter]][[data_label]]
      f <- formulas[[formula_label]][[formula_type_label]]
      listw <- list_of_listw[[quarter]][[data_label]]
      
      cat.data_and_formula(quarter, formula_label, formula_type_label)
      # Calculate
      list_of_gwr.sel.stat[[quarter]][[formula_label]][[formula_type_label]] <-
        gwr.sel.stat(data, f, listw)
    }
  }
}

sink()
sink(type = 'message')

# Test F1 and F2
for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    for (formula_type_label in FORMULA_TYPE_LABELS) {
      stats <- list_of_gwr.sel.stat[[quarter]][[formula_label]][[formula_type_label]]
      for (i in 1:dim(stats)[1]) {
        F1.p <- stats$F1.p
        F2.p <- stats$F2.p
        if (as.numeric(F1.p) > 0.05 | as.numeric(F2.p) > 0.05) {
          print(paste(quarter, formula_label, formula_type_label, i))
        }
      }
    }
  }
}

# Export to excel
sheets <- list()
for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    sheet_label <- paste(quarter, formula_label, sep = '_')
    sheets[[sheet_label]] <- format_gwr_sel(
      list_of_gwr.sel.stat[[quarter]][[formula_label]][['Linear']],
      list_of_gwr.sel.stat[[quarter]][[formula_label]][['Polynomial']]
    )
  }
}
write.xlsx(sheets, file = GWR_SEL_EXCEL)

# Store image
rm(data, grid, listw, processed_dataset, stats, f, layer, nn_dist, F1.p, F2.p,
   data_label, formula_label, formula_type_label, quarter, i, sheet_label)
save.image(GWR_SEL_IMAGE)

# Fit Gwr models
list_of_gwr_paras <- list(
  Q1 = list(
    NonLocal.Day = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    NonLocal.Night = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    HousePrice = list(adapt = TRUE, method = 'aic', gweight = gwr.Gauss)
  ),
  Q2 = list(
    NonLocal.Day = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    NonLocal.Night = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    HousePrice = list(adapt = TRUE, method = 'cv', gweight = gwr.Gauss)
  ),
  Q3 = list(
    NonLocal.Day = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    NonLocal.Night = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    HousePrice = list(adapt = TRUE, method = 'cv', gweight = gwr.Gauss)
  ),
  Q4 = list(
    NonLocal.Day = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    NonLocal.Night = list(adapt = FALSE, method = 'aic', gweight = gwr.Gauss),
    HousePrice = list(adapt = TRUE, method = 'cv', gweight = gwr.Gauss)
  )
)

list_of_gwr_model <- list()
for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    cat.data_and_formula(quarter, formula_label, 'Linear')
    
    data_label <- str_split(formula_label, '\\.')[[1]][1]
    
    data <- processed_datasets[[quarter]][[data_label]]
    f <- formulas[[formula_label]][['Linear']]
    listw <- list_of_listw[[quarter]][[data_label]]
    
    adapt <- list_of_gwr_paras[[quarter]][[formula_label]][['adapt']]
    method <- list_of_gwr_paras[[quarter]][[formula_label]][['method']]
    gweight <- list_of_gwr_paras[[quarter]][[formula_label]][['gweight']]
    
    bandwidth <- gwr.sel(
      formula = f,
      data = data,
      adapt = adapt,
      method = method,
      gweight = gweight
    )
    
    if (adapt) {
      list_of_gwr_model[[quarter]][[formula_label]] <- gwr(
        formula = f, 
        data = data,
        adapt = bandwidth, 
        hatmatrix = T
      )
    } else {
      list_of_gwr_model[[quarter]][[formula_label]] <- gwr(
        formula = f, 
        data = data,
        bandwidth = bandwidth, 
        hatmatrix = T
      )
    }
  }
}

# Calculate Significance
for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    data_label <- str_split(formula_label, '\\.')[[1]][1]
    
    gwr_model <- list_of_gwr_model[[quarter]][[formula_label]]
    original_data <- processed_datasets[[quarter]][[data_label]]
    positive_sig <- ifelse(data_label == 'NonLocal', TRUE, FALSE)
    varname <- data_label
    
    list_of_gwr_model[[quarter]][[formula_label]] <- gwr.sig(
      gwr_model, 
      original_data, 
      positive_sig, 
      varname
    )
  }
}

# Merge quarter results
list_of_gwr_sig_sum <- list()
for (formula_label in FORMULA_LABELS) {
  quarter_merged <- data.frame()
  positive <- ifelse(formula_label == 'HousePrice', FALSE, TRUE)
  for (quarter in QUARTERS) {
    sp <- list_of_gwr_model[[quarter]][[formula_label]]$SDF
    extracted_gwr_sig <- extract_gwr_sig(sp, positive)
    extracted_gwr_sig[['quarter']] <- quarter
    quarter_merged <- rbind(quarter_merged, extracted_gwr_sig)
  }
  if (positive) {
    quarter_merged <- quarter_merged %>%
      group_by(GridID) %>%
      summarise(
        sig_count = sum(positive_sig)
      )
  } else {
    quarter_merged <- quarter_merged %>%
      group_by(GridID) %>%
      summarise(
        sig_count = sum(negative_sig)
      )
  }
  list_of_gwr_sig_sum[[formula_label]] <- grid_shape %>%
    full_join(quarter_merged, by = 'GridID')
}

# Export to gdb
for (formula_label in FORMULA_LABELS) {
  json_file <- paste0(TEMPLATE_OUT_LAYER, '_Merged_', formula_label, '.geojson')
  quarter_merged <- list_of_gwr_sig_sum[[formula_label]]
  write_sf(quarter_merged, json_file)
}

for (quarter in QUARTERS) {
  for (formula_label in FORMULA_LABELS) {
    json_file <- paste0(TEMPLATE_OUT_LAYER, quarter, '_', formula_label, '.geojson')
    gwr_model <- list_of_gwr_model[[quarter]][[formula_label]]$SDF %>%
      st_as_sf()
    write_sf(gwr_model, json_file)
  }
}


# Store image
rm(original_data, gwr_model, listw, f, nn_dist, varname, 
   positive_sig, adapt, bandwidth, method, gweight, 
   data_label, formula_label, formula_type_label, quarter)
save.image(GWR_MODEL_IMAGE)

# Plot

processed_datasets[['Q1']][['NonLocal']]@data %>%
  ggplot(aes(x = NonLocal, y = PolNight)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_smooth(method = 'lm')

