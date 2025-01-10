library(tidyverse)
library(sf)
library(latex2exp)
library(ggpubr)

# Load Font
windowsFonts(
  Constantia = windowsFont('Constantia')
)

# Define constants
IN_DSN <- './Gwr.gdb'

# Define ids
QUARTERS <- paste0('Q', 1:4)
DATA_LABELS <- c('NonLocal', 'HousePrice')
FORMULA_LABELS <- c('NonLocal.Day', 'NonLocal.Night')

# Define templates
TEMPLATE_IN_LAYER <- 'Gwr_2021'

# Define the container
processed_datasets <- list()

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
      mutate(
        NonLocal = NonLocal * 100
      )
    processed_datasets[[quarter]][[data_label]] <- processed_dataset
  }
}

expression <- c(
  PolDay = TeX('$PM_{2.5}\\ Concentration\\ (\\mu g/m^3)$'),
  PolNight = TeX('$PM_{2.5}\\ Concentration\\ (\\mu g/m^3)$'),
  NonLocal = 'Non Local (%)'
)

plots <- list()
var_pairs <- list(
  c('NonLocal', 'PolDay'),
  c('NonLocal', 'PolNight')
)

for (quarter in QUARTERS) {
  for (pair in var_pairs) {
    dataset <- processed_datasets[[quarter]][[pair[1]]]
    plots[[paste(pair[1], pair[2], quarter)]] <- 
      ggplot(
        dataset, 
        aes(
          x = .data[[pair[1]]], 
          y = .data[[pair[2]]]
        )
      ) +
      geom_point(
        fill = '#f1b8b1', 
        color = '#f1b8b1', 
        alpha = 0.7, 
        size = 1
      ) +
      geom_smooth(
        method = 'lm', 
        formula = y ~ I(x ^ 3) + I(x ^ 2) + x, 
        se = FALSE,
        color = '#707070', 
        linewidth = 0.75
      ) +
      geom_smooth(
        method = 'lm', 
        formula = y ~ x, 
        se = FALSE, 
        color = '#707070', 
        linewidth = 0.75
      ) +
      labs(
        x = expression[pair[1]], 
        y = expression[pair[2]]
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(family = 'Constantia'),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      )
  }
}

ggarrange(
  plotlist = plots,
  ncol = 4,
  nrow = 2,
  labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'),
  font.label = list(size = 16, family = 'Constantia')
) %>%
 ggsave(
   'Figures/Figure S3.png', 
   ., 
   height = unit(6.5, 'cm'), 
   width = unit(13, 'cm')
  )  

