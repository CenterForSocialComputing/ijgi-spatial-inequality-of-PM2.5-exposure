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
QUARTERS <- c('Annual')
DATA_LABELS <- c('NonLocal', 'HousePrice')
FORMULA_LABELS <- c('NonLocal.Day', 'NonLocal.Night', 'HousePrice')

# Define templates
TEMPLATE_IN_LAYER <- 'Gwr_2021'

# Define the container
processed_datasets <- list()

# Start
layer <- TEMPLATE_IN_LAYER
grid_shape <- read_sf(dsn = IN_DSN, layer = layer) %>%
  dplyr::select(GridID)

# Read dataset
grid <- read_sf(dsn = IN_DSN, layer = layer)

# Process dataset
for (data_label in DATA_LABELS) {
  # Drop missing value
  processed_dataset <- grid %>%
    dplyr::filter(!is.na(.data[[data_label]])) %>%
    mutate(
      NonLocal = NonLocal * 100
    )
  processed_datasets[[data_label]] <- processed_dataset
}

expression <- c(
  PolDay = TeX('$PM_{2.5}\\ Concentration\\ (\\mu g/m^3)$'),
  PolNight = TeX('$PM_{2.5}\\ Concentration\\ (\\mu g/m^3)$'),
  NonLocal = 'Non Local (%)',
  HousePrice = 'House Price (CNY)'
)

plots <- list()
var_pairs <- list(
  c('NonLocal', 'PolDay'),
  c('NonLocal', 'PolNight'),
  c('HousePrice', 'PolNight')
)

for (quarter in QUARTERS) {
  for (pair in var_pairs) {
    dataset <- processed_datasets[[pair[1]]]
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
        axis.text = element_text(size = 10),
        plot.margin = grid::unit(c(5.5, 16.5, 5.5, 5.5), "points")
      )
  }
}

ggarrange(
  plotlist = plots,
  ncol = 3,
  labels = c('A', 'B', 'C'),
  font.label = list(size = 16, family = 'Constantia')
) %>%
 ggsave(
   'Figures/Figure 5.png', 
   ., 
   height = unit(3.25, 'cm'), 
   width = unit(9.75, 'cm')
  )  


