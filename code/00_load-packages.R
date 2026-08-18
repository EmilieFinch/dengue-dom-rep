# -----------------------------------------------------------------------------------------------------------
# Load all packages
# -----------------------------------------------------------------------------------------------------------

#devtools::install_github("nathoze/Rsero")

# Load pacman if needed
if (!require("pacman")) install.packages("pacman")

# Load all packages used in this repo
pacman::p_load(
  Rsero, showtext, here, dplyr, sf, ggplot2, geofacet, ggpubr, scales, MetBrewer, 
  readxl, lubridate, RColorBrewer, tidyr, cowplot, stringr, qs, scoringutils, 
  data.table, INLA, tibble, hydroGOF, withr, purrr, zoo
)

# For plotting
dir.create(here("figures"))

plot_font <- "Open Sans"
font_add_google(plot_font)
showtext_opts(dpi = 300)
showtext_auto()