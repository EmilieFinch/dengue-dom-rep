# -----------------------------------------------------------------------------------------------------------
# Generate manuscript plots of dengue and climate data for Dominican Republic (2013-2023)
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))

# Create output folder
output_folder <- here("figures")

# Set up plot font
plot_font <- "Open Sans"
font_add_google(plot_font)
showtext_opts(dpi = 300)
showtext_auto()

# Read in data

## Shape files
shp_adm2 <- readRDS(here("data", "adm2-shp-clean.rds"))
shp_adm1 <- readRDS(here("data", "adm1-shp-clean.rds"))

## Cliamte and dengue data
adm2_data <- readRDS(here("data", "adm2-climate-epi.rds"))
adm1_data <- readRDS(here("data", "adm1-climate-epi.rds"))

grid <- read_xlsx(here("data", "adm1-grid.xlsx"))

# Figure 1 --------------------------------------------------------------------------------------------------------------

map_fig <- adm2_data |> 
  filter(date >= as.Date("2013-03-31") & date <= as.Date("2023-03-26")) |> 
  group_by(epiweek, adm2_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  ungroup() |> 
  mutate(year_index = case_when(epiweek == 53 ~ lag(year_index), T ~ year_index)) |>  # To fix week 53
  group_by(year_index, adm2_name) |> 
  summarise(
    total_reported = sum(reported),
    population     = population[which.min(date)],   # denominator = start-of-season year
    incidence      = total_reported / population * 1e5,
    .groups        = "drop") |>
  left_join(shp_adm2, by = "adm2_name") |> 
  group_by(adm2_name, geometry) |> 
  summarise(incidence = mean(incidence), .groups = "keep") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill = log(incidence + 1))) +
  scale_fill_distiller(name = "Mean\nannual \nDIR\n(log)", palette = "Oranges", direction = 1) +
  theme_void() +
  theme(
    legend.title = element_text(size = 14, family = plot_font), legend.text = element_text(size = 13, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.position = "right")

caribs <- st_read(here("data", "caribis", "caribis.shp"))

coords <- as.data.frame(st_centroid(caribs$geometry)) |> 
  separate(geometry, into = c("long", "lat"), sep = " ") |> 
  mutate(long = as.numeric(gsub("[^0-9.-]", "", long)), lat = as.numeric(gsub("[^0-9.-]", "", lat)))

caribs <- caribs |> 
  mutate(long = coords$long, lat = coords$lat)

plot_carib <- ggplot() +
  geom_sf(data = caribs, fill = "#969696",  color = "#969696", size = 0.2) +  
  geom_rect(aes(xmin = -72.20307, xmax = -68.12224, ymin = 17.2 , ymax =20.3), fill = NA, col = "#08519c", size = 0.4) +
  guides(x = "none", y = "none") +
  coord_sf() +
  labs(x = NULL, y = NULL) +
  theme_void() + theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.01,0.01)) +
  theme(panel.border = element_rect(color = "black", 
                              fill = NA, 
                              size = 0.7))

map_inset <- ggdraw() +
  draw_plot(map_fig) +
  draw_plot(plot_carib, x = 0.60, y = 0.69, width = 0.3, height = 0.3)

cases_ts <- adm1_data |> 
  filter(date >= as.Date("2013-03-31") & date <= as.Date("2023-03-26")) |> 
  group_by(date) |> 
  summarise(reported = sum(reported)) |> 
  ggplot() +
  geom_col(aes(x = date, y = reported), col = "#969696", width = 5) +
  labs(x = "Date", y = "Dengue cases") +
  theme_classic2() +
  scale_x_date(labels = label_date_short(), breaks = "1 years") +
  scale_y_continuous(breaks = c(0,200,400,600,800,1000)) +
  theme(
    axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 14, family = plot_font),
    axis.title.y = element_text(size = 14, family = plot_font), axis.title.x = element_text(size = 14, family = plot_font),
    legend.title = element_text(size = 14, family = plot_font), legend.text = element_text(size = 14, family = plot_font),
    legend.key.height = unit(1.2, "cm"))


grid <- grid |> 
  mutate(name = case_when(name == "Santiago Rodriguez" ~ "Santiago Rod.", 
                          name == "Hermanas Mirabal" ~ "Hermanas Mir.",
                          name == "Maria Trinidad Sanchez" ~ "Maria Trin. San.",
                          name == "San Pedro de Macoris" ~ "San Pedro de Mac.",
                          T ~ name))

incidence_tile <- adm1_data |> 
  filter(date >= as.Date("2013-04-01") & date <= as.Date("2023-03-26")) |> 
  mutate(month = month(date), incidence = reported/population*100000) |>
  mutate(adm1_name = case_when(adm1_name == "Santiago Rodriguez" ~ "Santiago Rod.", 
                          adm1_name == "Hermanas Mirabal" ~ "Hermanas Mir.",
                          adm1_name == "Maria Trinidad Sanchez" ~ "Maria Trin. San.",
                          adm1_name == "San Pedro de Macoris" ~ "San Pedro de Mac.",
                          T ~ adm1_name))  |> 
  group_by(year, month, adm1_name) |> 
  summarise(incidence = mean(incidence)) |> 
  group_by(month, adm1_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  group_by(adm1_name, year_index) |> 
  mutate(month_index = 1:n()) |> 
  ungroup() |> 
  ggplot(aes(x = month_index, y = as.factor(year_index), fill = log(incidence + 1))) + 
  geom_raster() +
  ylab(NULL) + 
  xlab(NULL) +
  scale_fill_gradientn(name = "Monthly\nDIR\n(log)", colours=brewer.pal(9, "Oranges")) +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Apr", "Jul", "Oct","Jan")) +
  scale_y_discrete(breaks = seq(1,10, by = 1), labels = c("2013/14", "", "2015/16", "", "2017/18", "", "2019/20", "", "2021/22", "")) +
  facet_geo(~ adm1_name, grid = grid) +
  #ggtitle("c") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"), 
        axis.text.x = element_text(size = 11, family = plot_font), 
        axis.text.y = element_text(size = 12, family = plot_font), 
        axis.title=element_blank(), 
        legend.title = element_text(size = 15, family = plot_font), 
        legend.position = "right",
        legend.key.size = unit(0.75, "cm"), 
        legend.text = element_text(size = 13, family = plot_font), 
        legend.spacing.x = unit(0.8, "cm"),
        strip.text.x = element_text(size = 9.5, family = plot_font))


incidence_tile_grob <- get_geofacet_grob(incidence_tile)

fig_1_a <- ggarrange(map_inset, cases_ts, nrow = 1, labels = c("a", "b"), font.label = list(family = plot_font, size = 24))

fig_1 <- ggarrange(fig_1_a, incidence_tile_grob, labels = c("", "c"), heights = c(1,1.55), font.label = list(family = plot_font, size = 28), nrow = 2)
ggsave(here(output_folder, "figure_1.png"), fig_1,  width = 450, height = 450, units = "mm", bg="white", dpi = 300)

# Supplementary Figures 1-3: climate facets -----------------------------------------------------------------------

tasmean_province <- adm1_data |> 
  filter(date >= as.Date("2013-04-01") & date <= as.Date("2023-03-26")) |> 
  mutate(month = month(date), incidence = reported/population*100000) |>
  mutate(adm1_name = case_when(adm1_name == "Santiago Rodriguez" ~ "Santiago Rod.", 
                               adm1_name == "Hermanas Mirabal" ~ "Hermanas Mir.",
                               adm1_name == "Maria Trinidad Sanchez" ~ "Maria Trin. San.",
                               adm1_name == "San Pedro de Macoris" ~ "San Pedro de Mac.",
                               T ~ adm1_name))  |> 
  group_by(year, month, adm1_name) |> 
  summarise(tas = mean(tas)) |> 
  group_by(month, adm1_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  group_by(adm1_name, year_index) |> 
  mutate(month_index = 1:n()) |> 
  ungroup() |> 
  ggplot(aes(x = month_index, y = as.factor(year_index), fill = tas)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("") +
  scale_fill_gradientn(name = "Mean \ntemperature (°C)", colours=brewer.pal(9, "OrRd")) +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Apr", "Jul", "Oct","Jan")) +
  scale_y_discrete(breaks = seq(1,10, by = 1), labels = c("2013/14", "", "2015/16", "", "2017/18", "", "2019/20", "", "2021/22", "")) +
  facet_geo(~ adm1_name, grid = grid) +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"), 
        axis.text.x = element_text(size = 11, family = plot_font),
        axis.text.y = element_text(size = 12, family = plot_font), 
        axis.title=element_blank(), legend.title = element_text(size = 15, family = plot_font), 
        legend.position = "right",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size = 12, family = plot_font), 
        legend.spacing.x = unit(0.8, "cm"),
        strip.text.x = element_text(size = 12, family = plot_font))

ggsave(here(output_folder, "sfig_1.png"), tasmean_province,  width = 480, height = 300, units = "mm", bg="white", dpi = 300)

prlr_province <-  adm1_data |> 
  filter(date >= as.Date("2013-04-01") & date <= as.Date("2023-03-26")) |> 
  mutate(month = month(date), incidence = reported/population*100000) |>
  mutate(adm1_name = case_when(adm1_name == "Santiago Rodriguez" ~ "Santiago Rod.", 
                               adm1_name == "Hermanas Mirabal" ~ "Hermanas Mir.",
                               adm1_name == "Maria Trinidad Sanchez" ~ "Maria Trin. San.",
                               adm1_name == "San Pedro de Macoris" ~ "San Pedro de Mac.",
                               T ~ adm1_name))  |> 
  group_by(year, month, adm1_name) |> 
  summarise(prlr = sum(prlr)) |> 
  group_by(month, adm1_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  group_by(adm1_name, year_index) |> 
  mutate(month_index = 1:n()) |> 
  ungroup() |> 
  ggplot(aes(x = month_index, y = as.factor(year_index), fill = prlr)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("") +
  scale_fill_gradientn(name = "Precipitation \n(mm)", colours=brewer.pal(9, "Blues")) +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Apr", "Jul", "Oct","Jan")) +
  scale_y_discrete(breaks = seq(1,10, by = 1), labels = c("2013/14", "", "2015/16", "", "2017/18", "", "2019/20", "", "2021/22", "")) +
  facet_geo(~ adm1_name, grid = grid) +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"),
        axis.text.x = element_text(size = 11, family = plot_font),
        axis.text.y = element_text(size = 12, family = plot_font), 
        axis.title=element_blank(), 
        legend.title = element_text(size = 15, family = plot_font), 
        legend.position = "right",
        legend.key.size = unit(0.75, "cm"), 
        legend.text = element_text(size = 12, family = plot_font), 
        legend.spacing.x = unit(0.8, "cm"), 
        strip.text.x = element_text(size = 12, family = plot_font))

ggsave(here(output_folder, "sfig_2.png"), prlr_province,  width = 480, height = 300, units = "mm", bg="white", dpi = 300)


relhum_province <- adm1_data |> 
  filter(date >= as.Date("2013-04-01") & date <= as.Date("2023-03-26")) |> 
  mutate(month = month(date), incidence = reported/population*100000) |>
  mutate(adm1_name = case_when(adm1_name == "Santiago Rodriguez" ~ "Santiago Rod.", 
                               adm1_name == "Hermanas Mirabal" ~ "Hermanas Mir.",
                               adm1_name == "Maria Trinidad Sanchez" ~ "Maria Trin. San.",
                               adm1_name == "San Pedro de Macoris" ~ "San Pedro de Mac.",
                               T ~ adm1_name))  |> 
  group_by(year, month, adm1_name) |> 
  summarise(hurs = mean(hurs)) |> 
  group_by(month, adm1_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  group_by(adm1_name, year_index) |> 
  mutate(month_index = 1:n()) |> 
  ungroup() |> 
  ggplot(aes(x = month_index, y = as.factor(year_index), fill = hurs)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("") +
  scale_fill_gradientn(name = "Relative \nhumidity (%)", colours=brewer.pal(9, "Purples")) +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Apr", "Jul", "Oct","Jan")) +
  scale_y_discrete(breaks = seq(1,10, by = 1), labels = c("2013/14", "", "2015/16", "", "2017/18", "", "2019/20", "", "2021/22", "")) +
  facet_geo(~ adm1_name, grid = grid) +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"), 
        axis.text.x = element_text(size = 11, family = plot_font), 
        axis.text.y = element_text(size = 12, family = plot_font), 
        axis.title=element_blank(), 
        legend.title = element_text(size = 15, family = plot_font),
        legend.position = "right",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size = 12, family = plot_font), 
        legend.spacing.x = unit(0.8, "cm"), 
        strip.text.x = element_text(size = 12, family = plot_font))

ggsave(here(output_folder, "sfig_3.png"), relhum_province,  width = 480, height = 300, units = "mm", bg="white", dpi = 300)

