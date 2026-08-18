# -----------------------------------------------------------------------------------------------------------
# Visualise selected INLA model output
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))
source(here("code", "inla-loop_fn.R"))
source(here("code", "utils_fn.R"))
source(here("code", "create-lagged-data_fn.R"))

model_output <- qread(here("output", "selected-models-output.qs"))

# Prepare data
adm2_data <- readRDS(here("data", "adm2-climate-epi.rds"))
df_model <- lag_data(adm2_data)
shp_adm2 <- readRDS(here("data", "adm2-shp-clean.rds"))

df_model <- df_model |> 
  arrange(adm1_name, adm2_name, date) |> 
  group_by(adm2_name) |> 
  mutate(adm2_id = cur_group_id()) |> 
  ungroup() |> 
  mutate(month = month(date), incidence = reported/population*100000) |> 
  filter(date >= as.Date("2013-03-31") & date <= as.Date("2023-03-26")) |> 
  group_by(epiweek, adm2_name) |>
  # Define dengue season as running from epiweek 14 (beginning of April) until epiweek 13 of the following year
  mutate(year_index = 1:n()) |> 
  ungroup() |> 
  mutate(year_index = case_when(epiweek == 53 ~ lag(year_index), T ~ year_index)) |> # To fix week 53
  arrange(date, adm2_name) |> 
  group_by(adm2_name) |> 
  mutate(date_index = row_number()) |> 
  ungroup() |> 
  group_by(year_index, adm2_name) |> 
  mutate(cumulative_reported = cumsum(reported), estimated_infections = (cumulative_reported * 1/0.0747), estimated_prop_inf = estimated_infections/ population, estimated_sus = 1 - estimated_prop_inf) |> 
  mutate(estimated_sus = case_when(estimated_sus > 1 ~ 1L,
                                   estimated_sus < 0 ~ 0L,
                                   T ~ estimated_sus)) |> 
  group_by(adm2_name) |> 
  mutate(susceptibility = lag(estimated_sus),
         lag_weight_cases = weight_cases(reported),
         log_susceptibility = log(susceptibility + 0.001), 
         log_weight_cases = log(lag_weight_cases + 0.001)) |> 
  ungroup()

# Add outbreak thresholds --------------------------------------------------------------------------------

## Use a seasonal moving 75th percentile by adm2
## For a given month, year and adm we define an outbreak threshold of the 75th percentile of weekly cases in that month
## using all years up to (but not including) the given year

year_month <- df_model |>
  group_by(year_index, month) |>
  filter(year_index >= 5) |>
  summarise(.groups = "keep")

thresholds <- data.frame()
for(adm_id in 1: max(df_model$adm2_id)){ 
  
  threshold_adm <- purrr::map2_df(year_month$month, year_month$year_index, calculate_thresholds, data_input = df_model,adm_input = adm_id, quantile = 0.75, .progress = TRUE)
  thresholds <- rbind(thresholds, threshold_adm)
}
df_model <- df_model |> left_join(thresholds, by = c("year_index", "month", "adm2_id"))

# Here we define an outbreak week where the number of cases is > seasonal moving 75th percentile
# threshold, and an outbreak year as having more than 12 outbreak weeks

df_model <- df_model |>
  mutate(outbreak_week = case_when(reported > threshold ~ 1, TRUE ~ 0)) |> 
  mutate(outbreak_year = case_when(year_index == 1 | year_index == 3 | year_index == 7 | year_index == 10 ~ 1,
                                   T ~ 0))


# Figure 3 - effect estimate plot
mod_num <- 5
cat(paste0("Generating results for formula: ", model_output$adeq_stats$form[mod_num]))

nino_effect <- ggplot(model_output$random_effects[[mod_num]]$`inla.group(nino34_12_wk_avg_6, n = 9)`) + 
  geom_line(aes(x=ID, y = exp(mean)), col = "#3E7CA0") +
  geom_ribbon(aes(ymin= exp(`0.025quant`), ymax =exp(`0.975quant`), x =ID), alpha = 0.4, fill = "#3E7CA0") +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_y_continuous(limits = c(0,2.5), breaks = seq(0 ,2.5, by = 0.5)) +
  scale_x_continuous(limits = c(-1.7, 2.45)) +
  labs(x = NULL, y = "Relative risk") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10, family = plot_font), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

nino_density_plot <- ggplot(df_model) +
  geom_density(aes(y = nino34_12_wk_avg_6), fill = "#3E7CA0", alpha = 0.4, size = 0.1) +
  coord_flip() +
  labs(y = "\nNiño 3.4 SSTA (12 week average with 6 week lag)", x = "Density") +
  scale_y_continuous(limits = c(-1.7, 2.65)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 10, family = plot_font, color = "white"), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

nino_plots <- ggarrange(nino_effect, nino_density_plot, nrow = 2)

temp_effect <- ggplot(model_output$random_effects[[mod_num]]$`inla.group(tasmax_scale_8_wk_avg_4)`) + 
  geom_line(aes(x=ID + mean(df_model$tasmax, na.rm = TRUE), y = exp(mean)), col = "#3E7CA0") +
  geom_ribbon(aes(ymin= exp(`0.025quant`), ymax =exp(`0.975quant`), x =ID + mean(df_model$tasmax)), alpha = 0.4, fill = "#3E7CA0") +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_x_continuous(limits = c(18, 36)) +
  scale_y_continuous(limits = c(0,2.7), breaks = seq(0 ,2.5, by = 0.5)) +
  labs(x = NULL, y = "Relative risk") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10, family = plot_font), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

temp_density_plot <- ggplot(df_model) +
  geom_density(aes(y = tasmax_scale_8_wk_avg_4  + mean(df_model$tasmax)), fill = "#3E7CA0", alpha = 0.4, size = 0.1) +
  coord_flip() +
  scale_y_continuous(limits = c(18, 36)) +
  scale_x_continuous(limits = c(0,0.23), breaks = seq(0 ,0.2, by = 0.1)) +
  labs(y = "\n Maximum temperature (°C, 8 week average with 4 week lag)", x = "Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 10, family = plot_font, color = "white"), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

temp_plots <- ggarrange(temp_effect, temp_density_plot, nrow = 2)

hurs_effect <- ggplot(model_output$random_effects[[mod_num]]$`inla.group(hurs_scale_10_wk_avg_2)`) + 
  geom_line(aes(x=ID + mean(df_model$hurs, na.rm = TRUE), y = exp(mean)), col = "#3E7CA0") +
  geom_ribbon(aes(ymin= exp(`0.025quant`), ymax =exp(`0.975quant`), x =ID + mean(df_model$hurs)), alpha = 0.4, fill = "#3E7CA0") +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(x = NULL, y = "Relative risk") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10, family = plot_font), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

hurs_density_plot <- ggplot(df_model) +
  geom_density(aes(y = hurs_scale_10_wk_avg_2  + mean(df_model$hurs)), fill = "#3E7CA0", alpha = 0.4, size = 0.1) +
  coord_flip() +
  labs(y = "\n Relative humidity (%, 10 week average with 2 week lag)", x = "Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 10, family = plot_font, color = "white"), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

hurs_plots <- ggarrange(hurs_effect, hurs_density_plot, nrow = 2)

spi_effect <- ggplot(model_output$random_effects[[mod_num]]$`inla.group(spi_53)`) + 
  geom_line(aes(x=ID , y = exp(mean)), col = "#3E7CA0") +
  geom_ribbon(aes(ymin= exp(`0.025quant`), ymax =exp(`0.975quant`), x =ID), alpha = 0.4, fill = "#3E7CA0") +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(x = NULL, y = "Relative risk") +
  scale_y_continuous(limits = c(0,1.75), breaks = seq(0 ,1.75, by = 0.5)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10, family = plot_font), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

spi_density_plot <- ggplot(df_model) +
  geom_density(aes(y = spi_53), fill = "#3E7CA0", alpha = 0.4, size = 0.1) +
  coord_flip() +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0,0.5, by = 0.1)) +
  labs(y = "\n SPI (12 month timescale)", x = "Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 10, family = plot_font, color = "white"), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 10, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font))

spi_plots <- ggarrange(spi_effect, spi_density_plot, nrow = 2)
climate_plots <- ggarrange(nino_plots, temp_plots, hurs_plots, spi_plots, labels = c("a", "b", "c", "d"), font.label = list(family = plot_font, size = 18))

# Figure 3 - changes in random effects
spatial_re <- wrangle_random(model_output$random_effects[[mod_num]]$adm2_id, "spatial", is_bym = TRUE) |> 
  mutate(model = "climate_foi") |> 
  rbind(wrangle_random(model_output$random_effects[[2]]$adm2_id, "spatial", is_bym = TRUE) |> mutate(model = "climate")) |> 
  rbind(wrangle_random(model_output$random_effects[[1]]$adm2_id, "spatial", is_bym = TRUE) |> mutate(model = "baseline")) |> 
  rbind(wrangle_random(model_output$random_effects[[6]]$adm2_id, "spatial", is_bym = TRUE) |> mutate(model = "mechanistic")) |> 
  rename(year_index = group) |>
  left_join(df_model |> dplyr::select(adm2_name, adm2_id, year_index, outbreak_year) |> mutate(adm2_id = as.numeric(adm2_id)) |> rename(value = adm2_id) |> unique()) |> 
  mutate(model = factor(model, levels = c("baseline", "climate", "mechanistic", "climate_foi"), labels = c("Baseline", "Climate-only", "Epidemic-only", "Full model")))

# Average yearly random effect
year_re <- spatial_re |> 
  filter(component == "uv_joint") |> 
  group_by(year_index, model) |> 
  summarise(median_summary = median(median), lower = quantile(median, 0.025, na.rm = TRUE), upper = quantile(median, 0.975)) |> 
  ggplot() +
  geom_col(aes(x = year_index + 2012, y = median_summary, fill = model), position = "dodge") +
  scale_x_continuous(breaks = seq(2013,2023, by = 1)) + 
  labs(x = "Dengue season", y = "Marginal effect", fill = "") +
  scale_fill_manual(values = brewer.pal(4, "Dark2")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 10, family = plot_font), axis.title = element_text(size = 12, family = plot_font),
        legend.text = element_text(size = 12, family = plot_font), legend.title = element_text(vjust = 0.1, size = 12, family = plot_font), legend.position = "bottom")

fig_effects <- ggarrange(climate_plots, year_re, nrow = 2, heights = c(1.4,1), labels = c("", "e"),font.label = list(family = plot_font, size = 18))
ggsave(here(output_folder, "fig_3.png"), fig_effects,  width = 300, height = 300, units = "mm", bg="white", dpi = 300)

# Supplementary Figure - comparing spatial differences
# Change in spatial effect - by outbreak or non-outbreak year

difference_outbreak <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 1) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Full model`) - abs(`Baseline`)) |> # looking at the absolute change in distance from 0 
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Outbreak years") +
  #scale_fill_gradientn(colours = met.brewer("Cassatt2", 9)) +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

difference_non_outbreak <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 0) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Full model`) - abs(`Baseline`)) |>  # looking at the absolute change in distance from 0
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Non-outbreak years") +
  #scale_fill_gradientn(colours = met.brewer("Cassatt2", 9)) +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

difference_outbreak_climate <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 1) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Climate-only`) - abs(`Baseline`)) |>  # looking at the absolute change in distance from 0 
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Outbreak years") +
  #scale_fill_gradientn(colours = met.brewer("Cassatt2", 9)) +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  theme_void() +
  guides(fill = "none") +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

difference_non_outbreak_climate <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 0) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Climate-only`) - abs(`Baseline`)) |>  # looking at the absolute change in distance from 0
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Non-outbreak years") +
  guides(fill = "none") +
  #scale_fill_gradientn(colours = met.brewer("Cassatt2", 9)) +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

difference_outbreak_mechanism <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 1) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Epidemic-only`) - abs(`Baseline`)) |> # looking at the absolute change in distance from 0 
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Outbreak years") +
  #scale_fill_gradientn(colours = met.brewer("Cassatt2", 9)) +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  guides(fill = "none") +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

difference_non_outbreak_mechanism <- spatial_re |> 
  group_by(adm2_name, component, model, outbreak_year) |> 
  summarise(mean = mean(median)) |> 
  ungroup() |> 
  filter(component == "uv_joint" & outbreak_year == 0) |> 
  pivot_wider(values_from = "mean", names_from = "model") |> 
  mutate(diff = abs(`Epidemic-only`) - abs(`Baseline`)) |>  # looking at the absolute change in distance from 0
  left_join(shp_adm2, by = "adm2_name") |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill  = diff)) + 
  labs(fill = "  Difference in \nmarginal effect", x = "Non-outbreak years") +
  scale_fill_gradient2(low = "#40004b", mid = "white", high = "#00441b", limits = c(-1.4, 1.4)) +
  guides(fill = "none") +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.text = element_text(size = 10, family = plot_font),
    legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"), legend.position = "bottom", axis.title.x = element_text(size = 16, family = plot_font))

row1 <- ggarrange(difference_outbreak_climate, difference_non_outbreak_climate, common.legend = TRUE, legend = "bottom", nrow =1, ncol = 2, labels = c("a", "b"),font.label = list(family = plot_font, size = 18))
row1 <- annotate_figure(row1, left = text_grob("Climate-only", color = "black", rot = 90, size = 20, family = plot_font))
row2 <- ggarrange(difference_outbreak_mechanism, difference_non_outbreak_mechanism, common.legend = TRUE, legend = "bottom", nrow =1, ncol = 2, labels = c( "c", "d"),font.label = list(family = plot_font, size = 18))
row2 <- annotate_figure(row2, left = text_grob("Epidemic-only", color = "black", rot = 90, size = 20, family = plot_font))
row3 <- ggarrange(difference_outbreak, difference_non_outbreak, common.legend = TRUE, legend = "bottom", nrow =1, ncol = 2, labels = c("e", "f"),font.label = list(family = plot_font, size = 18))
row3 <- annotate_figure(row3, left = text_grob("Full model", color = "black", rot = 90, size = 20, family = plot_font))

spatial_compare_plot <- ggarrange(row1, row2, row3, nrow = 3, common.legend = TRUE, heights = c(1,1,1.1))
ggsave(here(output_folder, "sfig_5.png"), spatial_compare_plot,  width = 350, height = 350, units = "mm", bg="white", dpi = 300)
