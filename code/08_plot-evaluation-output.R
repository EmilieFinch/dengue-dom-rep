# -----------------------------------------------------------------------------------------------------------
# Visualise model evaluation results (TSCV and k-fold)
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))
source(here("code", "utils_fn.R"))

grid <- read_xlsx(here("data", "adm1-grid.xlsx"))

# Model output to plot
model_name <- "full-model"


# Figure 4 - TSCV predictions for full model ---------------------------------------------------------------

tscv_predictions <- qread(here("output", "tscv", model_name, paste0("tscv-predictions_", model_name,".qs")))

fig_4 <- ggplot(tscv_predictions$adm1_preds) + 
  geom_line(aes(y = reported, x = date), col = "#737373", size = 0.15) +
  geom_line(aes(x = date, y = as.integer(median)), col = "#016c59", size = 0.3) +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975), alpha = 0.3, fill = "#016c59") +
  labs(x = NULL, y = NULL) +
  scale_x_date(limits = c(as.Date("2015-04-05"), as.Date("2023-03-26"))) +
  facet_geo(~ adm1_name, grid = grid, scale = "free_y") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_blank(), legend.title = element_text(size = 15, family = plot_font), legend.position = "right",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

ggsave(here(output_folder, "fig_4.png"), fig_4,  width = 450, height = 300, units = "mm", bg="white", dpi = 300)

# Supplementary Figure 8 -  K-fold predictions for full model

kfold_predictions <- qread(here("output", "kfold", model_name, paste0("kfold-predictions_", model_name,".qs")))

supp_fig_kfold <- ggplot(kfold_predictions$adm1_preds) + 
  geom_line(aes(y = reported, x = date), col = "#737373", size = 0.15) +
  geom_line(aes(x = date, y = as.integer(median)), col = "#016c59", size = 0.3) +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975), alpha = 0.3, fill = "#016c59") +
  labs(x = NULL, y = NULL) +
  scale_x_date(limits = c(as.Date("2015-04-05"), as.Date("2023-03-26"))) +
  facet_geo(~ adm1_name, grid = grid, scale = "free_y") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, family = plot_font, face = "bold"), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_blank(), legend.title = element_text(size = 15, family = plot_font), legend.position = "right",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

ggsave(here(output_folder, "sfig_8.png"), supp_fig_kfold,  width = 450, height = 300, units = "mm", bg="white", dpi = 300)


# Figure 5 - Influence of covariates of covariates on temporal and spatial predictions of dengue incidence
mod_names <- c("full-model", "without-nino", "without-spi", "without-tas", "without-hurs", "without-sus", "without-foi",  "climate-only", "mechanistic-only", "baseline")

tscv_crps <- data.frame()
for(mn in mod_names){
  tscv_crps <- tscv_crps |> 
    bind_rows(qread(here("output", "tscv",  mn, paste0("tscv-output_", mn, ".qs")))$crps |> 
                mutate(model = mn))

}

if(!file.exists(here("output", "kfold",  "crps-values.csv"))){
  kfold_crps <- data.frame()
  for(mn in mod_names){
    print(mn)
    for(rep in 1:10){
      print(rep)
      kfold_crps <- kfold_crps |> 
        rbind(qread(here("output", "kfold",  mn, paste0("kfold-output_", mn, ".qs")))[[rep]]$crps |> 
                mutate(model = mn, rep = rep))
    }
  }
  write.csv(kfold_crps, file = here("output", "kfold",  "crps-values.csv"))
  
} else {
  kfold_crps <- read.csv(here("output", "kfold",  "crps-values.csv"))
}

tscv_influence <- tscv_crps |> 
  mutate(full_crps = tscv_crps |>  filter(model == 'full-model') |> dplyr::pull(crps_natural)) |> 
  mutate(full_crps_log = tscv_crps |>  filter(model == 'full-model') |> dplyr::pull(crps_log)) |> 
  mutate(crps_natural_difference = (crps_natural - full_crps)/full_crps*100, crps_log_difference = (crps_log - full_crps_log)/full_crps_log*100) |> 
  filter(model!= "full-model" & model != "climate-only" & model != "mechanistic-only") |>  
  mutate(class = case_when(model == "without-sus" | model == "without-foi" ~ "Epidemic",
                           model == "baseline" ~ "Baseline",
                           T ~ "Climate")) |> 
  mutate(model = factor(model, levels = c("without-hurs","without-tas","without-spi","without-nino","without-sus","without-foi","baseline"), labels = c("Relative humidity", "Maximum temperature", "SPI-12", "Niño 3.4", "Susceptibility", "Weighted cases", "Baseline (all)"))) |>
  ggplot() +
  ggtitle("TSCV (temporal)") +
  geom_col(aes(x = model, y = crps_natural_difference, fill = class), width = 0.5) + 
  labs(y = "% increase in CRPS relative to full model", x = "Covariate excluded", fill = NULL) +
  scale_fill_manual(values = brewer.pal(3, "Dark2")) +
  geom_hline(aes(yintercept = 0), lty = "longdash", size = 0.5) +
  scale_y_continuous(limits = c(-0.8, 15.1)) +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 12, family = plot_font, hjust = 0.5), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_text(size = 12, family = plot_font), legend.title = element_text(size = 15, family = plot_font), legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))


kfold_influence <- kfold_crps |> 
  inner_join(kfold_crps |>  filter(model == 'full-model') |> rename(crps_full = crps_natural, crps_log_full = crps_log) |> dplyr::select(-model), by = "rep") |> 
  mutate(crps_difference = (crps_natural - crps_full)/crps_full*100) |> 
  group_by(model) |> 
  summarise(median = median(crps_difference), q_0.025 = quantile(crps_difference, 0.025), q_0.975 = quantile(crps_difference, 0.975)) |> 
  filter(model!= "full-model" & model != "climate-only" & model != "mechanistic-only") |>  
  mutate(class = case_when(model == "without-sus" | model == "without-foi" ~ "Epidemic",
                           model == "baseline" ~ "Baseline",
                           T ~ "Climate")) |> 
  mutate(model = factor(model, levels = c("without-hurs","without-tas","without-spi","without-nino","without-sus","without-foi","baseline"), labels = c("Relative humidity", "Maximum temperature", "SPI-12", "Niño 3.4", "Susceptibility", "Weighted cases", "Baseline (all)"))) |> 
  ggplot() +
  ggtitle("K-fold (spatial)") +
  geom_point(aes(x = model, y = median, col = class)) + 
  geom_linerange(aes(x = model, ymin = q_0.025, ymax = q_0.975, col = class), size = 0.4) +
  geom_hline(aes(yintercept = 0), lty = "longdash", size = 0.5) +
  labs(y = "% increase in CRPS relative to full model", x = NULL, color = NULL) +
  scale_color_manual(values = brewer.pal(3, "Dark2")) +
  scale_y_continuous(limits = c(-0.8, 15.1)) +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 12, family = plot_font, hjust = 0.5), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_text(size = 12, family = plot_font), legend.title = element_text(size = 15, family = plot_font), legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

fig_5 <- ggarrange(tscv_influence, kfold_influence, common.legend = TRUE, legend = "bottom", nrow = 1, labels = c("a", "b"), font.label = list(family = plot_font, size = 15))
ggsave(here(output_folder, "fig_5.png"), fig_5,  width = 300, height = 150, units = "mm", bg="white", dpi = 300)

# Supplementary Table of CRPS values

kfold_table <- kfold_crps |> 
  group_by(rep) |> 
  mutate(max_crps = max(crps_natural), crpss = (1 - crps_natural/max_crps)*100) |> 
  ungroup() |> 
  group_by(model) |> 
  summarise(median = median(crps_natural), q_0.025 = quantile(crps_natural, 0.025), q_0.975 = quantile(crps_natural, 0.975), median_crpss = median(crpss), q_0.025_crpss = quantile(crpss, 0.025), q_0.975_crpss = quantile(crpss, 0.975)) |> 
  mutate(crps = paste0(round(median,3), " 95% CrI: ", round(q_0.025,3), "-", round(q_0.975,3)), crpss = paste0(round(median_crpss,3), " 95% CrI: ", round(q_0.025_crpss,3), "-", round(q_0.975_crpss,3))) |> 
  dplyr::select(model, crps, crpss) |> 
  mutate(model = factor(model, levels = c("full-model", "climate-only", "mechanistic-only", "without-nino", "without-spi", "without-tas", "without-hurs", "without-sus", "without-foi", "baseline"))) |> 
  arrange(model)

write.csv(kfold_table, here(output_folder, "kfold-crps_table.csv"))        

tscv_table <- tscv_crps |> 
  mutate(crps = round(crps_natural, 3)) |> 
  dplyr::select(model, crps) |> 
  mutate(crpss = round((1 - crps/max(crps))*100, 1)) |> 
  mutate(model = factor(model, levels = c("full-model", "climate-only", "mechanistic-only", "without-nino", "without-spi", "without-tas", "without-hurs", "without-sus", "without-foi", "baseline"))) |> 
  arrange(model)

write.csv(tscv_table, here(output_folder, "tscv-crps_table.csv"))        

# Scoring plot
tscv_scoring <- tscv_predictions$adm2_preds |> 
  rbind(qread(here("output", "tscv",  "baseline", paste0("tscv-predictions_", "baseline",".qs")))$adm2_preds)


tscv_scoring <- tscv_scoring |> 
  rename(q_0.5 = median, model = mod_name) |> 
  mutate(target_end_date = date, true_value = reported) |>
  pivot_longer(cols = starts_with("q_"), names_to = "quantile", values_to = "prediction") |> 
  mutate(quantile = as.numeric(str_remove(quantile, "q_")))

# add_coverage() was removed from scoringutils - get_coverage() is its v2 replacement,
# but it works on the forecast object rather than the scored output, so it's computed
# separately here and joined back in rather than chained through score()

tscv_forecast <- tscv_scoring |>
  as_forecast_quantile(observed = "true_value", predicted = "prediction", quantile_level = "quantile") |>
  transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
  filter(!is.na(predicted))

get_coverage_wide <- function(forecast, by) {
  get_coverage(forecast, by = by) |>
    filter(interval_range %in% c(50, 95)) |>
    distinct(across(all_of(by)), interval_range, interval_coverage) |>
    pivot_wider(names_from = interval_range, values_from = interval_coverage, names_prefix = "coverage_")
}

scores <- tscv_forecast |>
  score() |>
  summarise_scores(by = c("scale", "model"), na.rm = TRUE) |>
  left_join(get_coverage_wide(tscv_forecast, by = c("scale", "model")), by = c("scale", "model")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

scores_week <- tscv_forecast |>
  score() |>
  summarise_scores(by = c("scale", "date", "model"), na.rm = TRUE) |>
  left_join(get_coverage_wide(tscv_forecast, by = c("scale", "date", "model")), by = c("scale", "date", "model")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# Plots over time

## Coverage plot
coverage_plot <- scores_week |> 
  filter(scale == "natural") |> 
  rename("50%" = coverage_50, "95%" = coverage_95) |> 
  pivot_longer(cols = ends_with("%"), names_to = "coverage", values_to = "value") |> 
  ggplot() +
  geom_line(aes(x = date, y = value, linetype = coverage, col = model)) + 
  geom_hline(aes(yintercept = 0.5), lty = "longdash") +
  geom_hline(aes(yintercept = 0.95), lty = "longdash") +
  scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
  scale_color_manual(values = c("#1B9E77", "#7570B3"), labels = c("Baseline", "Full model")) +
  labs(x = "Date", y = "Interval coverage", col = NULL, linetype = "Coverage level") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, family = plot_font, hjust = 0.5), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_text(size = 12, family = plot_font), legend.title = element_text(size = 12, family = plot_font), legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

## Bias plot
bias_plot <- scores_week |> 
  filter(scale == "natural") |> 
  ggplot() +
  geom_line(aes(x = date, y = bias, col = model)) + 
  geom_hline(aes(yintercept = 0), lty = "longdash") +
  scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
  scale_color_manual(values = c("#1B9E77", "#7570B3"), labels = c("Baseline", "Full model")) +
  labs(x = "Date", y = "Bias", col = NULL) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, family = plot_font, hjust = 0.5), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_text(size = 12, family = plot_font), legend.title = element_text(size = 15, family = plot_font), legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

## WIS plot
wis_plot <- scores_week |> 
  filter(scale == "natural") |> 
  ggplot() +
  geom_line(aes(x = date, y = wis, col = model)) +
  scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
  scale_color_manual(values = c("#1B9E77", "#7570B3"), labels = c("Baseline", "Full model")) +
  labs(x = "Date", y = "Weighted interval score", col = NULL) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, family = plot_font, hjust = 0.5), axis.text.x = element_text(size = 10, family = plot_font), axis.text.y = element_text(size = 11, family = plot_font), 
        axis.title=element_text(size = 12, family = plot_font), legend.title = element_text(size = 15, family = plot_font), legend.position = "bottom",
        legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 12, family = plot_font), legend.spacing.x = unit(0.8, "cm"), strip.text.x = element_text(size = 8, family = plot_font, margin = margin(0.2,0.2,0.2,0.2, "cm")))

score_plots <- ggarrange(coverage_plot, bias_plot, wis_plot, nrow = 3, heights = c(1.2, 1, 1, 1), common.legend = TRUE, legend = "bottom", labels = c("a", "b", "c"), font.label = list(family = plot_font, size = 15))
ggsave(here(output_folder, "sfig_7.png"), score_plots,  width = 360, height = 240, units = "mm", bg="white", dpi = 300)

