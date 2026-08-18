# -----------------------------------------------------------------------------------------------------------
# Generate catalytic model figures
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))
source(here("code", "helper-rsero_fn.R"))

# Load model fits and compute information criteria
fit_constant <- readRDS(here("output", "rsero_fit-constant.rds"))
constant_criteria <- compute_dic_waic(fit_constant)

fit_constant_serorev <- readRDS(here("output", "rsero_fit-constant_serorev.rds"))
serorev_criteria <- compute_dic_waic(fit_constant_serorev)

criteria_table <-data.frame(model = c("constant", "serorev"), 
           waic = c(round(constant_criteria$WAIC,2), round(serorev_criteria$WAIC,2)),
           dic = c(round(constant_criteria$DIC,2), round(serorev_criteria$DIC,2)))

write.csv(criteria_table, here(output_folder, "catalytic-criteria.csv"))

# Figure 3 -------------------------------------------------------------------------------------------------------------
## Panel b - model fit to data

model_fit_plot <- plot_fit(fit_constant, age_class = 5)[[1]]
ggsave(here(output_folder, "figure_2.png"), model_fit_plot,  width = 250, height = 200, units = "mm", bg="white", dpi = 300)

# Supplementary figure 4
## Trace plot
trace_plot <- rstan::traceplot(fit_constant$fit, pars = c("annual_foi")) +
  labs(col = "Chain") +
  theme(
    legend.position = "bottom",
    #  legend.spacing.x = unit(0, "mm"),
    legend.background = element_rect(fill = "white", colour = "white"),
    legend.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"),
    legend.text = element_text(family = plot_font, size = 12), 
    legend.title = element_text(family = plot_font, size = 12),
    axis.title = element_text(size = 12, family = plot_font), axis.text.y = element_text(size = 12, family = plot_font),
    axis.text.x = element_text(size = 12, family = plot_font), strip.text.x = element_text(size = 12, family = plot_font)) 

ggsave(here(output_folder, "sfig_4.png"), trace_plot,  width = 400, height = 200, units = "mm", bg="white", dpi = 300)
