# -----------------------------------------------------------------------------------------------------------
# Run k-fold spatial cross-validation model evaluation
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))

c_args = commandArgs(trailingOnly = TRUE)

forms = c_args[[1]]
mod_name = c_args[[2]]
n_reps = as.numeric(c_args[[3]])
k_folds = as.numeric(c_args[[4]])

cat("K-fold evaluation for: ", mod_name, ". Running ", n_reps, "reps for ", k_folds, "folds.")

start_date <- as.character(Sys.Date())

sessionInfo()

source(here("code", "utils_fn.R"))
source(here("code", "evaluate-kfold_fn.R"))
source(here("code", "create-lagged-data_fn.R"))

adm2_data <- readRDS(here("data", "adm2-climate-epi.rds"))
df_model <- lag_data(adm2_data)
nb_graph <- inla.read.graph(filename = here("data", "adm2_nb-graph.map"))

output_path <- here("output", "kfold", start_date, mod_name)

ifelse(!dir.exists(output_path), dir.create(output_path, recursive = TRUE), FALSE)

# Prepare data frame ----------------------------------------------------------------------------------------

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

## Use a seasonal 75th percentile by adm2
## For a given month and adm we define an outbreak threshold of the 75th percentile of weekly cases in that month
## using all years (note different from thresholds for tscv)

df_model <-df_model |> 
  group_by(month, adm2_id) |> 
  mutate(threshold = quantile(reported, 0.75, na.rm = TRUE)) |>
  ungroup()

# Here we define an outbreak week where the number of cases is > seasonal 75th percentile
# threshold, and an outbreak year as having more than 12 outbreak weeks

df_model <- df_model |>
  mutate(outbreak_week = case_when(reported > threshold ~ 1, TRUE ~ 0)) 

# Run spatial k-fold evaluation -------------------------------------------------------------------------------------------------------------------

kfold_out <- list()

for(nn in 1:n_reps){
  adm2_folds <- df_model |> select(adm2_id) |> distinct()
  adm2_folds$kfold <- kfold_func(adm2_folds, k = 5) # allocate each adm2 a fold
  write.csv(adm2_folds, here(output_path, paste0("kfold-evaluation_adm2-folds_", mod_name, "_rep-", nn, ".csv")))

  df_eval <- df_model |>
    left_join(adm2_folds, by = "adm2_id")

evaluate_kfold(data_input = df_eval,
                              nb_graph = nb_graph,
                              forms = forms,
                              mod_name = mod_name,
                              k_folds = k_folds,
                              rep = nn,
                              output_path = output_path)


}

# Read in each rep and save together 

for(nr in 1:n_reps){ 
kfold_out[[nr]] <- qread(here(output_path, paste0("kfold-output",  "_",mod_name, "-rep-",nr, ".qs")))
}

qsave(kfold_out, here(output_path, paste0("kfold-output",  "_",mod_name, ".qs")))

# Process output and plot predictions ---------------------------------------------------------------------------------------------------


adm2_preds <- get_adm2_preds(kfold_out, df_model, mod_name)
cat("Admin 2 predictions generated.")
print(head(adm2_preds))
adm1_preds <- get_adm1_preds(kfold_out, df_model, mod_name)
cat("Admin 1 predictions generated.")
print(head(adm1_preds))
rm(kfold_out)

preds_out <- list(adm2_preds = adm2_preds, adm1_preds = adm1_preds)
qsave(preds_out, here(output_path, paste0("kfold-predictions",  "_",mod_name,".qs")))
cat("Admin predictions saved.")

adm2_pred_plot <- 
  ggplot(adm2_preds) + 
  geom_line(aes(y = reported, x = date), col = "#000000", size = 0.3) +
  geom_line(aes(x = date, y = median), col = "#7570B3") +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975), alpha = 0.3, fill = "#7570B3") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ adm2_name, scale = "free_y") +
  theme_classic2() +
  theme(
    axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14, color = "white"),
    legend.title = element_text(size = 14), legend.text = element_text(size = 14),
    legend.key.height = unit(1.2, "cm"), legend.position = "none")

adm1_pred_plot <- 
  ggplot(adm1_preds) + 
  geom_line(aes(y = reported, x = date), col = "#000000", size = 0.3) +
  geom_line(aes(x = date, y = median), col = "#7570B3") +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975), alpha = 0.3, fill = "#7570B3") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ adm1_name, scale = "free_y") +
  theme_classic2() +
  theme(
    axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14, color = "white"),
    legend.title = element_text(size = 14), legend.text = element_text(size = 14),
    legend.key.height = unit(1.2, "cm"), legend.position = "none")

ggsave(here(output_path, paste0("kfold-adm2-predictions-plot",  "_",mod_name, "_", ".jpeg")), adm2_pred_plot, width = 36, height = 24)
ggsave(here(output_path, paste0("kfold-adm1-predictions-plot",  "_",mod_name, "_", ".jpeg")), adm1_pred_plot, width = 36, height = 24)


# Score predictions -------------------------------------------------------------------------------

cat("Scoring predictions.")
score_predictions(preds_input = adm2_preds, 
                  data_input = df_model,
                  output_path = output_path, 
                  mod_name = mod_name,
                  adm_level = "adm2_name")
cat("Finished scoring predictions.")


