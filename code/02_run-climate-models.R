# -----------------------------------------------------------------------------------------------------------
# Run INLA space-time climate models
# -----------------------------------------------------------------------------------------------------------

source(here("code", "00_load-packages.R"))

# Load data ------------------------------------------------------------------------------------------------

source(here("code", "inla-loop_fn.R"))
source(here("code", "utils_fn.R"))
source(here("code", "create-lagged-data_fn.R"))

# Load shape files for adjacency matrix
shp_adm2 <- readRDS(here("data", "adm2-shp-clean.rds"))

# Load and lag data
adm2_data <- readRDS(here("data", "adm2-climate-epi.rds"))
df_model <- lag_data(adm2_data)


# Model set up ---------------------------------------------------------------------------------------------

# Load graph
if(!file.exists(here("data", "adm2_nb-graph.map"))){ 
  nb <- poly2nb(shp_adm2) # compute the adjacency matrix
  nb2INLA(here("data", "adm2_nb-graph.map"), nb)
}

nb_graph <- inla.read.graph(filename = here("data", "adm2_nb-graph.map"))

# Data set up

# Add ids for random effects
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

# Set up formulae for modelling ------------------------------------------------------------------------------------

climate_form <- "+f(inla.group(nino34_12_wk_avg_6, n = 9),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(tasmax_scale_8_wk_avg_4),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(spi_53),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(hurs_scale_10_wk_avg_2),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
foi_vars <- c("+log_weight_cases", "+offset(log_susceptibility)", "+log_weight_cases+offset(log_susceptibility)")
foi_forms <- expand.grid(foi_vars, climate_form) |> 
  mutate(forms = paste0(Var1, Var2)) |> 
  pull(forms)

forms <- c("", climate_form, foi_forms, paste0("+log_weight_cases+offset(log_susceptibility)"))

# Fit inla models -------------------------------------------------------------------------------------------------

models_out <- fit_inla(forms)

qsave(models_out, here("output", paste0("selected-models-output.qs")))
