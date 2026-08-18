# -----------------------------------------------------------------------------------------------------------
# Estimate force of infection (FOI) from serosurvey data
# Note that serology data is not included in this repository
# -----------------------------------------------------------------------------------------------------------

source(here("code", "helper-rsero_fn.R"))

scores <- NULL

# Load in serology data (not included)

# Wrangle to correct format
sero_dom_rep <- sero_elisa  |>  
  select(age, age_cat, age_cat_mid, den_cos_2022_mix1_4_ig_g_elisa_p_n, province) |> 
  rename(outcome = den_cos_2022_mix1_4_ig_g_elisa_p_n) |> 
  mutate(outcome = case_when(outcome == "Positive" ~ 1,
                             outcome == "Negative" | outcome == "Equivocal" | outcome == "QNS" | outcome == "*" ~ 0,
                             T ~ NA_integer_))

denv_rsero <- SeroData(age_at_sampling = sero_dom_rep$age, Y = sero_dom_rep$outcome, sampling_year = 2021)

# Test different model structures

# 1. Constant FOI -----------------------------------------------------------------------------------

mod_constant <- FOImodel(type = "constant")
fit_constant <- fit(data = denv_rsero, model = mod_constant, chains = 4, iter = 20000)
saveRDS(trim_serofit(fit_constant), here("output", "rsero_fit-constant.rds"))

## Check fit

plot_constant <- plot_fit(fit_constant, age_class = 5)
ggsave(filename = here("figures", "rsero_fit-constant.png"), plot_constant[[1]],width = 14, height = 7)

## Check convergence
rstan::traceplot(fit_constant$fit, pars = c("annual_foi", "rho")) # Seems fine
print(fit_constant$fit)

## Check estimated FOI
constant_foi <- plot(fit_constant, YLIM = 1) # Estimated FOI of 0.1, approx 12.2% attack rate per year
ggsave(filename = here("figures", "rsero_foi-constant.png"), constant_foi[[1]],width = 14, height = 7)
foi <- extract_fois(fit_constant)
attack_rates <- get_attack_rate(foi) 

# Plot posteriors
plot_posterior(fit_constant)

# Calculate WAIC
scores <- rbind(scores, compute_dic_waic(fit_constant)) |> 
  as.data.frame()
scores$mod <- "constant"

# 2. Constant FOI with seroreversion -----------------------------------------------------------------------------------

mod_constant_serorev <- FOImodel(type = "constant", seroreversion = 1,
                                 priorRho = 0.1)
fit_constant_serorev <- fit(data = denv_rsero, model = mod_constant_serorev, chains = 4, iter = 20000)
saveRDS(trim_serofit(fit_constant_serorev), here("output", "rsero_fit-constant_serorev.rds"))

## Check fit

plot_constant_serorev <- plot_fit(fit_constant_serorev, age_class = 5)
ggsave(filename = here("figures", "rsero_fit-constant_serorev.png"), plot_constant_serorev[[1]],width = 14, height = 7)

## Check convergence
rstan::traceplot(fit_constant_serorev$fit, pars = c("annual_foi", "rho")) # Seems fine
print(fit_constant_serorev$fit)

## Check estimated FOI
constant_serorev_foi <- plot(fit_constant_serorev, YLIM = 1) # Estimated FOI of 0.1, approx 12.2% attack rate per year
ggsave(filename = here("figures", "rsero_foi-constant_serorev.png"), constant_serorev_foi[[1]],width = 14, height = 7)
foi <- extract_fois(fit_constant_serorev)
attack_rates <- get_attack_rate(foi) 

params <- parameters_credible_intervals(fit_constant_serorev) |> 
  rownames_to_column() |> 
  rename(var = rowname)

# Calculate WAIC
scores <- rbind(scores, c(compute_dic_waic(fit_constant_serorev), c("mod" = "constant_serorev"))) 

