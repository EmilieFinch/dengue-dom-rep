# -----------------------------------------------------------------------------------------------------------
# Function to run INLA models
#' @param forms character string with INLA formula for covariates
#' @param df_model data frame for model selection, must include variables specified in forms
#' @param res string of baseline random effects to include
#' @param nb_graph neighbourhood graph for INLA bym2 random effects
#' @param models_include logical, return INLA model output
#' @return List including:
#'   adeq_stats: data frame of adequacy statistics
#'   fitted_vals: data frame of fitted values
#'   fixed_effects: data frame of fixed effects
#'   random_effects: list of random effects
# -----------------------------------------------------------------------------------------------------------

fit_inla <- function(forms) { 
  
  # Model set up ---------------------------------------------------------------------------------------------
  
  models_include = FALSE
  
  tmpdr <- withr::local_tempdir()
  inla.setOption(working.directory = normalizePath(tmpdr))
  
  # Load graph
  # Set up random effects
  
  eweek_re <- " + f(epiweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
  adm2_year_re <- "+f(adm2_id, model = 'bym2', graph = nb_graph, scale.model = TRUE, hyper = precision_prior, replicate = year_index)"
  
  res <- paste0(eweek_re, adm2_year_re)
  
  # Define objects for model output
  
  fits <- NULL
  params <- NULL
  random <- NULL
  adeq_stats <- NULL
  mods <- NULL
  
  # Set up prior for random effects
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  # Run  model 
  
  for(form in 1:length(forms)){ 
    
    start_time <- Sys.time()
    
    
    mod <- inla(as.formula(paste0("reported~1", res, forms[form])),
                family = "nbinomial",
                control.inla = list(strategy = "adaptive"),
                control.predictor = list(link = 1, compute = TRUE),
                control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
                control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
                num.threads = 40,
                verbose = FALSE,
                data = df_model
    )
    
    #fits[[form]] <- mod$summary.fitted.values 
    #fits[[form]] <- cbind(fits, adm2_name = df_model$adm2_name, date = df_model$date, observed = df_model$reported)
    params[[form]] <- w_params(mod$summary.fixed)
    random[[form]] <- mod$summary.random
    mods[[form]] <- mod
    
    # Extract model adequacy criteria and save it in the table. Both MAE and Rsq are calculated relative to the base model (REs only)
    # MAE
    mean_fit <- hydroGOF::mae(mod$summary.fitted.values$mean, df_model$reported, na.rm = T)
    median_fit <- hydroGOF::mae(mod$summary.fitted.values$`0.5quant`, df_model$reported, na.rm = T)
    
    crps_mod <- crps(mod, df_model)
    print(paste0("CRPS calculated."))
    
    # Create adequacy table
    add <- data.frame(
      form = forms[form],
      dic = mod$dic$dic,
      deviance_mean = mod$dic$deviance.mean,
      waic = mod$waic$waic,
      cpo = -mean(log(mod$cpo$cpo)),
      mean_mae = mean_fit,
      median_mae = median_fit,
      crps = crps_mod
    )
    
    adeq_stats <- bind_rows(adeq_stats, add)
    
    qsave(list(adeq_stats = adeq_stats, fixed_effects = params, random_effects = random), here("output", paste0("model-selection_temp.qs")))
    
    
    print(Sys.time() - start_time)
    print(forms[form])
    print(paste0("Covariate model completed for ", form, "/", length(forms), "."))
  }
  
  if(models_include == TRUE){ 
    return(list(
      adeq_stats = adeq_stats,
     # fitted_vals = fits,
      fixed_effects = params,
      random_effects = random,
      mod = mods
    )) } else {
      return(list(
        adeq_stats = adeq_stats,
        fitted_vals = fits,
        fixed_effects = params,
        random_effects = random
      ))
    }
}


crps <- function(crps_mod, crps_data) {
  # Sample from the posterior
  s <- 1000 # Number of samples
  xx <- inla.posterior.sample(1000, crps_mod)
  # Extract values of interest
  xx_s <- inla.posterior.sample.eval(
    function(...) {
      c(
        theta[1], # This is the size parameter of the negative binomial distribution (overdispersion)
        Predictor
      )
    },
    xx
  )
  # Create posterior predictive sample
  y_pred <- matrix(NA, nrow(crps_data), s)
  for (s_idx in 1:s) {
    xx_sample <- xx_s[, s_idx]
    y_pred[, s_idx] <- rnbinom(nrow(crps_data),
                               mu = exp(xx_sample[-1]), # Predicted means
                               size = xx_sample[1]
    ) # Overdispersion parameter
  }
  
  # Calculate CRPS
  crps_val <- median(scoringutils::crps_sample(crps_data$reported, y_pred))
  
  return(crps_val)
}

process_output <- function(model_output){
  adeq_stats <- NULL
  fitted_vals <- list()
  fixed_effects <- list()
  random_effects <- list()
  for(ml in 1:length(model_output)){
    if(is.list(model_output[[ml]])){ 
    adeq_stats <- rbind(adeq_stats, model_output[[ml]]$adeq_stats)
    fitted_vals[[ml]] <- model_output[[ml]]$fitted_vals
    fixed_effects[[ml]] <- model_output[[ml]]$fixed_effects
    random_effects[[ml]] <- model_output[[ml]]$random_effects
    }
  }
  return(list(adeq_stats = adeq_stats, fitted_vals = fitted_vals, fixed_effects = fixed_effects, random_effects = random_effects))
}

rsq <- function(processed_output){
  rsqs <- NULL
  for(mds in 1:nrow(processed_output$adeq_stats)){
    dev <- processed_output$adeq_stats[mds,]$deviance_mean
    nulldev <- processed_output$adeq_stats[processed_output$adeq_stats$form == "",]$deviance_mean
    n <- nrow(processed_output$fitted_vals[[mds]])
    x <- round(1 - exp((-2 / n) * ((dev / -2) - (nulldev / -2))), 3)
    rsqs <- c(rsqs, x)
  }
  return(rsqs)
}


# Helper function to wrangle model parameters

w_params <- function(x) {
  x %>%
    rownames_to_column(var = "var")
}

# Adapted from Gibb et al
wrangle_random <- function(summary_random, effect_name, is_bym = FALSE, transform = FALSE){
  # extract model effect
  rf = summary_random %>%
    dplyr::rename("value"=1, "lower"=4, "median"=5, "upper"=6)
  
  # label by grouping factor (if not replicated, group is 1 for all observations)
  rf$group = rep(1:as.vector(table(rf$value)[1]), each=n_distinct(rf$value))
  
  # partition BYM into u and v components
  if(is_bym){
    rf$component = rep(c("uv_joint", "u_besag"), each=n_distinct(rf$value)/2)
    rf$value = rep(1:(n_distinct(rf$value)/2), n_distinct(rf$group)*2)
  }
  
  # back transform if specified
  if(transform == TRUE){
    rf[ , 2:7 ] = exp(rf[ , 2:7])
  }
  
  # name and return
  rf$effect = effect_name
  return(rf)
  
}