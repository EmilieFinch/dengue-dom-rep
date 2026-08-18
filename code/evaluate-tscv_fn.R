# -----------------------------------------------------------------------------------------------------------
# Function to run time-series cross validation with INLA
#' @param data_input data frame for modelling
#' @param nb_graph neighbourhood graph for spatial random effects
#' @param forms character string with INLA formula for covariates
#' @param mod_name character string for short hand model name
#' @param output_path character string of output file path
# -----------------------------------------------------------------------------------------------------------

evaluate_tscv <- function(data_input,
                          nb_graph,
                          forms,
                          mod_name,
                          output_path){
  
  # Define objects for model output
  posterior_samples <- NULL
  outbreak_probs <- NULL
  
  
  # Set up parameters
  s <- 1000
  pred_full <- array(dim = c(522, 155, 1000)) # Dimensions are: week, adm2_id, samples
  tscv_out <- list()
  
  # Set up random effects
  eweek_re <- " + f(epiweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
  adm2_year_re <- "+f(adm2_id, model = 'bym2', graph = nb_graph, scale.model = TRUE, hyper = precision_prior, replicate = year_index)"
  res <- paste0(eweek_re, adm2_year_re) # What to do for baseline model - seasonal and spatial only?
  
  # Set up prior for random effects
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  cat(paste0("Running TSCV evaluation with formula:", forms, res, "\n"))
  
  start_week <- data_input |> 
    arrange(date, adm2_name) |> # check arranged correctly
    filter(year_index == 3) |> # Use first 2 years to train 
    slice(1) |> 
    pull(date_index)
  
  for (wk in start_week:max(data_input$date_index)){ # 230
    cat(paste0("Generating predictions for week: ", wk))
    
    df_pred <- data_input |> 
      filter(date_index <= wk) |> 
      mutate(reported = case_when(date_index == wk ~ NA_integer_, TRUE ~ reported)) # set cases to NA for prediction week
    
    if(wk != 372){ 
    cat("Starting model fitting.")
    mod  <- inla(as.formula(paste0("reported~1", res, forms)),
                       family = "nbinomial",
                       #offset = log(population / 100000),
                       control.inla = list(strategy = "adaptive"),
                       control.predictor = list(link = 1, compute = TRUE),
                       control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
                       control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
                       num.threads = 10,
                       verbose = FALSE,
                       data = df_pred)
    
    cat("Model fitting completed.")
    
    inla_samples <- inla.posterior.sample(s, mod)
    rm(mod)
    cat("Samples generated from posterior.") 
    
    # Get idx for each adm2 in week of interest
    idx <- df_pred |> 
      mutate(row_number = row_number()) |> 
      filter(date_index == wk) |> 
      pull(row_number)
    
     samples <- inla.posterior.sample.eval(
      function(...) {
        c(
          theta[1], # This is the size parameter of the negative binomial distribution (overdispersion)
          Predictor
        )
      },
      inla_samples
    )
     cat("Necessary samples extracted.")
     rm(inla_samples)
     
    # Extract prediction for each adm
    size_param <- samples[c(1),]
    adm_preds <- matrix(NA, 155, s)
    for(adm in idx){ # For each adm2
      adm_num <- which(idx == adm)
      mean_samples <- samples[adm + 1, ] # Select size parameter and predictor of interest, note this is the n+1th row as the first row is the size parameter
     for(s_id in 1:s){ # For each sample
       adm_preds[adm_num,s_id] <- rnbinom(1,
                                     mu = exp(mean_samples[s_id]),
                                     size = size_param[s_id])
      if (is.na(adm_preds[adm_num, s_id])) {
        message(paste0("NA prediction generated with mu = ", exp(mean_samples[s_id]), " and size = ", size_param[s_id]))
      }
         }
      } 
    pred_full[wk, ,] <- adm_preds 
    rm(adm_preds, samples)
    cat(paste0("Progress: ", wk, "/", max(data_input$date_index), "\n"))
    }
  } 
  tscv_out[["post_samples"]] <- pred_full
  qsave(tscv_out, here(output_path, paste0("tscv-output",  "_",mod_name, ".qs")))
  
  # Calculate outbreak probabilities ---------------------------------------------------------------------------------------------
  cat("Generating model outbreak probabilities.")
  
   outbreak_probs <- matrix(NA, 522,155)
  
  for (wo in start_week:max(data_input$date_index)) {
    idx <- data_input |> 
      mutate(row_number = row_number()) |> 
      filter(date_index == wo) |> 
      pull(row_number)
    for(adm in idx){ 
      adm_num <- which(idx == adm)
    if (!is.na(sum(pred_full[wo, adm_num, ]))) {
      outbreak_probs[wo,adm_num] <- length(pred_full[wo,adm_num,][pred_full[wo,adm_num,] > data_input$threshold[adm]]) / s # What proportion of samples have a prediction > threshold
    } else {
      outbreak_probs[wo, adm_num] <- NA
    }
    }
  }

  tscv_out[["crps"]] <- get_crps(pred_full, data_input)
  print(tscv_out[["crps"]])
  cat("CRPS calculated.")
  tscv_out[["post_samples"]] <- pred_full
  rm(pred_full)
  tscv_out[["outbreak_probs"]] <- outbreak_probs
  tscv_out[["formula"]] <- forms
  tscv_out[["mod_name"]] <- mod_name
  
  qsave(tscv_out, here(output_path, paste0("tscv-output",  "_",mod_name, ".qs")))
  
  return(tscv_out)
}