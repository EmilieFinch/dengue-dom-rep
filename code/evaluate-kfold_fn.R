# ----------------------------------------------------------------------------------------------------
# Function to run k-fold spatial cross validation with INLA
#' @param data_input data frame for modelling
#' @param nb_graph neighbourhood graph for spatial random effects
#' @param forms character string with INLA formula for covariates
#' @param mod_name character string for short hand model name
#' @param k_folds integer, number of folds
#' @param rep integer, which repetition
#' @param output_path character string of output file path

# ----------------------------------------------------------------------------------------------------

evaluate_kfold <- function(data_input,
                          nb_graph,
                          forms,
                          mod_name,
                          k_folds,
                          rep,
                          output_path){
  
  # Define objects for model output
  posterior_samples <- NULL
  outbreak_probs <- NULL
  
  # Set up parameters
  s <- 1000
  pred_full <- array(dim = c(522, 155, s))
  kfold_out <- list()
  
  # Set up random effects
  eweek_re <- " + f(epiweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
  adm2_year_re <- "+f(adm2_id, model = 'bym2', graph = nb_graph, scale.model = TRUE, hyper = precision_prior, replicate = year_index)"
  res <- paste0(eweek_re, adm2_year_re) # What to do for baseline model - seasonal and spatial only?
  
  # Set up prior for random effects
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  cat(paste0("Running K-fold evaluation with formula:", forms, res, "\n"))

  
  for (kf in 1:k_folds){ 
  cat(paste0("Running fold:", kf, "\n"))
    
    df_pred <- data_input |> 
      mutate(reported = case_when(kfold == kf ~ NA_integer_, TRUE ~ reported)) # set cases to NA for kfold municipalities
    
    # Get idx for selected k-fold adm2s to predict
    k_adms <- df_pred |> 
      mutate(row_number = row_number()) |> 
      filter(kfold == kf) |> 
      select(adm2_id, row_number)  
    
    cat("Starting model fitting.")
    
    mod  <- inla(as.formula(paste0("reported~1", res, forms)),
                 family = "nbinomial",
                 #offset = log(population / 100000),
                 control.inla = list(strategy = "adaptive"),
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
                 control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
                 num.threads = 20,
                 verbose = FALSE,
                 data = df_pred)
    
    cat("Model fitting completed.")
    
    inla_samples <- inla.posterior.sample(s, mod)
    rm(mod)
    
    cat("Samples generated from posterior.")
    
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
    wk_preds <- matrix(NA, 522, s)
    
    for(adm in unique(k_adms$adm2_id)) { # For each selected adm
      curr_adms <- k_adms |> 
        filter(adm2_id == adm)
    for(obs in curr_adms$row_number){ # Go down the rows associated with that adm (1 for each week)
      wk_num <- which(curr_adms$row_number == obs)
      mean_samples <- samples[obs + 1, ] # Select size parameter and predictor of interest, note this is the n+1th row as the first row is the size parameter
      for(s_id in 1:s){ # Calculate prediction for each sample
        wk_preds[wk_num,s_id] <- rnbinom(1,
                                           mu = exp(mean_samples[s_id]),
                                           size = size_param[s_id])
        if (is.na(wk_preds[wk_num, s_id])) {
          message(paste0("NA prediction generated with mu = ", exp(mean_samples[s_id]), " and size = ", size_param[s_id]))
        }
      }
    } 
      pred_full[,adm,] <- wk_preds 
      cat(paste0("Progress: Admin ID ", adm, " completed. \n"))
      
      
    }
    rm(wk_preds, samples)
    
  } 
  
  kfold_out[["post_samples"]] <- pred_full
  qsave(kfold_out, here(output_path, paste0("kfold-output",  "_",mod_name, "-rep-",rep, ".qs")))
  
  cat("Posterior samples saved.")
  
  # Calculate outbreak probabilities ---------------------------------------------------------------------------------------------
  
  cat("Generating model outbreak probabilities.")
  outbreak_probs <- matrix(NA, 522,155)
  
  for (wo in 1:max(data_input$date_index)) {
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
  
  kfold_out[["crps"]] <- get_crps(pred_full, data_input)
  cat("CRPS calculated.")
  kfold_out[["post_samples"]] <- pred_full
  rm(pred_full)
  kfold_out[["outbreak_probs"]] <- outbreak_probs
  rm(outbreak_probs)
  kfold_out[["formula"]] <- forms
  kfold_out[["mod_name"]] <- mod_name
  kfold_out[["rep"]] <- rep

  
  qsave(kfold_out, here(output_path, paste0("kfold-output",  "_",mod_name, "-rep-",rep, ".qs")))
  cat("Full output saved.")
  return(kfold_out)
  
}