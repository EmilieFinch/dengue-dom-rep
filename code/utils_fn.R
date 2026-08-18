# -----------------------------------------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# Calculate outbreak thresholds
#' @param data_input data frame for calculation
#' @param month_input numeric month to filter by
#' @param year_input numeric year to filter by
#' @param adm_input numeric administrative unit ID
#' @param quantile numeric quantile to calculate
# -----------------------------------------------------------------------------------------------------------

calculate_thresholds <- function(data_input, month_input, year_input, adm_input, quantile) {
  # Using only data prior to the year/month of interest
  threshold <- data_input |>
    filter(year_index < as.numeric(year_input), month == month_input, adm2_id == adm_input) |>
    mutate(threshold = quantile(reported, quantile, na.rm = TRUE)) |>
    pull(threshold)
  
  threshold_out <- data.frame(year_index = year_input, month = month_input, adm2_id = adm_input, threshold = unique(threshold))
  
  return(threshold_out)
}

# -----------------------------------------------------------------------------------------------------------
# Calculate weighted lagged cases
#' Discretised probability distribution estimated in Siraj et al 2017
#' Maximum serial interval of 5 weeks (consistent with mosquito lifespan)
#' @param x numeric vector of cases
# -----------------------------------------------------------------------------------------------------------

weight_cases <- function(x){
  weighted_cases <- lag(x,2)*0.2 + lag(x,3)*0.425 + lag(x,4)*0.25 + lag(x,5)*0.125
  return(weighted_cases)
}

# -----------------------------------------------------------------------------------------------------------
# K-fold partition function (adapted from dismo package)
#' @param x numeric vector of observations
#' @param k integer number of folds
#' @param by optional grouping variable
# -----------------------------------------------------------------------------------------------------------

kfold_func <- function (x, k = 5, by = NULL) 
{
  singlefold <- function(obs, k) {
    if (k == 1) {
      return(rep(1, obs))
    }
    else {
      i <- obs/k
      if (i < 1) {
        stop("insufficient records:", obs, ", with k=", 
             k)
      }
      i <- round(c(0, i * 1:(k - 1), obs))
      times = i[-1] - i[-length(i)]
      group <- c()
      for (j in 1:(length(times))) {
        group <- c(group, rep(j, times = times[j]))
      }
      r <- order(runif(obs))
      return(group[r])
    }
  }
  if (is.vector(x)) {
    if (length(x) == 1) {
      if (x > 1) {
        x <- 1:x
      }
    }
    obs <- length(x)
  }
  else if (inherits(x, "Spatial")) {
    if (inherits(x, "SpatialPoints")) {
      obs <- nrow(coordinates(x))
    }
    else {
      obs <- nrow(x@data)
    }
  }
  else {
    obs <- nrow(x)
  }
  if (is.null(by)) {
    return(singlefold(obs, k))
  }
  by = as.vector(as.matrix(by))
  if (length(by) != obs) {
    stop("by should be a vector with the same number of records as x")
  }
  un <- unique(by)
  group <- vector(length = obs)
  for (u in un) {
    i = which(by == u)
    kk = min(length(i), k)
    if (kk < k) 
      warning("lowered k for by group: ", u, "  because the number of observations was  ", 
              length(i))
    group[i] <- singlefold(length(i), kk)
  }
  return(group)
}


# Function to calculate CRPS
get_crps <- function(post_samples, data_input) {
  
  tryCatch( 
    { 
  preds <- reshape2::melt(post_samples) |> as.data.table()
  setnames(preds, new = c("date_index", "adm2_id", "sample", "prediction"))  
  
  data_input <- as.data.table(data_input)
  data_input <- data_input[, .(reported, adm2_id, date_index, date)]
  preds <- merge(preds, data_input, by = c("adm2_id", "date_index"), all.x = TRUE)

    scores <- preds|>
      mutate(target_end_date = date, true_value = reported) |>
      transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
      filter(!is.na(prediction)) |>
      score() |>
      summarise_scores(by = c("scale")) |> 
      select(scale, crps) |> 
      pivot_wider(names_from = scale, names_prefix = "crps_", values_from = crps)
    return(scores)
    
    },
  error = function(cond){NA})
}

# Extract posterior predictions from model evaluation

## For adm2 
get_adm2_preds  <- function(preds_input, data_input, mod_name) {
  
  adm2_preds <- data.table()
  for(rep in 1:length(preds_input)){ 
  preds <- reshape2::melt(preds_input[[rep]]$post_samples) |> as.data.table()
  setnames(preds, new = c("date_index", "adm2_id", "sample", "predicted"))  
 
  preds[,sample := sample + 1000*(rep-1)] # update sample number to account for multiple reps
  
  adm2_preds <-rbindlist(list(adm2_preds, preds))
  }
  
  adm2_preds <- adm2_preds[, .(
    median = median(predicted, na.rm = TRUE),
    q_0.025 = quantile(predicted, 0.025, na.rm = TRUE),
    q_0.975 = quantile(predicted, 0.975, na.rm = TRUE),
    q_0.25 = quantile(predicted, 0.25, na.rm = TRUE),
    q_0.75 = quantile(predicted, 0.75, na.rm = TRUE)
  ), by = .(adm2_id, date_index)]
  
  data_input <- as.data.table(data_input)
  data_input <- data_input[, .(reported, threshold, adm2_id, date_index, date, adm2_name, adm1_name)]
  
  # Performing left join
  adm2_preds <- merge(adm2_preds, data_input, by = c("adm2_id", "date_index"), all.x = TRUE)
  adm2_preds$mod_name <- mod_name
  
  
  return(adm2_preds)
}

## Extract posterior predictions and aggregate to adm1 level
get_adm1_preds <- function(preds_input, data_input, mod_name){
  
  adm1_preds <- data.table()
  for(rep in 1:length(preds_input)){ 
    preds <- reshape2::melt(preds_input[[rep]]$post_samples) |> as.data.table()
    setnames(preds, new = c("date_index", "adm2_id", "sample", "predicted"))  
    
    preds[,sample := sample + 1000*(rep-1)] # update sample number to account for multiple reps
    
    adm1_preds <-rbindlist(list(adm1_preds, preds))
  }
  
  adm_ids <- data_input |> select(adm2_id, adm1_name) |> unique() |> as.data.table()
  setkey(adm_ids, "adm2_id")
  setkey(adm1_preds, "adm2_id")
  
  adm1_preds <- adm1_preds[adm_ids]
  
  adm1_preds <- adm1_preds[, .(
    predicted = sum(predicted)
  ), by = .(adm1_name, sample, date_index)] 
  
  adm1_preds <- adm1_preds[, .(
    median = median(predicted, na.rm = TRUE),
    q_0.025 = quantile(predicted, 0.025, na.rm = TRUE),
    q_0.975 = quantile(predicted, 0.975, na.rm = TRUE),
    q_0.25 = quantile(predicted, 0.25, na.rm = TRUE),
    q_0.75 = quantile(predicted, 0.75, na.rm = TRUE)
    
  ), by = .(adm1_name, date_index)]
  
  # Add relevant data columns
  data_input <- as.data.table(data_input)
  data_input <- data_input[, .(
    reported = sum(reported, na.rm = TRUE)
  ), by = .(adm1_name, date, date_index)] 
  
  # Performing left join
  adm1_preds <- merge(adm1_preds, data_input, by = c("adm1_name", "date_index"), all.x = TRUE)
  adm1_preds$mod_name <- mod_name
  
  return(adm1_preds)
}

# Score predictions
score_predictions <- function(preds_input, data_input, output_path, mod_name, adm_level) {

  preds_input <- preds_input |>
    rename(q_0.5 = median, model = mod_name) |>
    mutate(target_end_date = date, true_value = reported) |>
    pivot_longer(cols = starts_with("q_"), names_to = "quantile", values_to = "prediction") |>
    mutate(quantile = as.numeric(str_remove(quantile, "q_")))

  scores <- preds_input |>
        transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
        filter(!is.na(prediction)) |>
        score() |>
        summarise_scores(fun = round, digits = 3)

  scores_week <- preds_input |>
    transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
    filter(!is.na(prediction)) |>
    score() |>
    add_coverage(by = c("scale", "date"), ranges = c(50,95)) |>
    summarise_scores(by = c("scale", "date"), na.rm = TRUE) |>
    summarise_scores(fun = round, digits = 3)

  scores_adm <- preds_input |>
    transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
    filter(!is.na(prediction)) |>
    score() |>
    add_coverage(by = c("scale", adm_level), ranges = c(50,95)) |>
    summarise_scores(by = c("scale", adm_level), na.rm = TRUE) |>
    summarise_scores(fun = round, digits = 3)

  # Plots over time
  
  ## Coverage plot
  
  coverage_plot <- scores_week |> 
    filter(scale == "natural") |> 
    rename("50%" = coverage_50, "95%" = coverage_95) |> 
    pivot_longer(cols = ends_with("%"), names_to = "coverage", values_to = "value") |> 
    ggplot() +
    geom_line(aes(x = date, y = value, linetype = coverage), col = "#016c59") + 
    geom_hline(aes(yintercept = 0.5), lty = "longdash") +
    geom_hline(aes(yintercept = 0.95), lty = "longdash") +
    scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
    labs(x = "Date", y = "Interval coverage") +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ## Bias plot
  
  bias_plot <- scores_week |> 
    filter(scale == "natural") |> 
    ggplot() +
    geom_line(aes(x = date, y = bias), col = "#016c59") + 
    geom_hline(aes(yintercept = 0), lty = "longdash") +
    scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
    labs(x = "Date", y = "Bias") +
    theme_classic()
  
  ## WIS plot
  
  wis_plot <- scores_week |> 
    filter(scale == "natural") |> 
    ggplot() +
    geom_line(aes(x = date, y = interval_score), col = "#016c59") + 
    scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
    labs(x = "Date", y = "Weighted interval score") +
    theme_classic()
  
  ## log WIS plot
  
  rwis_plot <- scores_week |> 
    filter(scale == "log") |> 
    ggplot() +
    geom_line(aes(x = date, y = interval_score), col = "#016c59") + 
    scale_x_date(date_breaks = "1 year", labels = scales::label_date_short()) +
    labs(x = "Date", y = "Relative weighted interval score") +
    theme_classic()
  
  # Plots by adm
  
  # WIS
  
  wis_plot_adm <- scores_adm |> 
    filter(scale == "natural") |> 
    ggplot() +
    geom_col(aes(x=!! rlang::sym(adm_level), y = interval_score), fill = "#016c59") + 
    labs(x = "", y = "Weighted interval score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95))
  
  # Bias
  
  bias_plot_adm <- scores_adm |> 
    filter(scale == "natural") |> 
    ggplot() +
    geom_col(aes(x=!! rlang::sym(adm_level), y = bias), fill = "#016c59") + 
    labs(x = "", y = "Bias") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95))
  
  # rWIS
  
  rwis_plot_adm <- scores_adm |> 
    filter(scale == "log") |> 
    ggplot() +
    geom_col(aes(x=!! rlang::sym(adm_level), y = interval_score), fill = "#016c59") + 
    labs(x = "", y = "Relative weighted interval score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95))
  
  score_plots <- ggarrange(coverage_plot, bias_plot, wis_plot, rwis_plot, nrow = 4, heights = c(1.2, 1, 1, 1))
  score_plots_adm <- ggarrange(bias_plot_adm, wis_plot_adm, rwis_plot_adm, nrow = 3)
  
  ggsave(here(output_path, paste0("scoring-plots_", mod_name,".png")), score_plots, width = 16, height = 12)
  ggsave(here(output_path, paste0("scoring-plots-adm_", mod_name, ".png")), score_plots_adm, width = 20, height = 15)
  write.csv(scores, here(output_path, paste0("scoring-table_", mod_name, ".csv")))
      
    }