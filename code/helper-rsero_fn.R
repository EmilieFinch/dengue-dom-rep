# Helper fns for RSero

# -----------------------------------------------------------------------------------------------------------
# Function to strip a fitted FOIfit object down to only the stanfit parameters needed for downstream
# plotting and model comparison (see 06_plot-catalytic.R), dropping the rest and the warmup draws.
# The full fit stores every parameter (including ones only needed while fitting, e.g. Flambda2, x, P1) for
# every warmup and post-warmup iteration, which is what makes the saved .rds files hundreds of MB.
#' @param FOIfit fitted FOIfit object as returned by Rsero::fit()
#' @param keep_pars character vector of stanfit parameter names to retain.
#'   Defaults to what compute_dic_waic() and plot_fit() need: "P" and "Likelihood",
#'   plus "annual_foi" for the convergence trace plot
#' @return FOIfit object with the same data/model but a much smaller $fit slot
# -----------------------------------------------------------------------------------------------------------

trim_serofit <- function(FOIfit, keep_pars = c("annual_foi", "P", "Likelihood")){
  sf <- FOIfit$fit
  fnames <- sf@sim$fnames_oi
  base <- sub("\\[.*\\]$", "", fnames)
  keep_idx <- which(base %in% keep_pars)
  ordered_keep_pars <- unique(base[keep_idx]) # must match the original flat storage order, not keep_pars order

  sf@sim$samples <- Map(function(chain, warmup_n){
    kept <- lapply(chain[keep_idx], function(v) v[(warmup_n + 1):length(v)])
    attrs <- attributes(chain)
    attrs$names <- attrs$names[keep_idx]
    if(!is.null(attrs$mean_pars)) attrs$mean_pars <- attrs$mean_pars[keep_idx]
    attrs$warmup2 <- 0
    attributes(kept) <- attrs
    kept
  }, sf@sim$samples, sf@sim$warmup2)

  sf@sim$fnames_oi <- fnames[keep_idx]
  sf@sim$dims_oi <- sf@sim$dims_oi[ordered_keep_pars]
  sf@sim$pars_oi <- ordered_keep_pars
  sf@sim$n_flatnames <- length(keep_idx)
  sf@model_pars <- ordered_keep_pars
  sf@par_dims <- sf@par_dims[ordered_keep_pars]
  sf@sim$n_save <- sf@sim$n_save - sf@sim$warmup2
  sf@sim$warmup2 <- rep(0, sf@sim$chains)
  sf@sim$warmup <- 0
  sf@sim$iter <- sf@sim$n_save[1]

  FOIfit$fit <- sf
  return(FOIfit)
}

# DIC/WAIC only, skipping Rsero::compute_information_criteria's PSIS_LOO step.
# That step calls loo::loo() on the raw stanfit, which requires a "log_lik"
# generated quantity - older Rsero model versions (e.g. fits saved before the
# Stan model was updated) only expose "Likelihood", so loo() errors with
# "no parameter log_lik" even though DIC/WAIC can still be computed manually.
compute_dic_waic <- function(FOIfit){
  estimated_parameters <- FOIfit$model$estimated_parameters
  chains <- rstan::extract(FOIfit$fit)
  sensitivity <- FOIfit$model$se
  specificity <- FOIfit$model$sp
  S <- nrow(chains$Likelihood)
  N <- FOIfit$data$N
  Y <- FOIfit$data$Y
  category <- FOIfit$data$categoryindex

  LogLikelihoods <- matrix(0, nrow = S, ncol = N)
  for(s in seq(1, S)){
    lk <- chains$Likelihood[s, ]
    for(i in seq(1, N)){
      LogLikelihoods[s, i] <- if(Y[i] == FALSE) log(1 - lk[i]) else log(lk[i])
    }
  }

  LogLikelihoodMean <- 0
  P <- colMeans(chains$P)
  for(i in seq(1, N)){
    age <- FOIfit$data$age[i]
    age_group <- FOIfit$data$age_group[i]
    cat <- category[i]
    p <- P[age, age_group, cat]
    LogLikelihoodMean <- LogLikelihoodMean + if(Y[i] == TRUE){
      log(sensitivity - p * (sensitivity + specificity - 1))
    } else {
      log(1 - sensitivity + p * (sensitivity + specificity - 1))
    }
  }

  LP <- rowSums(LogLikelihoods)
  AIC <- -2 * max(LP) + 2 * estimated_parameters
  Dbar <- -2 * mean(LP)
  Dthetabar <- -2 * LogLikelihoodMean
  pD <- Dbar - Dthetabar
  DIC <- pD + Dbar
  V <- Rsero:::ColVar(exp(LogLikelihoods))
  pwaic <- sum(V)
  lpd <- sum(log(colSums(exp(LogLikelihoods)) / S))
  WAIC <- -2 * (lpd - pwaic)

  information_criteria <- list(AIC = AIC, DIC = DIC, pD = pD, Dbar = Dbar,
                               WAIC = WAIC, pwaic = pwaic, lpd = lpd,
                               k = estimated_parameters, MLE = max(LP))
  class(information_criteria) <- "information_criteria"
  return(information_criteria)
}

extract_fois <- function(fit){
  annual_foi <- data.frame()
  chains <- rstan::extract(fit$fit)
  L1 <- chains$lambda
  
  for(cat in fit$data$unique.categories){ 
    temp <- NULL
    w = which( fit$data$category==cat, arr.ind = TRUE)[,1]
    
    d = fit$data$categoryindex[w]
    p1=proportions.index(d)
    
    M=dim(chains$P)[1] 
    L=matrix(0, nrow = M, ncol=fit$data$A)
    # weighted average of the force of infection for each subcategory
    
    for(i in 1:length(p1$index)){
      L =  L+  p1$prop[i]*chains$Flambda[, p1$index[i]]*L1 
    }
    
    # compute mean and 95% credible interval  - note one for each year
    par_out <- apply(L, 2, function(x)c(mean(x), quantile(x, probs=c(0.025, 0.5, 0.975))))
    temp$mean <- par_out[1,1]
    temp$median <- par_out[3,1]
    temp$conf_0_025 <- par_out[2,1]
    temp$conf_0_975 <- par_out[4,1]
    temp$group <- cat
    annual_foi <- rbind(annual_foi, temp)
  }
  return(annual_foi)
}

get_attack_rate <- function(fois){
  attack_rates<- fois %>% 
    filter(group == group) %>% 
    mutate(mean = 1 - exp(1)^-mean, median = 1 - exp(1)^-median, conf_0_025 = 1 -exp(1)^-conf_0_025, conf_0_975 = 1 -exp(1)^-conf_0_975) %>% 
    mutate(mean = mean*100, median = median*100, conf_0_025 = conf_0_025*100, conf_0_975 = conf_0_975*100)
  return(attack_rates)
}


plot_fit <- function (FOIfit, individual_samples = 0, age_class = 10, YLIM = 1, 
                      ...) 
{
  plots <- NULL
  data = FOIfit$data
  chains <- rstan::extract(FOIfit$fit)
  se = FOIfit$model$se
  sp = FOIfit$model$sp
  A <- FOIfit$data$A
  latest_sampling_year <- max(FOIfit$data$sampling_year)
  years <- seq(1, A)
  index.plot = 0
  unique.categories = data$unique.categories
  sorted.year = sort.int(unique(FOIfit$data$sampling_year), 
                         index.return = TRUE)
  Y = 0
  for (sampling_year in sorted.year$x) {
    Y = Y + 1
    for (cat in unique.categories) {
      if (length(unique.categories) == 1) {
        title = NULL
      }
      if (length(unique.categories) > 1) {
        title = paste0(cat)
      }
      index.plot = index.plot + 1
      age_group = data$age_group[which(data$sampling_year == 
                                         sampling_year)][1]
      w = which(data$sampling_year == sampling_year & data$category == 
                  cat, arr.ind = TRUE)[, 1]
      subdat = subset(data, sub = w)
      P = chains$P[, , sorted.year$ix[Y], 1]
      d = data$categoryindex[w]
      p1 = proportions.index(d)
      M = dim(chains$P)[1]
      Pinf = matrix(0, nrow = M, ncol = FOIfit$data$A)
      for (i in 1:length(p1$index)) {
        Pinf = Pinf + p1$prop[i] * (se - (se + sp - 1) * 
                                      chains$P[, , sorted.year$ix[Y], p1$index[i]])
      }
      par_out <- apply(Pinf, 2, function(x) c(mean(x), 
                                              quantile(x, probs = c(0.025, 0.975))))
      par_out[par_out > YLIM] = YLIM
      years.plotted = seq(latest_sampling_year - sampling_year + 
                            1, dim(chains$P)[2])
      years.plotted.normal = years.plotted - min(years.plotted) + 
        1
      meanFit <- data.frame(x = years.plotted.normal, y = par_out[1, 
                                                                  years.plotted])
      xpoly <- (c(years.plotted.normal, rev(years.plotted.normal)))
      ypoly <- c(par_out[3, years.plotted], rev(par_out[2, 
                                                        years.plotted]))
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      histdata <- sero.age.groups(dat = subdat, age_class = age_class, 
                                  YLIM = YLIM)
      histdata$labels <- str_replace_all(histdata$labels, " ", "")
      histdata$labels <- str_replace_all(histdata$labels, ">=", "")
      histdata$labels[max(nrow(histdata))] <- paste0(histdata$labels[max(nrow(histdata))], "+")
      histdata$age[max(nrow(histdata))] <- 82.5
      
      histdata$labels_text <- as.character(histdata$labels)
      last_row_index <- tail(which(complete.cases(histdata)), 
                             1)
      if (length(last_row_index) > 0) {
        histdata = histdata[1:last_row_index, ]
      }
      max.age = max(histdata$age)
      DataEnvelope = subset(DataEnvelope, x <= max.age)
      meanFit = subset(meanFit, x <= max.age)
      p <- ggplot2::ggplot() + ggplot2::geom_polygon(data = DataEnvelope, 
                                                     ggplot2::aes(x, y), fill = "#016c59", alpha = 0.35) + ggplot2::geom_line(data = meanFit, 
                                                                                                                              ggplot2::aes(x = x, y = y), size = 0.8, color = "#016c59")
      if (individual_samples > 0) {
        Index_samples <- sample(nrow(Pinf), individual_samples)
        for (i in Index_samples) {
          ind_foi <- data.frame(x = years.plotted.normal, 
                                y = Pinf[i, years.plotted])
          p <- p + ggplot2::geom_line(data = ind_foi, 
                                      ggplot2::aes(x = x, y = y), size = 0.8, colour = "#bbbbbb", 
                                      alpha = 0.6)
        }
      }
      p <- p + scale_x_continuous(breaks = histdata$age, 
                                  labels = histdata$labels_text) + geom_point(data = histdata, 
                                                                              aes(x = age, y = mean), size = 1) + geom_segment(data = histdata, 
                                                                                                                     aes(x = age, y = lower, xend = age, yend = upper), size = 0.4) + 
        ggplot2::xlab("Age (years)") + ggplot2::ylab("Seroprevalence") + 
        ggtitle(title) +
        theme_classic() + theme(axis.text.x = element_text(family = plot_font, 
                                                           size = 11), 
                                axis.text.y = element_text(family = plot_font, size = 12), 
                                axis.title = element_text(family = plot_font, size = 13),
                                text = element_text(family = plot_font, size = 12),
                                title = element_text(family = plot_font, size = 12)) + ylim(0, YLIM) 
      plots[[index.plot]] <- p
      plots[[index.plot]]$category <- cat
      plots[[index.plot]]$year <- sampling_year
    }
  }
  return(plots)
}
