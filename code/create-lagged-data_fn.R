# -----------------------------------------------------------------------------------------------------------
# Function to create running averages/totals and lags (2 to 20 weeks) of climate covariates
#' @param df data frame with climate variables (tas, tasmin, tasmax, prlr, hurs, huss, absh, nino34, spi) for model selection
#' @return data frame with added running average/total and lagged climate variables
#' Running averages/totals are tagged in variable names with e.g. _4_wk_avg or _4_wk_total
#' Number of weeks a variable is lagged by is tagged in the variable name with e.g. _4 for a four week lag
#' E.g. hurs_scale_4_wk_avg_8 is a four week running average of scaled relative humidity, at an 8 week lag
# -----------------------------------------------------------------------------------------------------------

lag_data <- function(df){
  df <- df |> 
    mutate(tas_scale = tas - mean(tas), tasmin_scale = tasmin - mean(tasmin), tasmax_scale = tasmax - mean(tasmax),
           prlr_scale = prlr - mean(prlr), hurs_scale = hurs - mean(hurs), huss_scale = huss - mean(huss), absh_scale = absh - mean(absh)) 

vars_avg <- c("tas_scale", "tasmin_scale", "tasmax_scale", "hurs_scale", "huss_scale", "absh_scale", "nino34")
vars_sum <- c("prlr_scale")

df <- df |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_avg)), rollmean, k = 4, fill = NA, align = "right"),
                     c(paste0(c(vars_avg), "_4_wk_avg")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_avg)), rollmean, k = 6, fill = NA, align = "right"),
                     c(paste0(c(vars_avg), "_6_wk_avg")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_avg)), rollmean, k = 8, fill = NA, align = "right"),
                   c(paste0(c(vars_avg), "_8_wk_avg")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_avg)), rollmean, k = 10, fill = NA, align = "right"),
                    c(paste0(c(vars_avg), "_10_wk_avg")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_avg)), rollmean, k = 12, fill = NA, align = "right"),
                     c(paste0(c(vars_avg), "_12_wk_avg")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_sum)), rollsum, k = 4, fill = NA, align = "right"),
                     c(paste0(c(vars_sum), "_4_wk_total")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_sum)), rollsum, k = 6, fill = NA, align = "right"),
                     c(paste0(c(vars_sum), "_6_wk_total")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_sum)), rollsum, k = 8, fill = NA, align = "right"),
                     c(paste0(c(vars_sum), "_8_wk_total")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_sum)), rollsum, k = 10, fill = NA, align = "right"),
                     c(paste0(c(vars_sum), "_10_wk_total")))) |> 
  bind_cols(setNames(lapply(df |> st_drop_geometry() |> dplyr::select(all_of(vars_sum)), rollsum, k = 12, fill = NA, align = "right"),
                     c(paste0(c(vars_sum), "_12_wk_total"))))

df <- df |> 
   bind_cols(setNames(shift(df$tas_scale,seq(2,16, by = 2)), c(paste0("tas_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$tasmin_scale,seq(2,16, by = 2)), c(paste0("tasmin_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$tasmax_scale,seq(2,16, by = 2)), c(paste0("tasmax_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$prlr_scale,seq(2,16, by = 2)), c(paste0("prlr_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$hurs_scale,seq(2,16, by = 2)), c(paste0("hurs_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$huss_scale,seq(2,16, by = 2)), c(paste0("huss_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$absh_scale,seq(2,16, by = 2)), c(paste0("absh_scale_", seq(2,16, by = 2))))) |> 
   bind_cols(setNames(shift(df$spi_4,seq(2,20, by = 2)), c(paste0("spi_4_", seq(2,20, by = 2))))) |> 
   bind_cols(setNames(shift(df$spi_12,seq(2,20, by = 2)), c(paste0("spi_12_", seq(2,20, by = 2))))) |> 
   bind_cols(setNames(shift(df$spi_24,seq(2,20, by = 2)), c(paste0("spi_24_", seq(2,20, by = 2))))) |> 
   bind_cols(setNames(shift(df$spi_53,seq(2,20, by = 2)), c(paste0("spi_53_", seq(2,20, by = 2))))) |> 
   bind_cols(setNames(shift(df$nino34,seq(2,20, by = 2)), c(paste0("nino34_", seq(2,20, by = 2))))) |> 
  
  bind_cols(setNames(shift(df$tas_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("tas_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmin_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("tasmin_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmax_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("tasmax_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$prlr_scale_4_wk_total,seq(2,12, by = 2)), c(paste0("prlr_scale_4_wk_total_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$hurs_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("hurs_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$huss_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("huss_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$absh_scale_4_wk_avg,seq(2,12, by = 2)), c(paste0("absh_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
  bind_cols(setNames(shift(df$nino34_4_wk_avg,seq(2,16, by = 2)), c(paste0("nino34_4_wk_avg", seq(2,16, by = 2))))) |> 
  
  bind_cols(setNames(shift(df$tas_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("tas_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmin_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("tasmin_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmax_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("tasmax_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$prlr_scale_6_wk_total,seq(2,10, by = 2)), c(paste0("prlr_scale_6_wk_total_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$hurs_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("hurs_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$huss_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("huss_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$absh_scale_6_wk_avg,seq(2,10, by = 2)), c(paste0("absh_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
  bind_cols(setNames(shift(df$nino34_6_wk_avg,seq(2,14, by = 2)), c(paste0("nino34_6_wk_avg_", seq(2,14, by = 2))))) |> 
  
  bind_cols(setNames(shift(df$tas_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("tas_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmin_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("tasmin_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmax_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("tasmax_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$prlr_scale_8_wk_total,seq(2,8, by = 2)), c(paste0("prlr_scale_8_wk_total_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$hurs_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("hurs_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$huss_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("huss_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$absh_scale_8_wk_avg,seq(2,8, by = 2)), c(paste0("absh_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
  bind_cols(setNames(shift(df$nino34_8_wk_avg,seq(2,12, by = 2)), c(paste0("nino34_8_wk_avg_", seq(2,12, by = 2))))) |> 
  
  bind_cols(setNames(shift(df$tas_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("tas_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmin_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("tasmin_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmax_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("tasmax_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$prlr_scale_10_wk_total,seq(2,6, by = 2)), c(paste0("prlr_scale_10_wk_total_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$hurs_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("hurs_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$huss_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("huss_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$absh_scale_10_wk_avg,seq(2,6, by = 2)), c(paste0("absh_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
  bind_cols(setNames(shift(df$nino34_10_wk_avg,seq(2,10, by = 2)), c(paste0("nino34_10_wk_avg_", seq(2,10, by = 2))))) |> 
  
  bind_cols(setNames(shift(df$tas_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("tas_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmin_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("tasmin_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$tasmax_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("tasmax_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$prlr_scale_12_wk_total,seq(2,4, by = 2)), c(paste0("prlr_scale_12_wk_total_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$hurs_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("hurs_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$huss_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("huss_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
  bind_cols(setNames(shift(df$absh_scale_12_wk_avg,seq(2,4, by = 2)), c(paste0("absh_scale_12_wk_avg_", seq(2,4, by = 2)))))  |> 
  bind_cols(setNames(shift(df$nino34_12_wk_avg,seq(2,8, by = 2)), c(paste0("nino34_12_wk_avg_", seq(2,8, by = 2)))))  
  
  return(df)
}