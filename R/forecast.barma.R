#' Forecast a barma Model
#'
#' @description
#' S3 method for producing forecasts from a fitted `"barma"` model object.
#'
#' @details
#' This function computes dynamic, multi-step-ahead point forecasts, where 
#' future predictor-scale errors are assumed to be zero.
#' 
#' The autoregressive part of the forecast uses the conditional expectation
#' \eqn{E[g(y_t) | F_{t-1}] = \eta_t}, so future \eqn{g(y)} values are 
#' replaced by their forecasted \eqn{\eta_t} values.
#'
#' @param object A fitted model object of class `"barma"`.
#' @param h The number of steps to forecast ahead (forecast horizon).
#'   Default is 6.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A `ts` object containing the point forecasts for $h$ steps ahead.
#'
#' @importFrom forecast forecast
#' @export
#' @method forecast barma
forecast.barma <- function(object, h = 6, ...) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Extract Model Components ---
  # ------------------------------------------------------------------------- #
  
  # --- Get parameters ---
  alpha  <- object$alpha
  varphi <- object$varphi
  theta  <- object$theta
  
  # --- Get link function ---
  # link_structure <- do.call("make_link_structure", list(object$link))
  link_structure <- make_link_structure(object$link)
  linkinv <- link_structure$linkinv
  linkfun <- link_structure$linkfun
  
  # --- Get lags ---
  ar_lags <- object$ar_lags
  ma_lags <- object$ma_lags
  has_ar  <- length(ar_lags) > 0
  has_ma  <- length(ma_lags) > 0
  
  # --- Get historical data ---
  y <- object$y
  n_obs <- object$n_obs
  
  # Get g(y)
  ynew <- linkfun(y) 
  
  # Get historical errors (predictor scale)
  errorhat <- object$errorhat
  
  # ------------------------------------------------------------------------- #
  # --- 2. Setup Padded Vectors for Forecasting ---
  # ------------------------------------------------------------------------- #
  
  # Create vectors to hold historical data + future forecasts
  # ynew_padded holds g(y_t) for t <= n_obs and eta_t for t > n_obs
  ynew_padded  <- c(ynew, rep(NA_real_, h))
  error_padded <- c(errorhat, rep(NA_real_, h))
  
  # Create vector to store the final forecasts (response scale)
  forecast_values <- rep(NA_real_, h)
  
  # ------------------------------------------------------------------------- #
  # --- 3. Iterate and Compute Forecasts ---
  # ------------------------------------------------------------------------- #
  
  for (i in 1:h) {
    
    t <- n_obs + i # The time step we are forecasting
    
    # --- Calculate eta[t] ---
    eta_forecast <- alpha
    
    if (has_ar) {
      # Use past values of ynew_padded (which are g(y) or past eta forecasts)
      eta_forecast <- eta_forecast + 
        as.numeric(crossprod(varphi, ynew_padded[t - ar_lags]))
    }
    
    if (has_ma) {
      # Use past values of error_padded (historical errors)
      eta_forecast <- eta_forecast + 
        as.numeric(crossprod(theta, error_padded[t - ma_lags]))
    }
    
    # --- Store forecast ---
    # Convert eta forecast to mu forecast
    mu_forecast <- linkinv(eta_forecast)
    forecast_values[i] <- mu_forecast
    
    # --- Update padded vectors for the *next* iteration ---
    # The expected value of g(y_t) is eta_t
    ynew_padded[t] <- eta_forecast
    
    # The expected value of future errors is 0
    error_padded[t] <- 0 
  }
  
  # ------------------------------------------------------------------------- #
  # --- 4. Format and Return as ts Object ---
  # ------------------------------------------------------------------------- #
  
  # Get time properties from original series
  y_ts <- stats::ts(y)
  ts_start <- stats::tsp(y_ts)[2] + (1 / stats::frequency(y_ts))
  ts_freq <- stats::frequency(y_ts)
  
  forecast_ts <- stats::ts(
    forecast_values,
    start = ts_start,
    frequency = ts_freq
  )
  
  return(forecast_ts)
}