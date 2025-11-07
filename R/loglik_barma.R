#' @title Log-Likelihood for a BARMA Model (Internal Version)
#' @description This function computes the log-likelihood of the BARMA model.
#'  It is designed to be called with named parameters, making it
#'  easy to test and debug.
#'
#' @param y A numeric vector representing the time series data, with values in
#'   (0, 1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#'   Can be NA or NULL if no AR component.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#'   Can be NA or NULL if no MA component.
#' @param alpha The intercept term.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#'   Use numeric(0) or empty vector if no AR component.
#' @param theta A numeric vector of moving average (MA) parameters.
#'   Use numeric(0) or empty vector if no MA component.
#' @param phi The precision parameter of the BARMA model (must be positive).
#' @param link A character string for the link function (e.g., "logit").
#'
#' @return The conditional log-likelihood value of the BARMA model.
#'
#' @importFrom stats dbeta
#' 
#' @keywords internal
loglik_barma <- function(y, ar, ma, alpha, varphi, theta, phi, link) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Validate Precision Parameter ---
  # ------------------------------------------------------------------------- #
  if (phi <= 0 || !is.finite(phi)) {
    return(-Inf)
  }
  
  # ------------------------------------------------------------------------- #
  # --- 2. Determine Model Structure ---
  # ------------------------------------------------------------------------- #
  # Handle cases where ar/ma might be NA, NULL, or empty
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  
  n_ar_params <- length(varphi)
  n_ma_params <- length(theta)
  
  # Get lag specifications (handle empty cases)
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  # Validate parameter-lag consistency
  if (has_ar && n_ar_params != length(ar_lags)) {
    stop("Mismatch between 'ar' lags and 'varphi' parameters.")
  }
  if (has_ma && n_ma_params != length(ma_lags)) {
    stop("Mismatch between 'ma' lags and 'theta' parameters.")
  }
  
  # ------------------------------------------------------------------------- #
  # --- 3. Setup Link Functions and Time Series Properties ---
  # ------------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  ynew <- linkfun(y)
  
  # Determine maximum lag
  ar_order <- if (has_ar) max(ar_lags) else 0
  ma_order <- if (has_ma) max(ma_lags) else 0
  max_lag  <- max(ar_order, ma_order)
  n_obs <- length(y)
  
  # Check for sufficient observations
  if (n_obs <= max_lag) {
    warning("Insufficient observations for the specified lag structure")
    return(-Inf)
  }
  
  # ------------------------------------------------------------------------- #
  # --- 4. Calculate Error and Predictor Iteratively (Optimized) ---
  # ------------------------------------------------------------------------- #
  error <- rep(0, n_obs)
  eta   <- rep(NA_real_, n_obs)
  
  for (t in (max_lag + 1):n_obs) {
    
    eta[t] <- alpha
    
    if (has_ar) {
      eta[t] <- eta[t] + as.numeric(crossprod(varphi, ynew[t - ar_lags]))
    }
    
    if (has_ma) {
      eta[t] <- eta[t] + as.numeric(crossprod(theta, error[t - ma_lags]))
    }
    
    error[t] <- ynew[t] - eta[t]
  }
  
  # ------------------------------------------------------------------------- #
  # --- 5. Calculate the Final Log-Likelihood ---
  # ------------------------------------------------------------------------- #
  eta_effective <- eta[(max_lag + 1):n_obs]
  y_effective   <- y[(max_lag + 1):n_obs]
  mu_effective <- linkinv(eta = eta_effective)
  
  # Check for valid mu values
  if (any(mu_effective <= 0 | mu_effective >= 1 | !is.finite(mu_effective))) {
    return(-Inf)
  }
  
  # Compute log-likelihood using beta distribution
  ll_terms <- dbeta(y_effective,
                    shape1 = mu_effective * phi,
                    shape2 = (1 - mu_effective) * phi,
                    log = TRUE)
  
  # Check for numerical issues
  if (any(!is.finite(ll_terms))) {
    return(-Inf)
  }
  
  final_loglik <- sum(ll_terms)
  return(final_loglik)
}