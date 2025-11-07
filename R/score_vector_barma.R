#' @title Score Vector for a BARMA Model (Internal Version)
#' @description This function computes the score vector (the gradient of the
#'   log-likelihood) for the BARMA model. It is designed to be called with
#'   named parameters, making it easy to test and debug.
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
#' @return A numeric vector representing the score for each parameter.
#'
#' @keywords internal
#'
score_vector_barma <- function(y, ar, ma, alpha, varphi, theta, phi, link) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Validate Precision Parameter ---
  # ------------------------------------------------------------------------- #
  if (phi <= 0 || !is.finite(phi)) {
    warning("phi must be positive and finite; returning zero gradient")
    # Return zero gradient with correct length
    n_params <- 1 + length(varphi) + length(theta) + 1
    return(rep(0, n_params))
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
  mu_eta_fun <- link_structure$mu.eta
  ynew <- linkfun(y)
  
  # Determine maximum lag
  ar_order <- if (has_ar) max(ar_lags) else 0
  ma_order <- if (has_ma) max(ma_lags) else 0
  max_lag  <- max(ar_order, ma_order)
  n_obs <- length(y)
  
  # Check for sufficient observations
  if (n_obs <= max_lag) {
    warning("Insufficient observations for the specified lag structure")
    n_params <- 1 + n_ar_params + n_ma_params + 1
    return(rep(0, n_params))
  }
  
  # ------------------------------------------------------------------------- #
  # --- 4. Recursive Calculation of Predictor, Error, and Derivatives ---
  # ------------------------------------------------------------------------- #
  error <- rep(0, n_obs)
  eta   <- rep(NA_real_, n_obs)
  d_eta_d_alpha  <- rep(0, n_obs)
  d_eta_d_varphi <- if (has_ar) matrix(0, nrow = n_obs, ncol = n_ar_params) else matrix(0, nrow = n_obs, ncol = 0)
  d_eta_d_theta  <- if (has_ma) matrix(0, nrow = n_obs, ncol = n_ma_params) else matrix(0, nrow = n_obs, ncol = 0)
  
  for (t in (max_lag + 1):n_obs) {
    
    # --- Part A: Compute Linear Predictor and Score Error ---
    # ----------------------------------------------------------------------- #
    eta[t] <- alpha
    
    if (has_ar) {
      eta[t] <- eta[t] + as.numeric(crossprod(varphi, ynew[t - ar_lags]))
    }
    
    if (has_ma) {
      eta[t] <- eta[t] + as.numeric(crossprod(theta, error[t - ma_lags]))
    }
    
    # The score error is defined on the predictor scale
    error[t] <- ynew[t] - eta[t]
    
    # --- Part B: Compute Recursive Derivatives of the Linear Predictor ---
    # ----------------------------------------------------------------------- #
    # Derivative w.r.t. alpha
    d_eta_d_alpha[t] <- 1
    if (has_ma) {
      # Recursion: d(eta_t)/d(alpha) = 1 - sum(theta * d(eta_{t-j})/d(alpha))
      d_eta_d_alpha[t] <- 1 - as.numeric(crossprod(theta,
                                                   d_eta_d_alpha[t - ma_lags]))
    }
    
    # Derivative w.r.t. AR coefficients (varphi)
    if (has_ar) {
      # Base term: d(eta_t)/d(varphi_i) starts with ynew_{t-lag_i}
      d_eta_d_varphi[t, ] <- ynew[t - ar_lags]
      if (has_ma) {
        # Recursive term from the MA component
        d_eta_d_varphi[t, ] <- d_eta_d_varphi[t, ] -
          as.numeric(crossprod(theta, d_eta_d_varphi[t - ma_lags, , drop = FALSE]))
      }
    }
    
    # Derivative w.r.t. MA coefficients (theta)
    if (has_ma) {
      # Base term: d(eta_t)/d(theta_j) starts with error_{t-lag_j}
      d_eta_d_theta[t, ] <- error[t - ma_lags]
      # Recursive term from the MA component
      d_eta_d_theta[t, ] <- d_eta_d_theta[t, ] -
        as.numeric(crossprod(theta, d_eta_d_theta[t - ma_lags, , drop = FALSE]))
    }
  }
  
  # ------------------------------------------------------------------------- #
  # --- 5. Calculate Score Vector Components ---
  # ------------------------------------------------------------------------- #
  eta_effective <- eta[(max_lag + 1):n_obs]
  y_effective   <- y[(max_lag + 1):n_obs]
  mu_effective  <- linkinv(eta = eta_effective)
  
  # Bounds check for numerical stability
  if (any(mu_effective <= 0 | mu_effective >= 1 | !is.finite(mu_effective))) {
    warning("mu_effective out of bounds; returning zero gradient")
    n_params <- 1 + n_ar_params + n_ma_params + 1
    return(rep(0, n_params))
  }
  
  # Common component for scores of alpha, AR, and MA
  mu_eta_val <- mu_eta_fun(eta = eta_effective)
  ystar <- linkfun(y_effective)
  mustar <- digamma(mu_effective * phi) - digamma((1 - mu_effective) * phi)
  common_term <- mu_eta_val * (ystar - mustar)
  
  # Score for alpha
  score_alpha <- as.numeric(phi * crossprod(
    d_eta_d_alpha[(max_lag + 1):n_obs], common_term))
  
  # Score for AR coefficients
  score_varphi <- if (has_ar) {
    as.numeric(phi * crossprod(
      d_eta_d_varphi[(max_lag + 1):n_obs, , drop = FALSE], common_term))
  } else {
    numeric(0)
  }
  
  # Score for MA coefficients
  score_theta <- if (has_ma) {
    as.numeric(phi * crossprod(
      d_eta_d_theta[(max_lag + 1):n_obs, , drop = FALSE], common_term))
  } else {
    numeric(0)
  }
  
  # Score for the precision parameter phi
  score_phi <- sum(
    mu_effective * (ystar - mustar) +
      log(1 - y_effective) -
      digamma((1 - mu_effective) * phi) +
      digamma(phi)
  )
  
  # ------------------------------------------------------------------------- #
  # --- 6. Assemble and Return the Final Score Vector ---
  # ------------------------------------------------------------------------- #
  final_score_vector <- c(score_alpha, score_varphi, score_theta, score_phi)
  
  # Check for non-finite values
  if (any(!is.finite(final_score_vector))) {
    warning("Non-finite values in score vector; returning zeros")
    n_params <- 1 + n_ar_params + n_ma_params + 1
    return(rep(0, n_params))
  }
  
  return(as.numeric(final_score_vector))
}