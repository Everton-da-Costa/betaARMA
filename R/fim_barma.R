#' @title FIM and Fitted Values for a BARMA Model (Internal)
#' @description This function computes the Fisher information matrix of the
#'   Beta Autoregressive Moving Average (BARMA) model. It also efficiently
#'   returns auxiliary values like fitted values and residuals required
#'   by the main 'barma' function.
#'
#' @param y A numeric vector representing the time series data, in (0, 1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#' @param alpha The intercept term.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#' @param theta A numeric vector of moving average (MA) parameters.
#' @param phi The precision parameter of the BARMA model.
#' @param link The link function to be used.
#'
#' @return A list with the following elements:
#'   - 'fisher_info_mat': The Fisher information matrix.
#'   - 'fitted_ts': The fitted values (conditional mean) as a 'ts' object.
#'   - 'muhat_effective': The effective vector of fitted conditional means.
#'   - 'etahat_full': The full vector (NA-padded) of the linear predictor.
#'   - 'errorhat_full': The full vector of errors on the predictor scale.
#' @keywords internal
fim_barma <- function(y,
                      ar,
                      ma,
                      alpha,
                      varphi,
                      theta,
                      phi,
                      link) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Setup Link Functions and Model Structure ---
  # ------------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu_eta_fun <- link_structure$mu.eta
  ynew <- linkfun(y)
  
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  ar_order <- if (has_ar) max(ar_lags) else 0L
  ma_order <- if (has_ma) max(ma_lags) else 0L
  max_lag  <- max(ar_order, ma_order)
  
  n_obs <- length(y)
  
  # ------------------------------------------------------------------------- #
  # --- 2. Recursive Calculation of Predictor, Error, and Derivatives ---
  # ------------------------------------------------------------------------- #
  error   <- rep(0, n_obs)
  eta     <- rep(NA_real_, n_obs)
  d_eta_d_alpha  <- rep(0, n_obs)
  d_eta_d_varphi <- matrix(0, nrow = n_obs, ncol = n_ar_params)
  d_eta_d_theta  <- matrix(0, nrow = n_obs, ncol = n_ma_params)
  
  for (t in (max_lag + 1):n_obs) {
    # Part A: Compute Linear Predictor and Score Error
    eta[t] <- alpha
    if (has_ar) {
      eta[t] <- eta[t] + as.numeric(crossprod(varphi, ynew[t - ar_lags]))
    }
    if (has_ma) {
      eta[t] <- eta[t] + as.numeric(crossprod(theta, error[t - ma_lags]))
    }
    error[t] <- ynew[t] - eta[t]
    
    # Part B: Compute Recursive Derivatives
    d_eta_d_alpha[t] <- 1
    if (has_ma) {
      d_eta_d_alpha[t] <- 1 - as.numeric(crossprod(theta,
                                                   d_eta_d_alpha[t - ma_lags]))
    }
    if (has_ar) {
      d_eta_d_varphi[t, ] <- ynew[t - ar_lags]
      if (has_ma) {
        d_eta_d_varphi[t, ] <- d_eta_d_varphi[t, ] -
          as.numeric(crossprod(theta, 
                               d_eta_d_varphi[t - ma_lags, , drop = FALSE]))
      }
    }
    if (has_ma) {
      d_eta_d_theta[t, ] <- error[t - ma_lags]
      d_eta_d_theta[t, ] <- d_eta_d_theta[t, ] -
        as.numeric(crossprod(theta, 
                             d_eta_d_theta[t - ma_lags, , drop = FALSE]))
    }
  }
  
  # ------------------------------------------------------------------------- #
  # --- 3. Get Effective (non-NA) Observations ---
  # ------------------------------------------------------------------------- #
  eta_effective <- eta[(max_lag + 1):n_obs]
  mu_effective  <- linkinv(eta = eta_effective)
  
  if (any(mu_effective <= 0 | mu_effective >= 1 | !is.finite(mu_effective))) {
    warning("mu_effective out of bounds during FIM calculation")
    # Return an empty matrix of the correct size
    n_params <- 1 + n_ar_params + n_ma_params + 1
    return(list(
      fisher_info_mat = matrix(NA_real_, nrow = n_params, ncol = n_params)
    ))
  }
  
  # Get effective derivatives (s, P, and R from erratum)
  s_vec <- d_eta_d_alpha[(max_lag + 1):n_obs]
  P_mat <- d_eta_d_varphi[(max_lag + 1):n_obs, , drop = FALSE]
  R_mat <- d_eta_d_theta[(max_lag + 1):n_obs, , drop = FALSE]
  
  # ------------------------------------------------------------------------- #
  # --- 4. Calculate FIM Component Vectors (Efficiently) ---
  # ------------------------------------------------------------------------- #
  trigamma_p <- trigamma(mu_effective * phi)
  trigamma_q <- trigamma((1 - mu_effective) * phi)
  mu_eta_val <- mu_eta_fun(eta = eta_effective)
  
  # w_t = phi * {psi'(mu*phi) + psi'((1-mu)*phi)} * (mu_eta)^2
  w_t_vec <- phi * (trigamma_p + trigamma_q) * mu_eta_val^2
  
  # c_t = phi * {psi'(mu*phi)*mu - psi'((1-mu)*phi)*(1-mu)}
  c_t_vec <- phi * (trigamma_p * mu_effective - 
                      trigamma_q * (1 - mu_effective))
  
  # d_t = psi'(mu*phi)*mu^2 + psi'((1-mu)*phi)*(1-mu)^2 - psi'(phi)
  d_t_vec <- trigamma_p * mu_effective^2 + 
    trigamma_q * (1 - mu_effective)^2 - 
    trigamma(phi)
  
  # ------------------------------------------------------------------------- #
  # --- 5. Calculate FIM Blocks ---
  # ------------------------------------------------------------------------- #
  K_aa <- phi * as.numeric(crossprod(s_vec, w_t_vec * s_vec))
  
  K_pa <- if (has_ar) {
    phi * as.numeric(crossprod(P_mat, w_t_vec * s_vec))
  } else {
    numeric(0)
  }
  
  K_ta <- if (has_ma) {
    phi * as.numeric(crossprod(R_mat, w_t_vec * s_vec))
  } else {
    numeric(0)
  }
  
  K_ap <- as.numeric(crossprod(s_vec, mu_eta_val * c_t_vec))
  
  K_pp <- if (has_ar) {
    phi * crossprod(P_mat, P_mat * w_t_vec)
  } else {
    matrix(0, nrow = 0, ncol = 0)
  }
  
  K_tt <- if (has_ma) {
    phi * crossprod(R_mat, R_mat * w_t_vec)
  } else {
    matrix(0, nrow = 0, ncol = 0)
  }
  
  K_pt <- if (has_ar && has_ma) {
    phi * crossprod(P_mat, R_mat * w_t_vec)
  } else {
    matrix(0, nrow = n_ar_params, ncol = n_ma_params)
  }
  
  K_pphi <- if (has_ar) {
    as.numeric(crossprod(P_mat, mu_eta_val * c_t_vec))
  } else {
    numeric(0)
  }
  
  K_tphi <- if (has_ma) {
    as.numeric(crossprod(R_mat, mu_eta_val * c_t_vec))
  } else {
    numeric(0)
  }
  
  K_phiphi <- sum(d_t_vec)
  
  # ------------------------------------------------------------------------- #
  # --- 6. Assemble and Return the Final FIM and Other Values ---
  # ------------------------------------------------------------------------- #
  n_params <- 1 + n_ar_params + n_ma_params + 1
  fim <- matrix(NA_real_, nrow = n_params, ncol = n_params)
  
  # Define indices
  idx_a <- 1
  idx_p <- if (n_ar_params > 0) 2:(1 + n_ar_params) else integer(0)
  idx_t <- if (n_ma_params > 0) {
    (2 + n_ar_params):(1 + n_ar_params + n_ma_params)
  } else {
    integer(0)
  }
  idx_phi <- n_params
  
  # Assign blocks (FIM is symmetric)
  fim[idx_a, idx_a] <- K_aa
  fim[idx_phi, idx_phi] <- K_phiphi
  
  fim[idx_phi, idx_a] <- K_ap
  fim[idx_a, idx_phi] <- K_ap
  
  if (has_ar) {
    fim[idx_p, idx_p] <- K_pp
    fim[idx_p, idx_a] <- K_pa
    fim[idx_a, idx_p] <- K_pa
    fim[idx_p, idx_phi] <- K_pphi
    fim[idx_phi, idx_p] <- K_pphi
  }
  
  if (has_ma) {
    fim[idx_t, idx_t] <- K_tt
    fim[idx_t, idx_a] <- K_ta
    fim[idx_a, idx_t] <- K_ta
    fim[idx_t, idx_phi] <- K_tphi
    fim[idx_phi, idx_t] <- K_tphi
  }
  
  if (has_ar && has_ma) {
    fim[idx_p, idx_t] <- K_pt
    fim[idx_t, idx_p] <- t(K_pt)
  }
  
  # Name the rows and columns
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta  <- if (has_ma) paste0("theta", ma_lags) else character(0)
  
  names_fim <- c("alpha", names_varphi, names_theta, "phi")
  colnames(fim) <- names_fim
  rownames(fim) <- names_fim
  
  # Prepare fitted values as a time series object
  fitted_ts <- ts(c(rep(NA, max_lag), mu_effective),
                  start = start(y),
                  frequency = frequency(y))
  
  # Prepare output list
  output_list <- list(
    fisher_info_mat = fim,
    fitted_ts = fitted_ts,
    muhat_effective = mu_effective,
    etahat_full = eta,
    errorhat_full = error
  )
  
  return(output_list)
}