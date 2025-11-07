#' Fit Beta Autoregressive Moving Average (βARMA) Models
#'
#' @description
#' Fits Beta Autoregressive Moving Average (βARMA) models for time series 
#' valued in (0, 1) via standard Maximum Likelihood Estimation (MLE).
#'
#' @details
#' This function fits the βARMA(p,q) model as proposed by Rocha & 
#' Cribari-Neto (2009, with erratum 2017). It serves as the main wrapper 
#' for the optimization process, calling helper functions for the 
#' log-likelihood, score vector, and Fisher Information Matrix.
#'
#' The model is specified via the `ar` and `ma` arguments:
#' \itemize{
#'   \item βARMA(p,q): `ar` and `ma` are specified.
#'   \item βAR(p): Specify `ar` and set `ma = NA`.
#'   \item βMA(q): Specify `ma` and set `ar = NA`.
#' }
#'
#' The optimization is performed using the "BFGS" algorithm via the
#' `optim` function.
#'
#' @author
#' Everton da Costa (everton.ecosta@ufpe.br);
#' Francisco Cribari-Neto (francisco.cribari@ufpe.br);
#'
#' @note
#' The original version of this function was developed by Fabio M. Bayer
#'  (bayer@ufsm.br). Substantially modified and improved by Everton da Costa
#'  (everton.ecosta@ufpe.br). With suggestions and contributions from
#'  Francisco Cribari Neto (francisco.cribari@ufpe.br).
#'
#' @references
#' Rocha, A.V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. *TEST*, 18(3), 529-545.
#' <doi:10.1007/s11749-008-0113-2>
#'
#' Rocha, A.V., & Cribari-Neto, F. (2017). Erratum to: Beta autoregressive
#' moving average models. *TEST*, 26, 451–459.
#' <doi:10.1007/s11749-017-0528-4>
#'
#' @param y A time series object (`ts`) with values in (0, 1).
#' @param ar A numeric vector of autoregressive (AR) lags (e.g., `c(1, 2)`).
#' @param ma A numeric vector of moving average (MA) lags (e.g., `1`).
#' @param link The link function to connect the mean to the linear
#'   predictor. Default is `"logit"`.
#'
#' @importFrom stats is.ts pnorm optim
#'
#' @return
#' A list containing the fitted model results:
#' \item{model}{A summary table of coefficients, standard errors, z-values,
#'   and p-values.}
#' \item{coef}{A named vector of all estimated coefficients.}
#' \item{vcov}{The variance-covariance matrix of the coefficients.}
#' \item{fitted}{The in-sample fitted mean values (`muhat`).}
#' \item{loglik}{The value of the log-likelihood function at the optimum.}
#' \item{aic, bic, hq}{Akaike, Bayesian, and Hannan-Quinn information
#'   criteria.}
#' \item{conv}{Convergence code from `optim` (`0` indicates success).}
#' \item{phi, alpha, varphi, theta}{Individual estimated parameters.}
#' \item{etahat}{The estimated linear predictor values.}
#' \item{muhat}{The estimated mean values (same as `fitted`).}
#' \item{errorhat}{The estimated errors on the predictor scale.}
#' \item{fisher_info_mat}{The observed Fisher Information Matrix.}
#' \item{start_values}{Initial values used in the optimization.}
#' \item{status_inv_inf_matrix}{Indicates if the Fisher matrix was
#'   invertible (`0` = success, `1` = failure).}
#' \item{opt}{The raw output object from the `optim` call.}
#' @export
barma <- function(y,
                  ar = NA, ma = NA,
                  link = "logit") {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Validate Input Data ---
  # ------------------------------------------------------------------------- #
  if (min(y) <= 0 || max(y) >= 1) stop("y values must be in (0,1)!")
  if (is.ts(y) == FALSE) stop("Data must be a time-series (ts) object")
  
  # --- Store the function call ---
  z <- list(call = match.call())
  
  # ------------------------------------------------------------------------- #
  # --- 2. Determine Model Structure ---
  # ------------------------------------------------------------------------- #
  has_ar <- !is.null(ar) && !any(is.na(ar)) && length(ar) > 0
  has_ma <- !is.null(ma) && !any(is.na(ma)) && length(ma) > 0
  
  # Use integer(0) for empty lags, as expected by helper functions
  ar_lags <- if (has_ar) ar else integer(0)
  ma_lags <- if (has_ma) ma else integer(0)
  
  n_ar_params <- length(ar_lags)
  n_ma_params <- length(ma_lags)
  
  # Create parameter names for the output
  names_varphi <- if (has_ar) paste0("varphi", ar_lags) else character(0)
  names_theta  <- if (has_ma) paste0("theta", ma_lags) else character(0)
  
  # ------------------------------------------------------------------------- #
  # --- 3. Setup Time Series Properties ---
  # ------------------------------------------------------------------------- #
  ar_order <- if (has_ar) max(ar_lags) else 0L
  ma_order <- if (has_ma) max(ma_lags) else 0L
  max_lag  <- max(ar_order, ma_order)
  
  n_obs <- length(y)
  
  if (n_obs <= max_lag) {
    stop("Insufficient observations for the specified lag structure")
  }
  
  # ------------------------------------------------------------------------- #
  # --- 4. Setup Link Function & Transformed Data ---
  # ------------------------------------------------------------------------- #
  # Using make_link_structure() available in the package environment
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  ynew <- linkfun(y)
  
  # ------------------------------------------------------------------------- #
  # --- 5. Setup Parameter Indices ---
  # ------------------------------------------------------------------------- #
  idx_alpha <- 1
  idx_varphi <- if (has_ar) {
    2:(idx_alpha + n_ar_params)
  } else {
    integer(0)
  }
  
  idx_theta <- if (has_ma) {
    idx_theta_start <- max(c(idx_alpha, idx_varphi)) + 1
    idx_theta_end <- max(c(idx_alpha, idx_varphi)) + n_ma_params
    
  idx_theta <- idx_theta_start:idx_theta_end
  } else {
    integer(0)
  }
  
  idx_phi <- max(c(idx_alpha, idx_varphi, idx_theta)) + 1
  n_params <- idx_phi
  
  # ------------------------------------------------------------------------- #
  # --- 6. Get Initial Values ---
  # ------------------------------------------------------------------------- #
  # start_values(): helper function
  start_values <- start_values(y, link = link, ar = ar_lags, ma = ma_lags)
  
  # Ensure start_values has the correct length
  if (length(start_values) != n_params) {
    warning("Length of start_values does not match number of parameters.")
  }
  
  
  # ------------------------------------------------------------------------- #
  # --- 7. Optimization (MLE) ---
  # ------------------------------------------------------------------------- #
  opt <- optim(
    par = start_values,
    fn = function(x) {
      (-1) * loglik_barma(
        y = y,
        ar = ar_lags,
        ma = ma_lags,
        alpha = x[idx_alpha],
        varphi = x[idx_varphi],
        theta = x[idx_theta],
        phi = x[idx_phi],
        link = link
      )
    },
    gr = function(x) {
      (-1) * score_vector_barma(
        y = y,
        ar = ar_lags,
        ma = ma_lags,
        alpha = x[idx_alpha],
        varphi = x[idx_varphi],
        theta = x[idx_theta],
        phi = x[idx_phi],
        link = link
      )
    },
    method = "BFGS"
  )
  
  if (opt$conv != 0) {
    warning("FUNCTION DID NOT CONVERGE!")
  }
  
  z$conv <- opt$convergence
  z$opt <- opt
  z$loglik <- -1 * opt$value
  
  # ------------------------------------------------------------------------- #
  # --- 8. Extract Final Parameters ---
  # ------------------------------------------------------------------------- #
  coef <- opt$par[1:n_params]
  names(coef) <- c("alpha", names_varphi, names_theta, "phi")
  z$coef <- coef
  
  # Extract by name
  z$alpha  <- coef["alpha"]
  z$varphi <- if (has_ar) coef[names_varphi] else numeric(0)
  z$theta  <- if (has_ma) coef[names_theta] else numeric(0)
  z$phi    <- coef["phi"]
  
  # ------------------------------------------------------------------------- #
  # --- 9. Compute FIM and Fitted Values (Efficiently) ---
  # ------------------------------------------------------------------------- #
  
  # This one call gets the FIM and all fitted components
  fim_results <- fim_barma(
    y = y,
    ar = ar_lags,
    ma = ma_lags,
    alpha = z$alpha,
    varphi = z$varphi,
    theta = z$theta,
    phi = z$phi,
    link = link
  )
  
  # Unpack the results
  z$fisher_info_mat <- fim_results$fisher_info_mat
  z$fitted   <- fim_results$fitted_ts      # 'ts' object, NA-padded
  z$muhat    <- fim_results$fitted_ts      # Alias for fitted      
  z$etahat   <- fim_results$etahat_full    # Full vector, NA-padded    
  z$errorhat <- fim_results$errorhat_full  # Full vector, 0-padded
  
  # ------------------------------------------------------------------------- #
  # --- 10. Store Raw Components for Summary ---
  # ------------------------------------------------------------------------- #
  z$y <- y
  z$link <- link
  z$start_values <- start_values
  z$n_params <- n_params
  z$n_obs <- n_obs
  z$max_lag <- max_lag
  z$ar_lags <- ar_lags # Stored for residuals()
  z$ma_lags <- ma_lags # Stored for residuals()
  
  # --- Assign class and return ---
  class(z) <- "barma"
  
  return(z)
}