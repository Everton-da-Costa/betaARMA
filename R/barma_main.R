#' Fit Beta Autoregressive Moving Average (βARMA) Models
#'
#' @description
#' Fits a comprehensive suite of Beta Autoregressive Moving Average (βARMA)
#' models for time series valued in (0, 1). Its core feature is the
#' integration of a Ridge (L2) penalty to provide numerically stable
#' parameter estimates, addressing common convergence issues. The function
#' supports both penalized (PMLE) and standard (MLE) estimation.
#'
#' @details
#' The estimation of βARMA models can be susceptible to numerical
#' instability, which may lead to convergence failures or implausible
#' estimates. To overcome this, this function implements a penalized
#' maximum likelihood approach using a Ridge (L2) penalty, as introduced
#' by Cribari-Neto, Costa, & Fonseca (2025).
#'
#' By adding a penalty term to the log-likelihood, the function's
#' curvature is improved, leading to more reliable and robust
#' optimization. Standard MLE, as in Rocha & Cribari-Neto (2009), is
#' performed when `penalty = 0`.
#'
#' The model is specified via the `ar`, `ma`, and `X` arguments:
#' \itemize{
#'   \item βARMA(p,q): `ar` and `ma` are specified.
#'   \item βAR(p): Only `ar` is specified.
#'   \item βMA(q): Only `ma` is specified.
#'   \item βARMAX, βARX, βMAX: Models with exogenous variables `X`.
#' }
#'
#' The optimization is performed using the "BFGS" algorithm via the
#' `optim` function.
#'
#' \strong{Implementation Notes:}
#' This version is a substantial rewrite of the original code by Fabio M.
#' Bayer. The modifications, implemented by Everton da Costa, include:
#' \itemize{
#'   \item A complete rewrite of the code following the tidyverse style guide.
#'   \item A thorough review against the original article by Rocha &
#'     Cribari-Neto (2009, with erratum 2017).
#'   \item Code optimization to avoid redundant computations.
#'   \item A modular structure to allow functions to be run separately.
#' }
#'
#' @author
#' Everton da Costa (everton.ecosta@ufpe.br);
#' Francisco Cribari-Neto (francisco.cribari@ufpe.br);
#' Rodney V. Fonseca (rodneyfonseca@ufba.br).
#'
#' @note
#' The original version of this function was developed by Fabio M. Bayer
#'  (bayer@ufsm.br). Substantially modified and improved by Everton da Costa
#'  (everton.ecosta@ufpe.br). With suggestions and contributions from
#'  Rodney Fonceca (rodneyfonseca@ufba.br) and
#'  Francisco Cribari Neto (francisco.cribari@ufpe.br).
#'
#' @references
#' Cribari-Neto, F., Costa, E., & Fonseca, R. V. (2025). Numerical
#' stability enhancements in beta autoregressive moving average model
#' estimation. *Brazilian Journal of Probability and Statistics*,
#' Vol(Issue), pp-pp.
#'
#' Rocha, A.V., & Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. *TEST*, 18(3), 529-545.
#' <doi:10.1007/s11749-008-0113-2>
#'
#' @param y A time series object (`ts`) with values in (0, 1).
#' @param ar A numeric vector of autoregressive (AR) lags (e.g., `c(1, 2)`).
#' @param ma A numeric vector of moving average (MA) lags (e.g., `1`).
#' @param link The link function to connect the mean to the linear
#'   predictor. Default is `"logit"`.
#' @param penalty A non-negative scalar for the L2 (Ridge) penalty. Set to a
#'   positive value to enable penalized estimation, which is recommended
#'   for addressing numerical instability. Default is `0` (standard MLE).
#' @param h1 An integer specifying the forecast horizon. Default is `6`.
#' @param X An optional matrix of exogenous variables.
#' @param X_hat An optional matrix of future exogenous variables, required
#'   for forecasting when `X` is used.
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
#' \item{forecast}{A vector of out-of-sample forecasted values.}
#' \item{loglik}{The value of the log-likelihood function at the optimum.}
#' \item{aic, bic, hq}{Akaike, Bayesian, and Hannan-Quinn information
#'   criteria.}
#' \item{conv}{Convergence code from `optim` (`0` indicates success).}
#' \item{resid2}{Standardized residuals on the predictor scale.}
#' \item{phi, alpha, varphi, theta, beta}{Individual estimated parameters.}
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
                  link = "logit",
                  penalty = 0,
                  h1 = 6,
                  X = NA, X_hat = NA) {
  # Link functions
  link_structure <- make_link_structure(link)

  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  # Check-up time series
  if (min(y) <= 0 || max(y) >= 1) stop("OUT OF RANGE (0,1)!")
  if (is.ts(y) == FALSE) stop("Data might not be a time-series object")

  # Create variable for the AR and MA components
  has_ar <- !any(is.na(ar))
  has_ma <- !any(is.na(ma))
  has_X <- !any(is.na(X))
  use_penalty <- penalty != 0

  n <- length(y)

  # Define p, q (max lags) and p1, q1 (number of parameters)
  p <- ifelse(has_ar, max(ar), 0)
  q <- ifelse(has_ma, max(ma), 0)
  p1 <- ifelse(has_ar, length(ar), 0)
  q1 <- ifelse(has_ma, length(ma), 0)

  m <- max(p, q)
  ynew <- linkfun(y)

  names_beta <- colnames(X)
  y_prev <- c(rep(NA, (n + h1)))

  # Create variable names for the AR, MA and X components conditionally
  names_varphi <- if (has_ar) paste("varphi", ar, sep = "") else character(0)
  names_theta <- if (has_ma) paste("theta", ma, sep = "") else character(0)
  names_beta <- if (has_X) colnames(X) else character(0)

  # ========================================================================= #
  # BARMA model
  if (has_ar && has_ma && !has_X && !use_penalty) {
    # cat("Fitting βARMA model (MLE) \n")

    z <- list()

    m <- max(p, q, na.rm = T)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization using optim
    # ----------------------------------------------------------------------- #
    opt_arma <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_arma(
          y = y,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_arma(
          y = y,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          link = link
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_arma$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # convergence status
    z$conv <- opt_arma$convergence

    # estimates
    coef_arma <- opt_arma$par[1:(p1 + q1 + 2)]
    names(coef_arma) <- c("alpha", names_varphi, names_theta, "phi")

    # output estimates
    z$coef <- coef_arma

    # log likelihood value
    z$loglik <- -1 * opt_arma$value

    z$opt <- opt_arma

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_arma[1]
    varphi <- coef_arma[2:(p1 + 1)]
    theta <- coef_arma[(p1 + 2):(p1 + q1 + 1)]
    phi <- coef_arma[p1 + q1 + 2]

    z$alpha <- alpha
    z$varphi <- varphi
    z$theta <- theta
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARMA
    # ------------------------------------------------------------------------
    output_inf_matrix_arma <- inf_matrix_arma(
      y = y,
      ar = ar,
      ma = ma,
      alpha = z$alpha,
      varphi = z$varphi,
      theta = z$theta,
      phi = z$phi,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_arma$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arma$muhat
    z$fitted <- output_inf_matrix_arma$fitted
    z$etahat <- output_inf_matrix_arma$etahat
    z$errorhat <- output_inf_matrix_arma$errorhat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha +
        (varphi %*% ynew_prev[n + i - ar]) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      errorhat[n + i] <- 0
    }

  }

  # BARMA model - Ridge Regression Penalty
  if (has_ar && has_ma && !has_X && use_penalty) {
    # cat("Fitting βARMA model with Ridge Penalty (PMLE) \n")

    z <- list()
    m <- max(p, q, na.rm = T)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_arma_ridge <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_arma_ridge(
          y = y,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_arma_ridge(
          y = y,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          link = link,
          penalty = penalty
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_arma_ridge$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # optimization values
    z$opt <- opt_arma_ridge

    # convergence status
    z$conv <- opt_arma_ridge$convergence

    # estimates
    coef_arma_ridge <- opt_arma_ridge$par[1:(p1 + q1 + 2)]
    names(coef_arma_ridge) <- c("alpha", names_varphi, names_theta, "phi")

    # output estimates
    z$coef <- coef_arma_ridge

    # log likelihood value
    z$loglik <- -1 * opt_arma_ridge$value

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha <- coef_arma_ridge[1]
    varphi <- coef_arma_ridge[2:(p1 + 1)]
    theta <- coef_arma_ridge[(p1 + 2):(p1 + q1 + 1)]
    phi <- coef_arma_ridge[p1 + q1 + 2]

    # --------------------------------------------------------------------- #
    z$alpha <- alpha
    z$varphi <- varphi
    z$theta <- theta
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARMA
    # ------------------------------------------------------------------------
    output_inf_matrix_arma_ridge <- inf_matrix_arma_ridge(
      y = y,
      ar = ar,
      ma = ma,
      alpha = z$alpha,
      varphi = z$varphi,
      theta = z$theta,
      phi = z$phi,
      penalty = penalty,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_arma_ridge$fisher_info_mat_penalized

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arma_ridge$muhat
    z$fitted <- output_inf_matrix_arma_ridge$fitted
    z$etahat <- output_inf_matrix_arma_ridge$etahat
    z$errorhat <- output_inf_matrix_arma_ridge$errorhat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha +
        (varphi %*% ynew_prev[n + i - ar]) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      errorhat[n + i] <- 0
    }

  }

  # BARMAX model
  if (has_ar && has_ma && has_X && !use_penalty) {
    # cat("Fitting βARMAX model (MLE) \n")

    # -------------------------------------------------------------------------
    z <- list()
    z$model_name <- "MLE"

    m <- max(p, q, na.rm = TRUE)

    X <- as.matrix(X)
    X_hat <- as.matrix(X_hat)
    z$n_obs <- n

    # -------------------------------------------------------------------------
    # initial values
    # -------------------------------------------------------------------------
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    names_par <- c("alpha", names_varphi, names_theta, "phi", names_beta)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_armax(
          y = y,
          X = X,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          beta = x[(p1 + q1 + 3):length(x)],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_armax(
          y = y,
          X = X,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          beta = x[(p1 + q1 + 3):length(x)],
          link = link
        )
      },
      method = "BFGS"
    )

    if (opt$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # ----------------------------------------------------------------------- #

    # convergence status
    z$conv <- opt$convergence
    z$opt <- opt

    # log likelihood value
    z$loglik <- -1 * opt$value

    # ----------------------------------------------------------------------- #
    names_par <- c("alpha", names_varphi, names_theta, "phi", names_beta)

    # estimates
    coef <- opt$par
    names(coef) <- names_par
    z$coef <- coef

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha <- coef[1]
    varphi <- coef[2:(p1 + 1)]
    theta <- coef[(p1 + 2):(p1 + q1 + 1)]
    phi <- coef[p1 + q1 + 2]
    beta <- coef[(p1 + q1 + 3):length(coef)]

    # --------------------------------------------------------------------- #
    z$alpha <- alpha
    z$varphi <- varphi
    z$theta <- theta
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARMA
    # ------------------------------------------------------------------------
    output_inf_matrix_arma <- inf_matrix_armax(
      y = y,
      ar = ar,
      ma = ma,
      X = X,
      alpha = z$alpha,
      varphi = z$varphi,
      theta = z$theta,
      phi = z$phi,
      beta = z$beta,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_arma$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arma$muhat
    z$fitted <- output_inf_matrix_arma$fitted
    z$etahat <- output_inf_matrix_arma$etahat
    z$errorhat <- output_inf_matrix_arma$errorhat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat
    ynew_prev <- c(ynew, rep(NA, h1))

    X_prev <- rbind(X, X_hat)

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha +
        X_prev[n + i, ] %*% as.matrix(beta) +
        (varphi %*% (ynew_prev[n + i - ar] - X_prev[n + i - ar, ] %*% as.matrix(beta))) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0
    }

  }

  # BARMAX Ridge Regression model
  if (has_ar && has_ma && has_X && use_penalty) {
    # cat("Fitting βARMAX model with Ridge Penalty (PMLE) \n")

    # -------------------------------------------------------------------------
    z <- list()
    z$model_name <- "PMLE"

    m <- max(p, q, na.rm = TRUE)

    X <- as.matrix(X)
    X_hat <- as.matrix(X_hat)
    z$n_obs <- n

    # -------------------------------------------------------------------------
    # initial values
    # -------------------------------------------------------------------------
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    names_par <- c("alpha", names_varphi, names_theta, "phi", names_beta)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_armax_ridge(
          y = y,
          X = X,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          beta = x[(p1 + q1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_armax_ridge(
          y = y,
          X = X,
          ar = ar,
          ma = ma,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          theta = x[(p1 + 2):(p1 + q1 + 1)],
          phi = x[p1 + q1 + 2],
          beta = x[(p1 + q1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },
      method = "BFGS"
    )

    if (opt$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # ----------------------------------------------------------------------- #

    # convergence status
    z$conv <- opt$convergence
    z$opt <- opt

    # log likelihood value
    z$loglik <- -1 * opt$value

    # ----------------------------------------------------------------------- #
    names_par <- c("alpha", names_varphi, names_theta, "phi", names_beta)

    # estimates
    coef <- opt$par
    names(coef) <- names_par
    z$coef <- coef

    alpha <- coef[1]
    varphi <- coef[2:(p1 + 1)]
    theta <- coef[(p1 + 2):(p1 + q1 + 1)]
    phi <- coef[p1 + q1 + 2]
    beta <- coef[(p1 + q1 + 3):length(coef)]

    z$alpha <- alpha
    z$varphi <- varphi
    z$theta <- theta
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha <- coef[1]
    varphi <- coef[2:(p1 + 1)]
    theta <- coef[(p1 + 2):(p1 + q1 + 1)]
    phi <- coef[p1 + q1 + 2]
    beta <- coef[(p1 + q1 + 3):length(coef)]

    # --------------------------------------------------------------------- #
    z$alpha <- alpha
    z$varphi <- varphi
    z$theta <- theta
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARMA
    # ------------------------------------------------------------------------
    output_inf_matrix_arma_ridge <- inf_matrix_armax_ridge(
      y = y,
      ar = ar,
      ma = ma,
      X = X,
      alpha = z$alpha,
      varphi = z$varphi,
      theta = z$theta,
      phi = z$phi,
      beta = z$beta,
      link = link,
      penalty = penalty
    )

    fisher_info_mat <- output_inf_matrix_arma_ridge$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arma_ridge$muhat
    z$fitted <- output_inf_matrix_arma_ridge$fitted
    z$etahat <- output_inf_matrix_arma_ridge$etahat
    z$errorhat <- output_inf_matrix_arma_ridge$errorhat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat
    ynew_prev <- c(ynew, rep(NA, h1))

    X_prev <- rbind(X, X_hat)

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha +
        X_prev[n + i, ] %*% as.matrix(beta) +
        (varphi %*% (ynew_prev[n + i - ar] - X_prev[n + i - ar, ] %*% as.matrix(beta))) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0
    }

  }

  # ============================================================================
  # BAR Model
  if (has_ar && !has_ma && !has_X && !use_penalty) {
    # cat("Fitting βAR model (MLE) \n")

    z <- list()
    m <- max(p, na.rm = TRUE)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization using optim
    # ----------------------------------------------------------------------- #
    opt_ar <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_ar(
          y = y,
          ar = ar,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_ar(
          y = y,
          ar = ar,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          link = link
        )
      },

      # control = list(
      #   maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_ar$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$opt <- opt_ar

    # convergence status
    z$conv <- opt_ar$convergence

    # estimates
    coef_ar <- opt_ar$par
    names(coef_ar) <- c("alpha", names_varphi, "phi")

    # output estimates
    z$coef <- coef_ar

    # log likelihood value
    z$loglik <- -1 * opt_ar$value

    z$opt <- opt_ar

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_ar[1]
    varphi <- coef_ar[2:(p1 + 1)]
    phi <- coef_ar[p1 + 2]

    z$alpha <- alpha
    z$varphi <- varphi
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, AR
    # ------------------------------------------------------------------------
    output_inf_matrix_ar <- inf_matrix_ar(
      y = y,
      ar = ar,
      alpha = z$alpha,
      varphi = z$varphi,
      phi = z$phi,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_ar$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_ar$muhat
    z$fitted <- output_inf_matrix_ar$fitted
    z$etahat <- output_inf_matrix_ar$etahat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + varphi %*% ynew_prev[n + i - ar]
      y_prev[n + i] <- linkinv(ynew_prev[n + i])
    }

  }

  # BAR model - Ridge Regression Penalty
  if (has_ar && !has_ma && !has_X && use_penalty) {
    # cat("Fitting βAR model with Ridge Penalty (PMLE) \n")

    z <- list()
    m <- max(p, na.rm = T)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_ar_ridge <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_ar_ridge(
          y = y,
          ar = ar,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_ar_ridge(
          y = y,
          ar = ar,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          link = link,
          penalty = penalty
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_ar_ridge$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # optimization values
    z$opt <- opt_ar_ridge

    # convergence status
    z$conv <- opt_ar_ridge$convergence

    # estimates
    coef_ar_ridge <- opt_ar_ridge$par[1:(p1 + 2)]
    names(coef_ar_ridge) <- c("alpha", names_varphi, "phi")

    # output estimates
    z$coef <- coef_ar_ridge

    # log likelihood value
    z$loglik <- -1 * opt_ar_ridge$value

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha <- coef_ar_ridge[1]
    varphi <- coef_ar_ridge[2:(p1 + 1)]
    phi <- coef_ar_ridge[p1 + 2]

    # --------------------------------------------------------------------- #
    z$alpha <- alpha
    z$varphi <- varphi
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ar
    # ------------------------------------------------------------------------
    output_inf_matrix_ar_ridge <- inf_matrix_ar_ridge(
      y = y,
      ar = ar,
      alpha = z$alpha,
      varphi = z$varphi,
      phi = z$phi,
      penalty = penalty,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_ar_ridge$fisher_info_mat_penalized

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_ar_ridge$muhat
    z$fitted <- output_inf_matrix_ar_ridge$fitted
    z$etahat <- output_inf_matrix_ar_ridge$etahat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + varphi %*% ynew_prev[n + i - ar]
      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      # errorhat[n + i] <- 0
    }

  }

  # BARX Model
  if (has_ar && !has_ma && has_X && !use_penalty) {
    # cat("Fitting βARX model (MLE) \n")

    z <- list()
    m <- max(p, na.rm = TRUE)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization using optim
    # ----------------------------------------------------------------------- #
    opt_arx <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_arx(
          y = y,
          ar = ar,
          X = X,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          beta = x[(p1 + 3):length(x)],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_arx(
          y = y,
          ar = ar,
          X = X,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          beta = x[(p1 + 3):length(x)],
          link = link
        )
      },

      # control = list(
      #   maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_arx$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$opt <- opt_arx

    # convergence status
    z$conv <- opt_arx$convergence

    # estimates
    coef_arx <- opt_arx$par
    names(coef_arx) <- c("alpha", names_varphi, "phi", names_beta)

    # output estimates
    z$coef <- coef_arx

    # log likelihood value
    z$loglik <- -1 * opt_arx$value

    z$opt <- opt_arx

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_arx[1]
    varphi <- coef_arx[2:(p1 + 1)]
    phi <- coef_arx[p1 + 2]
    beta <- coef_arx[(p1 + 3):length(coef_arx)]

    z$alpha <- alpha
    z$varphi <- varphi
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, AR
    # ------------------------------------------------------------------------
    output_inf_matrix_arx <- inf_matrix_arx(
      y = y,
      ar = ar,
      X = X,
      alpha = z$alpha,
      varphi = z$varphi,
      phi = z$phi,
      beta = z$beta,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_arx$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arx$muhat
    z$fitted <- output_inf_matrix_arx$fitted
    z$etahat <- output_inf_matrix_arx$etahat

    # --------------------------------------------------------------------- #
    # Forecasting
    # --------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    X_prev <- rbind(X, X_hat)

    # ----------------------------------------------------------------------- #
    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + X_prev[n + i, ] %*%
        as.matrix(beta) +
        (varphi %*% (ynew_prev[n + i - ar] - X_prev[n + i - ar, ] %*%
          as.matrix(beta)))

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      # errorhat[n + i] <- 0
    }

  }

  # BARX Model - Ridge Regression Penalty
  if (has_ar && !has_ma && has_X && use_penalty) {
    # cat("Fitting βARX model with Ridge Penalty (PMLE) \n")

    z <- list()
    m <- max(p, na.rm = TRUE)

    X <- as.matrix(X)
    X_hat <- as.matrix(X_hat)
    z$n_obs <- n

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    names_par <- c("alpha", names_varphi, "phi", names_beta)
    names(start_values) <- names_par

    # ----------------------------------------------------------------------- #
    # optimization using optim
    # ----------------------------------------------------------------------- #
    opt_arx_ridge <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_arx_ridge(
          y = y,
          ar = ar,
          X = X,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          beta = x[(p1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_arx_ridge(
          y = y,
          ar = ar,
          X = X,
          alpha = x[1],
          varphi = x[2:(p1 + 1)],
          phi = x[p1 + 2],
          beta = x[(p1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },

      # control = list(
      #   maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_arx_ridge$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$opt <- opt_arx_ridge

    # convergence status
    z$conv <- opt_arx_ridge$convergence

    # estimates
    coef_arx_ridge <- opt_arx_ridge$par
    names(coef_arx_ridge) <- c("alpha", names_varphi, "phi", names_beta)

    # output estimates
    z$coef <- coef_arx_ridge

    # log likelihood value
    z$loglik <- -1 * opt_arx_ridge$value

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_arx_ridge[1]
    varphi <- coef_arx_ridge[2:(p1 + 1)]
    phi <- coef_arx_ridge[p1 + 2]
    beta <- coef_arx_ridge[(p1 + 3):length(coef_arx_ridge)]

    z$alpha <- alpha
    z$varphi <- varphi
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARX with Penalty
    # ------------------------------------------------------------------------
    output_inf_matrix_arx_ridge <- inf_matrix_arx_ridge(
      y = y,
      ar = ar,
      X = X,
      alpha = z$alpha,
      varphi = z$varphi,
      phi = z$phi,
      beta = z$beta,
      link = link,
      penalty = penalty
    )

    fisher_info_mat <- output_inf_matrix_arx_ridge$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arx_ridge$muhat
    z$fitted <- output_inf_matrix_arx_ridge$fitted
    z$etahat <- output_inf_matrix_arx_ridge$etahat

    # --------------------------------------------------------------------- #
    # Forecasting
    # --------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    X_prev <- rbind(X, X_hat)

    # ----------------------------------------------------------------------- #
    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + X_prev[n + i, ] %*%
        as.matrix(beta) +
        (varphi %*% (ynew_prev[n + i - ar] - X_prev[n + i - ar, ] %*%
                       as.matrix(beta)))

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      # errorhat[n + i] <- 0
    }
  }

  # ============================================================================
  # BMA model
  if (!has_ar && has_ma && !has_X && !use_penalty) {
    # cat("Fitting βMA model (MLE) \n")
    z <- list()
    m <- max(q, na.rm = TRUE)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    names_par_ma <- c("alpha", names_theta, "phi")

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_ma <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_ma(
          y = y,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_ma(
          y = y,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          link = link
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_ma$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # convergence status
    z$conv <- opt_ma$convergence

    # estimates
    coef_ma <- opt_ma$par
    names(coef_ma) <- c("alpha", names_theta, "phi")

    # output estimates
    z$coef <- coef_ma

    # log likelihood value
    z$loglik <- -1 * opt_ma$value

    z$opt <- opt_ma

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_ma[1]
    theta <- coef_ma[2:(q1 + 1)]
    phi <- coef_ma[q1 + 2]

    z$alpha <- alpha
    z$theta <- theta
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ma
    # ------------------------------------------------------------------------
    output_inf_matrix_ma <- inf_matrix_ma(
      y = y,
      ma = ma,
      alpha = z$alpha,
      theta = z$theta,
      phi = z$phi,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_ma$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_ma$muhat
    z$fitted <- output_inf_matrix_ma$fitted
    z$etahat <- output_inf_matrix_ma$etahat
    z$errorhat <- output_inf_matrix_ma$errorhat

    # --------------------------------------------------------------------- #
    # Forecasting
    # --------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    errorhat <- z$errorhat
    y_prev[1:n] <- z$fitted

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + theta %*% errorhat[n + i - ma]
      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0
    }

  }

  # BMA model - Ridge Regression Penalty
  if (!has_ar && has_ma && !has_X && use_penalty) {
    # cat("Fitting βMA model with Ridge Penalty (PMLE) \n")

    z <- list()
    m <- max(q, na.rm = T)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_ma_ridge <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_ma_ridge(
          y = y,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_ma_ridge(
          y = y,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          link = link,
          penalty = penalty
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_ma_ridge$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # optimization values
    z$opt <- opt_ma_ridge

    # convergence status
    z$conv <- opt_ma_ridge$convergence

    # estimates
    coef_ma_ridge <- opt_ma_ridge$par[1:(q1 + 2)]
    names(coef_ma_ridge) <- c("alpha", names_theta, "phi")

    # output estimates
    z$coef <- coef_ma_ridge

    # log likelihood value
    z$loglik <- -1 * opt_ma_ridge$value

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha <- coef_ma_ridge[1]
    theta <- coef_ma_ridge[2:(q1 + 1)]
    phi <- coef_ma_ridge[q1 + 2]

    z$alpha <- alpha
    z$theta <- theta
    z$phi <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ma
    # ------------------------------------------------------------------------
    output_inf_matrix_ma_ridge <- inf_matrix_ma_ridge(
      y = y,
      ma = ma,
      alpha = z$alpha,
      theta = z$theta,
      phi = z$phi,
      penalty = penalty,
      link = link
    )

    fisher_info_mat <- output_inf_matrix_ma_ridge$fisher_info_mat_penalized

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_ma_ridge$muhat
    z$fitted <- output_inf_matrix_ma_ridge$fitted
    z$etahat <- output_inf_matrix_ma_ridge$etahat
    z$errorhat <- output_inf_matrix_ma_ridge$errorhat

    # --------------------------------------------------------------------- #
    # Forecasting
    # --------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + theta %*% errorhat[n + i - ma]
      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0
    }

  }

  # BMAX model
  if (!has_ar && has_ma && has_X && !use_penalty) {
    # cat("Fitting βMAX model (MLE) \n")
    z <- list()
    m <- max(q, na.rm = TRUE)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)
    names_par <- c("alpha", names_theta, "phi", names_beta)
    names(start_values) <- names_par

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_max <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_max(
          y = y,
          X = X,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          beta = x[(q1 + 3):length(x)],
          link = link
        )
      },
      gr = function(x) {
        (-1) * score_vector_max(
          y = y,
          X = X,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          beta = x[(q1 + 3):length(x)],
          link = link
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_max$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$conv <- opt_max$convergence
    z$opt <- opt_max
    z$loglik <- -1 * opt_max$value

    # Estimates
    coef_max <- opt_max$par
    names(coef_max) <- names_par
    z$coef <- coef_max

    # Parameter extraction
    alpha <- coef_max[1]
    theta <- coef_max[2:(q1 + 1)]
    phi <- coef_max[q1 + 2]
    beta <- coef_max[(q1 + 3):length(coef_max)]

    z$alpha <- alpha
    z$theta <- theta
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ma
    # ------------------------------------------------------------------------
    output_inf_matrix_max <- inf_matrix_max(
      y = y,
      ma = ma,
      X = X,
      alpha = z$alpha,
      theta = z$theta,
      phi = z$phi,
      beta = z$beta,
      link = link
    )

    z$fisher_info_mat <- output_inf_matrix_max$fisher_info_mat
    z$muhat <- output_inf_matrix_max$muhat
    z$fitted <- output_inf_matrix_max$fitted
    z$etahat <- output_inf_matrix_max$etahat
    z$errorhat <- output_inf_matrix_max$errorhat

    # ----------------------------------------------------------------------- #
    # Forecasting
    # ----------------------------------------------------------------------- #
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat
    ynew_prev <- c(ynew, rep(NA, h1))
    X_prev <- rbind(X, X_hat)

    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha + X_prev[n + i, ] %*% as.matrix(beta) +
        theta %*% errorhat[n + i - ma]
      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0
    }

  }

  # BMAX model - Ridge Regression Penalty
  if (!has_ar && has_ma && has_X && use_penalty) {
    # cat("Fitting βMAX model with Ridge Penalty (PMLE) \n")
    z <- list()
    m <- max(q, na.rm = TRUE)

    # initial values
    start_values <- start_values(y, link = link, ar = ar, ma = ma, X = X)

    # ----------------------------------------------------------------------- #
    # optimization: optim
    # ----------------------------------------------------------------------- #
    opt_max <- optim(
      par = start_values,
      fn = function(x) {
        (-1) * loglik_max_ridge(
          y = y,
          X = X,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          beta = x[(q1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },
      gr = function(x) {
        (-1) * score_vector_max_ridge(
          y = y,
          X = X,
          ma = ma,
          alpha = x[1],
          theta = x[2:(q1 + 1)],
          phi = x[q1 + 2],
          beta = x[(q1 + 3):length(x)],
          link = link,
          penalty = penalty
        )
      },

      # control = list(
      #   # maxit = 1000,
      #   reltol = 1e-12
      # ),

      method = "BFGS"
    )

    if (opt_max$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # convergence status
    z$conv <- opt_max$convergence

    # estimates
    coef_max <- opt_max$par
    names(coef_max) <- c("alpha", names_theta, "phi", names_beta)

    # output estimates
    z$coef <- coef_max

    # log likelihood value
    z$loglik <- -1 * opt_max$value

    z$opt <- opt_max

    # estimates output
    # --------------------------------------------------------------------- #
    alpha <- coef_max[1]
    theta <- coef_max[2:(q1 + 1)]
    phi <- coef_max[q1 + 2]
    beta <- coef_max[(q1 + 3):length(coef_max)]

    z$alpha <- alpha
    z$theta <- theta
    z$phi <- phi
    z$beta <- beta

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ma
    # ------------------------------------------------------------------------
    output_inf_matrix_max <- inf_matrix_max_ridge(
      y = y,
      ma = ma,
      X = X,
      alpha = z$alpha,
      theta = z$theta,
      phi = z$phi,
      beta = z$beta,
      link = link,
      penalty = penalty
    )

    fisher_info_mat <- output_inf_matrix_max$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_max$muhat
    z$fitted <- output_inf_matrix_max$fitted
    z$etahat <- output_inf_matrix_max$etahat
    z$errorhat <- output_inf_matrix_max$errorhat

    # ----------------------------------------------------------------------- #
    # Forecasting
    # ----------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted
    errorhat <- z$errorhat

    X_prev <- rbind(X, X_hat)
    for (i in 1:h1) {
      ynew_prev[n + i] <- alpha +
        X_prev[n + i, ] %*% as.matrix(beta) + theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      errorhat[n + i] <- 0
    }

  }

  # ============================================================================
  # Final Values
  # ============================================================================

  # Inverse of Fisher information matrix
  vcov <- try(solve(fisher_info_mat, tol = 1e-20), silent = TRUE)

  # Check if the inverse matrix is possible ----------------------------- #
  if (inherits(vcov, "try-error")) {
    warning("FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
    # Set status_inv_inf_matrix to 1
    # status_inv_inf_matrix indicates if the Fisher's information matrix is
    # invertible:
    #             - 0 indicates a invertible matrix
    #             - 1 indicates a non-invertible matrix
    z$status_inv_inf_matrix <- 1
  } else {
    # output
    z$start_values <- start_values
    z$fisher_info_mat <- fisher_info_mat
    z$forecast <- y_prev[(n + 1):(n + h1)]

    z$status_inv_inf_matrix <- 0

    # output inverse of Fisher information matrix
    z$vcov <- vcov

    # ---------------------------------------------------------------------
    # Model presentation
    # ---------------------------------------------------------------------
    # Standard error
    stderror <- sqrt(diag(vcov))
    z_stderror <- stderror

    z_zstat <- abs(z$coef / stderror)
    z_pvalues <- 2 * (1 - pnorm(z_zstat))

    model_presentation <- cbind(
      round(z$coef, 4),
      round(z_stderror, 4),
      round(z_zstat, 4),
      round(z_pvalues, 4)
    )

    colnames(model_presentation) <- c(
      "Estimate",
      "Std. Error",
      "z value",
      "Pr(>|z|)"
    )

    z$stderror <- stderror
    z$model <- model_presentation

    z$link <- link

    # ---------------------------------------------------------------------
    # information criteria
    # ---------------------------------------------------------------------
    if (has_X) {
      aux_info1 <- -2 * z$loglik
      # aux_info2 <- p1 + q1 + 2 + length(beta)
      aux_info2 <- length(z$coef)
      log_n <- log(n)

      z$aic <- aux_info1 + 2 * aux_info2
      z$bic <- aux_info1 + log_n * aux_info2
      z$hq <- aux_info1 + log(log_n) * 2 * aux_info2
    } else {
      aux_info1 <- -2 * z$loglik
      # aux_info2 <- p1 + q1 + 2
      aux_info2 <- length(z$coef)
      log_n <- log(n)

      z$aic <- aux_info1 + 2 * aux_info2
      z$bic <- aux_info1 + log_n * aux_info2
      z$hq <- aux_info1 + log(log_n) * 2 * aux_info2
    }

    # ---------------------------------------------------------------------
    # Residuals in the predictor scale
    # ---------------------------------------------------------------------

    # Define objects necessary for computing residuals
    fitted_res <- z$fitted[(m + 1):n]
    etahat_res <- z$etahat[(m + 1):n]
    ynew_res <- ynew[(m + 1):n]
    phi_res <- z$phi

    # Computing residuals

    # V(mu_t) = mu_t(1 - mu_t)
    resid2_Vmu <- fitted_res * (1 - fitted_res)

    # V(mu_t) / (1 + phi) = mu_t(1 - mu_t) / (1 + phi)
    resid2_Vmu_phi <- resid2_Vmu / (1 + phi_res)

    # Note that 1 / g'(mu) = dmu_deta
    resid2_dmu_deta2 <- 1 / mu.eta(eta = etahat_res)^2

    # Var(g(y_t) − eta_t) = (g'(mu_t))^2 * V(mu_t) / (1 + phi)
    resid2_var_res_eta <- resid2_dmu_deta2 * resid2_Vmu_phi

    # g(y_t) − eta_t
    resid2_res_eta <- ynew_res - etahat_res

    # Final
    resid2_final <- resid2_res_eta / sqrt(resid2_var_res_eta)

    # output
    z$resid2 <- resid2_final
  }
  # check inv matrix - end -------------------------------------------- #

  return(z)
}
