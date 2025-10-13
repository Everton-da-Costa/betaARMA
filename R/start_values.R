#' Generate Initial Values for βARMA Model Estimation
#'
#' @description
#' This function calculates reasonable starting values for the parameters of
#' various Beta Autoregressive Moving Average (βARMA) models. The method is
#' based on the approach proposed by Ferrari & Cribari-Neto (2004) for
#' beta regression, adapted here for the time series context.
#'
#' @details
#' The function computes initial values by fitting a linear model
#' (`lm.fit`) to the link-transformed response variable `g(y)`.
#' This provides a computationally cheap and stable way to initialize the
#' main optimization algorithm.
#'
#' The specific procedure depends on the model's structure:
#' \itemize{
#'   \item For models with autoregressive (AR) terms, the lagged values of
#'     the transformed series `g(y)` are used as predictors to estimate
#'     the initial `alpha` (intercept) and `varphi` (AR) coefficients.
#'   \item If exogenous variables `X` are included, they are added as
#'     predictors in the linear model to obtain initial `beta` values.
#'   \item Moving average (MA) coefficients (`theta`) are initialized to zero,
#'     a standard practice in ARMA model estimation.
#'   \item The initial value for the precision parameter `phi` is derived
#'     from the variance of the residuals of the initial linear fit,
#'     following the methodology from Ferrari & Cribari-Neto (2004).
#'   \item For a pure βMA model (no AR or X components), a simpler method is
#'     used where `alpha` is the mean of `g(y)` and `phi` is based on the
#'     unconditional variance of `y`.
#' }
#'
#' This function is called internally by the main model fitting function
#' but is exported for standalone use and inspection.
#'
#' @author
#' Original code by Fabio M. Bayer (bayer@ufsm.br).
#' Substantially modified and improved by Everton da Costa
#' (everto_cost@gmail.com).
#'
#' @references
#' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for
#' modelling rates and proportions. *Journal of Applied Statistics*,
#' 31(7), 799-818. <doi:10.1080/0266476042000214501>
#'
#' @param y A numeric time series with values in the open interval (0, 1).
#' @param link A string specifying the link function for the mean, such as
#'   `"logit"`, `"probit"`, or `"cloglog"`.
#' @param ar A numeric vector of autoregressive (AR) lags. Defaults to `NA`
#'   for models without an AR component.
#' @param ma A numeric vector of moving average (MA) lags. Defaults to `NA`
#'   for models without an MA component.
#' @param X An optional numeric matrix or data frame of exogenous
#'   variables (regressors).
#'
#' @importFrom stats lm.fit fitted residuals var
#'
#' @return
#' A named numeric vector containing the initial values for the model
#' parameters (`alpha`, `varphi`, `theta`, `phi`, `beta`), ready to be
#' used by an optimization routine. Returns `NULL` if the model
#' specification is not recognized.
#'
#' @keywords external
start_values <- function(y, link,
                         ar = NA, ma = NA, X = NA) {
  # from:
  #       Beta Regression for Modelling Rates and Proportions
  #       Silvia Ferrari & Francisco Cribari-Neto
  #       p. 805

  # Link functions
  # ----------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  # d(mu)/d(eta)
  mu.eta  <- link_structure$mu.eta

  ynew <- linkfun(y)
  n <- length(y)

  # Determine model components presence
  has_ar <- !any(is.na(ar))
  has_ma <- !any(is.na(ma))
  has_X <- !is.null(X) && !all(is.na(X)) &&
    (is.matrix(X) || is.data.frame(X))

  # Define p, q (max lags) and p1, q1 (number of parameters)
  p <- ifelse(has_ar, max(ar), 0)
  q <- ifelse(has_ma, max(ma), 0)
  p1 <- ifelse(has_ar, length(ar), 0)
  q1 <- ifelse(has_ma, length(ma), 0)
  m <- max(p, q) # Max lag for initial values

  if (has_X) {
    X <- as.matrix(X) # Ensure X is a matrix if it's used
  }

  # Create variable names for the AR, MA and X components conditionally
  names_varphi <- if (has_ar) paste("varphi", ar, sep = "") else character(0)
  names_theta <- if (has_ma) paste("theta", ma, sep = "") else character(0)
  names_beta <- if (has_X) colnames(X) else character(0)

  # ========================================================================= #
  # BARMA initial values (has_ar, has_ma, !has_X)
  # ========================================================================= #
  if (has_ar && has_ma && !has_X) {

    # ----------------------------------------------------------------------- #
    # P: Matrix of lagged ynew values corresponding to AR lags
    P <- matrix(NA, nrow = n - m, ncol = p1)
    for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

    # ----------------------------------------------------------------------- #
    # Prepare data for initial LM fit: y_start ~ intercept + P
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, P)
    y_start <- linkfun(y[(m + 1):n])

    fit_start  <- lm.fit(x = x_start, y = y_start)

    mqo <- fit_start$coef

    alpha_start <- mqo[1]
    varphi_start <- mqo[-1] # All coefficients after intercept are AR

    # --------------------------------------------------- #
    # precision (phi) initial value calculation
    # --------------------------------------------------- #
    k  <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # ------------------------------------------------------------------------
    # theta: MA coefficients initialized to 0
    theta_start <- rep(0, q1)

    # final combined initial values vector
    start_value <- c(alpha_start, varphi_start, theta_start, phi_start)
    names(start_value) <- c("alpha", names_varphi, names_theta, "phi")

    return(start_value)

  }

  # ============================================================================
  # BAR initial values (has_ar, !has_ma, !has_X)
  # ============================================================================
  if (has_ar && !has_ma && !has_X) {

    # ----------------------------------------------------------------------- #
    # P: Matrix of lagged ynew values corresponding to AR lags
    P <- matrix(NA, nrow = n - m, ncol = p1)
    for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

    # ----------------------------------------------------------------------- #
    # Prepare data for initial LM fit: y_start ~ intercept + P
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, P)
    y_start <- linkfun(y[(m + 1):n])

    fit_start  <- lm.fit(x = x_start, y = y_start)

    mqo <- fit_start$coef

    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(p1 + 1)]

    # --------------------------------------------------- #
    # precision (phi) initial value calculation
    # --------------------------------------------------- #
    k  <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # ------------------------------------------------------------------------
    # final combined initial values vector
    start_value <- c(alpha_start, varphi_start, phi_start)
    names(start_value) <- c("alpha", names_varphi, "phi")

    return(start_value)

  }

  # ============================================================================
  # BMA initial values (!has_ar, has_ma, !has_X)
  # ============================================================================
  if (!has_ar && has_ma && !has_X) {

    # These initial values are simpler for pure MA models without regressors
    mean_y <- mean(y)

    # alpha: n^{-1} \sum_{t=1}^n g(y_t)
    alpha_start <- mean(linkfun(y))

    # theta: MA coefficients initialized to 0
    theta_start <- rep(0, q1)

    # initial value for phi: \bar{y}(1-\bar{y})/var(y)
    phi_start <- (mean_y * (1 - mean_y)) / var(y)

    # ------------------------------------------------------------------------
    # final combined initial values vector
    start_value <- c(alpha_start, theta_start, phi_start)
    names(start_value) <- c("alpha", names_theta, "phi")

    return(start_value)

  }

  # ========================================================================= #
  # BARMAX initial values (has_ar, has_ma, has_X)
  # ========================================================================= #
  if (has_ar && has_ma && has_X) {

    # ----------------------------------------------------------------------- #
    # P: Matrix of lagged ynew values corresponding to AR lags
    P <- matrix(NA, nrow = n - m, ncol = p1)
    for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

    # ----------------------------------------------------------------------- #
    # Prepare data for initial LM fit: y_start ~ intercept + P + X
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, P, X[(m + 1):n, , drop = FALSE])
    y_start <- linkfun(y[(m + 1):n])

    fit_start <- lm.fit(x = x_start, y = y_start)

    mqo <- c(fit_start$coef) # Coefficients from LM: (alpha, varphi, beta)

    # --------------------------------------------------- #
    # precision (phi) initial value calculation
    # --------------------------------------------------- #
    k  <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # initial values
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(p1 + 1)]
    theta_start <- rep(0, q1)
    beta_start <- mqo[(p1 + 2):length(mqo)]

    # final combined initial values vector
    start_value <-
      c(alpha_start, varphi_start, theta_start, phi_start, beta_start)

    names(start_value) <-
      c("alpha", names_varphi, names_theta, "phi", names_beta)

    return(start_value)

  }

  # ========================================================================= #
  # BARX initial values (has_ar, !has_ma, has_X)
  # ========================================================================= #
  if (has_ar && !has_ma && has_X) {

    # ----------------------------------------------------------------------- #
    # P: Matrix of lagged ynew values corresponding to AR lags
    P <- matrix(NA, nrow = n - m, ncol = p1)
    for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

    # ----------------------------------------------------------------------- #
    # Prepare data for initial LM fit: y_start ~ intercept + P + X
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, P, X[(m + 1):n, , drop = FALSE])
    y_start <- linkfun(y[(m + 1):n])

    fit_start <- lm.fit(x = x_start, y = y_start)

    mqo <- c(fit_start$coef) # Coefficients from LM: (alpha, varphi, beta)

    # --------------------------------------------------- #
    # precision (phi) initial value calculation
    # --------------------------------------------------- #
    k  <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # initial values
    alpha_start <- mqo[1]
    varphi_start <- mqo[2:(p1 + 1)]
    beta_start <- mqo[(p1 + 2):length(mqo)]

    # final combined initial values vector
    start_value <- c(alpha_start, varphi_start, phi_start, beta_start)
    names(start_value) <- c("alpha", names_varphi, "phi", names_beta)

    return(start_value)

  }

  # ========================================================================= #
  # BMAX initial values (!has_ar, has_ma, has_X)
  # ========================================================================= #
  if (!has_ar && has_ma && has_X) {

    # ----------------------------------------------------------------------- #
    # Prepare data for initial LM fit: y_start ~ intercept + X
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, X[(m + 1):n, , drop = FALSE])
    y_start <- linkfun(y[(m + 1):n])

    fit_start <- lm.fit(x = x_start, y = y_start)

    mqo <- fit_start$coef # Coefficients from LM: (alpha, beta)

    # --------------------------------------------------- #
    # precision (phi) initial value calculation
    # --------------------------------------------------- #
    k  <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # initial values
    alpha_start <- mqo[1]
    theta_start <- rep(0, q1)
    beta_start <- mqo[2:length(mqo)]

    # final combined initial values vector
    start_value <- c(alpha_start, theta_start, phi_start, beta_start)
    names(start_value) <- c("alpha", names_theta, "phi", names_beta)

    return(start_value)

  }

  # If no matching model configuration is found
  warning("No matching model configuration found for initial values.")
  return(NULL)

}
