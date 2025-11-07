#' @title Beta Autoregressive Moving Average (\eqn{\beta}ARMA) Model Estimation
#'
#' @description
#' Fits \eqn{\beta}ARMA, \eqn{\beta}AR, \eqn{\beta}MA, \eqn{\beta}ARMAX,
#' \eqn{\beta}ARX, or \eqn{\beta}MAX models for time series data
#' in (0,1) using conditional maximum likelihood (errors on predictor scale).
#'
#' @details
#' The `barma` function models time series `y` (values in (0,1)) like rates or
#' proportions. It uses the L-BFGS algorithm for optimization via the `lbfgs`
#' package.
#'
#' Model specification depends on `ar`, `ma`, and `X` (exogenous regressors):
#' \itemize{
#'   \item BARMA(p,q): `ar`, `ma` given.
#'   \item BAR(p): `ar` given, `ma` is `NA`.
#'   \item BMA(q): `ma` given, `ar` is `NA`.
#'   \item BARMAX(p,q,k): `ar`, `ma`, `X` (k regressors) given.
#'   \item BARX(p,k): `ar`, `X` given, `ma` is `NA`.
#'   \item BMAX(q,k): `ma`, `X` given, `ar` is `NA`.
#' }
#' \eqn{g(\mu_t)} is the link function applied to the mean \eqn{\mu_t}.
#' \eqn{\epsilon_t = g(y_t) - g(\mu_t)} are predictor-scale errors.
#' If `X` is used, `ar` or `ma` must also be specified.
#' Forecasts assume future predictor-scale errors are zero.
#'
#' \strong{Key Enhancements in this Implementation:}
#' This version incorporates several improvements over earlier iterations,
#' including (but not restricted to):
#' \itemize{
#'   \item Code refactoring following the tidyverse style guide for clarity.
#'   \item Rigorous validation of the implementation against the original
#'     methodology by Rocha & Cribari-Neto (2009, 2017).
#'   \item Adoption of the L-BFGS optimization algorithm (via the `lbfgs`
#'     package) for robust and efficient parameter estimation.
#'   \item Enhanced computational efficiency through the avoidance of
#'     redundant calculations within the estimation routines.
#' }
#'
#' @param y A `ts` object with values in (0,1).
#' @param ar Numeric vector of AR lags (e.g., `c(1,2)`). Default: `NA` (no AR).
#' @param ma Numeric vector of MA lags (e.g., `c(1)`). Default: `NA` (no MA).
#' @param link Link function: "logit" (default), "probit", "cloglog", "loglog".
#'             Connects mean \eqn{\mu_t \in (0,1)} to linear predictor
#'             \eqn{\eta_t}.
#' @param h1 Integer, forecast horizon. Default: 12.
#' @param X Optional. Numeric matrix or `ts` of exogenous regressors.
#'          Rows must match `length(y)`. Default: `NA`.
#' @param X_hat Optional. Future values for `X` if `h1 > 0`. Must have `h1` rows.
#'              Default: `NA`.
#'
#' @return A list containing:
#' \item{coef}{Named vector of estimated parameters
#' (\eqn{\alpha, \varphi, \theta, \phi, \beta}).}
#' \item{vcov}{Variance-covariance matrix of estimates.}
#' \item{loglik}{Maximized conditional log-likelihood.}
#' \item{aic, bic, hq}{Information criteria (AIC, BIC, Hannan-Quinn).}
#' \item{model}{Summary table: estimates, std. errors, z-values, p-values.}
#' \item{fitted}{`ts` object of fitted mean values \eqn{\hat{\mu}_t}.}
#' \item{forecast}{Vector of `h1` point forecasts for \eqn{y_t} if `h1 > 0`.}
#' \item{resid2}{Standardized residuals on the predictor scale.}
#' \item{link}{Link function used.}
#' \item{n_obs}{Number of observations in `y`.}
#' \item{phi, alpha, varphi, theta, beta}{Individual estimated parameters.}
#' \item{etahat, muhat, errorhat}{Estimated predictors, means, and
#' predictor-scale errors.}
#' \item{fisher_info_mat}{Observed Fisher Information Matrix.}
#' \item{start_values}{Initial parameters for optimization.}
#' \item{conv}{Convergence code from `lbfgs` (0 = success).}
#' \item{inv_inf_matrix}{Indicator for Fisher matrix invertibility (0 = success).}
#' \item{opt (`opt_arma`, etc.)}{Raw `lbfgs` optimization output.}
#'
#' @keywords timeseries regression beta BARMA ARMA
#' @export
barma_old <- function(y,
                      ar = NA, ma = NA,
                      link = "logit",
                      h1 = 12,
                      X = NA, X_hat = NA) {

  # ---------------------------------------------------------------------------
  # Link functions
  # ---------------------------------------------------------------------------

  # Check if the link argument is a character string or an expression
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      linktemp <- eval(link)
    }
  }

  # Set up the link function based on the provided link argument
  if (linktemp == "logit") {

    # Use make.link function for logit link (efficient C implementation)
    stats <- make.link("logit")

  } else if (linktemp == "probit") {

    # Use make.link function for probit link
    stats <- make.link("probit")

  } else if (linktemp == "cloglog") {

    # Use make.link function for cloglog link
    stats <- make.link("cloglog")

  } else if (linktemp == "loglog") {

    # Manually define linkfun, linkinv, and mu.eta for loglog link
    stats <- list()

    stats$linkfun <- function(mu) -log(-log(mu))
    stats$linkinv <- function(eta) exp(-exp(-eta))
    stats$mu.eta <- function(eta) exp(-exp(-eta)) * exp(-eta)

  } else {
    # If the provided link is not supported, raise an error
    stop(
      paste0(
        "The '", linktemp, "' link function is not available.", "\n",
        "Available options are: ", "logit, loglog, cloglog, and probit."
      )
    )
  }

  # Create a list with the link function details
  link1 <- structure(list(
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta = stats$mu.eta
  ))

  # Assign link function components to separate variables
  linkfun <- link1$linkfun
  linkinv <- link1$linkinv
  mu.eta <- link1$mu.eta

  # ---------------------------------------------------------------------------

  if (min(y) <= 0 || max(y) >= 1) stop("OUT OF RANGE (0,1)!")

  if (is.ts(y) == TRUE) {
    freq <- frequency(y)
  } else {
    stop("Data must be a time-series object")
  }

  maxit1 <- 1000

  if (any(is.na(ar)) == FALSE) names_varphi <- paste("varphi", ar, sep = "")
  if (any(is.na(ma)) == FALSE) names_theta  <- paste("theta", ma, sep = "")

  p <- max(ar)
  q <- max(ma)
  n <- length(y)

  p1 <- length(ar)
  q1 <- length(ma)

  ynew <- linkfun(y)
  ystar <- log(y / (1 - y))

  names_beta <- colnames(X)
  y_prev <- c(rep(NA, (n + h1)))

  # ---------------------------------------------------------------------------
  # Initializing the parameter values
  # ---------------------------------------------------------------------------
  startValues <- function(y, link, ar = NA, ma = NA, X = NA) {

    # ---------------------------------------------------------------------- #
    # article:
    #         Beta Regression for Modelling Rates and Proportions
    #         Silvia Ferrari & Francisco Cribari-Neto
    #         p. 805
    # ---------------------------------------------------------------------- #

    # ---------------------------------------------------------------------- #
    # BARMA initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == FALSE) &&
        any(is.na(ma) == FALSE) &&
        any(is.na(X) == TRUE)) {

      m <- max(p, q)

      P <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

      # -------------------------------------------------------------------- #
      n <- length(y)
      m <- max(p, q, na.rm = TRUE)

      x_inter <- matrix(1, nrow = n - m, ncol = 1)
      x_start <- cbind(x_inter, P)
      y_start <- linkfun(y[(m + 1):n])

      fit_start  <- lm.fit(x = x_start, y = y_start)

      mqo <- fit_start$coef

      alpha_start <- mqo[1]
      varphi_start <- mqo[-1]

      # --------------------------------------------------- #
      # precision
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
      # initial value of \theta estimate
      theta_start <- rep(0, q1)

      # final
      start_value <- c(alpha_start,
                       varphi_start,
                       theta_start,
                       phi_start)

      # output
      return(start_value)

    }


    # ---------------------------------------------------------------------- #
    # BAR initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == FALSE) &&
        any(is.na(ma) == TRUE) &&
        any(is.na(X) == TRUE)) {

      m <- max(p, na.rm = TRUE)

      # -------------------------------------------------------------------- #
      P <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

      x_inter <- matrix(1, nrow = n - m, ncol = 1)
      x_start <- cbind(x_inter, P)
      y_start <- linkfun(y[(m + 1):n])

      fit_start  <- lm.fit(x = x_start, y = y_start)

      mqo <- fit_start$coef

      alpha_start <- mqo[1]
      varphi_start <- mqo[-1]

      # --------------------------------------------------- #
      # precision
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
      # final
      start_value <- c(alpha_start,
                       varphi_start,
                       phi_start)

      # output
      return(start_value)

    }


    # ---------------------------------------------------------------------- #
    # BMA initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == TRUE) &&
        any(is.na(ma) == FALSE) &&
        any(is.na(X) == TRUE)) {
      
      ynew <- linkfun(y)

      mean_y <- mean(y)
      one_mean_y <- (1 - mean_y)

      alpha_start <- mean(ynew)
      theta_start <- rep(0, q1)

      phi_start <- (mean_y * one_mean_y) / var(y)

      # ------------------------------------------------------------------------
      start_value <- c(alpha_start, theta_start, phi_start)


      return(start_value)

    }


    # ---------------------------------------------------------------------- #
    # BARMAX initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == FALSE) &&
        any(is.na(ma) == FALSE) &&
        any(is.na(X) == FALSE)) {

      # -------------------------------------------------------------------- #
      m <- max(p, q, na.rm = TRUE)
      X <- as.matrix(X)

      P <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

      # -------------------------------------------------------------------- #
      x_inter <- matrix(1, nrow = n - m, ncol = 1)
      x_start <- cbind(x_inter, P, X[(m + 1):n, ])
      y_start <- linkfun(y[(m + 1):n])

      fit_start <- lm.fit(x = x_start, y = y_start)

      mqo <- c(fit_start$coef)

      # --------------------------------------------------- #
      # precision
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

      start_value <- c(alpha_start,
                       varphi_start,
                       theta_start,
                       phi_start,
                       beta_start)

      # output: initial values
      return(start_value)

    }

    # ---------------------------------------------------------------------- #
    # BARX initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == FALSE) &&
        any(is.na(ma) == TRUE) &&
        any(is.na(X) == FALSE)) {

      # -------------------------------------------------------------------- #
      m <- max(p, q, na.rm = TRUE)
      X <- as.matrix(X)

      P <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

      # -------------------------------------------------------------------- #
      x_inter <- matrix(1, nrow = n - m, ncol = 1)
      x_start <- cbind(x_inter, P, X[(m + 1):n, ])
      y_start <- linkfun(y[(m + 1):n])

      fit_start <- lm.fit(x = x_start, y = y_start)

      mqo <- c(fit_start$coef)

      # --------------------------------------------------- #
      # precision
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

      start_value <- c(alpha_start,
                       varphi_start,
                       phi_start,
                       beta_start)

      # output: initial values
      return(start_value)

    }

    # ---------------------------------------------------------------------- #
    # BMAX initial values
    # ---------------------------------------------------------------------- #
    if (any(is.na(ar) == TRUE) &&
        any(is.na(ma) == FALSE) &&
        any(is.na(X) == FALSE)) {

      m <- max(q, na.rm = TRUE)
      X <- as.matrix(X)

      # -------------------------------------------------------------------- #
      x_inter <- matrix(1, nrow = n - m, ncol = 1)
      x_start <- cbind(x_inter, X[(m + 1):n, ])
      y_start <- linkfun(y[(m + 1):n])

      fit_start <- lm.fit(x = x_start, y = y_start)

      mqo <- c(fit_start$coef)

      # beta initial value
      beta_start <- mqo[(p1 + 2):length(mqo)]

      # alpha, theta, and phi start
      mean_y <- mean(y)
      one_mean_y <- (1 - mean_y)

      alpha_start <- log(mean_y / one_mean_y)
      theta_start <- rep(0, q1)

      phi_start <- (mean_y * one_mean_y) / var(y)

      start_value <- c(alpha_start,
                       theta_start,
                       phi_start,
                       beta_start)

      # output: initial values
      return(start_value)

    }

  }


  # ====================================================================== #
  # BARMA model
  # ====================================================================== #
  if (any(is.na(ar) == FALSE) &&
      any(is.na(ma) == FALSE) &&
      any(is.na(X) == TRUE)) {

    print("\beta ARMA model")

    z <- list()

    m <- max(p, q, na.rm = T)
    z$n_obs <- n

    P <- matrix(nrow = n - m, ncol = p1)
    for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

    # initial values
    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    # -------------------------------------------------------------------- #
    # Log-likelihood ARMA
    # -------------------------------------------------------------------- #
    loglik_arma <- function(z) {

      # --------------------------------------------------------------------- #
      alpha  <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta  <- z[(p1 + 2):(p1 + q1 + 1)]
      phi    <- z[p1 + q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta   <- rep(NA, n)
      mu    <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t]    <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
        mu[t]     <- linkinv(eta[t])
        error[t]  <- ynew[t] - eta[t]
      }

      mu1  <- mu[(m + 1):n]
      y1   <- y[(m + 1):n]

      ll_terms_arma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      final <- -1 * sum(ll_terms_arma)

      return(final)

    }


    # -------------------------------------------------------------------- #
    # Score Vector ARMA
    # -------------------------------------------------------------------- #
    score_arma <- function(z) {

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta <- z[(p1 + 2):(p1 + q1 + 1)]
      phi  <- z[p1 + q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta   <- rep(NA, n)
      mu    <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t]    <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
        mu[t]     <- linkinv(eta[t])
        error[t]  <- ynew[t] - eta[t]
      }

      eta1 <- eta[(m + 1):n]
      mu1  <- mu[(m + 1):n]
      y1   <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      mT     <- diag(mu.eta(eta = eta1))
      ystar  <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

      R <- matrix(nrow = n - m, ncol = q1)
      for (t in 1:(n - m))  R[t, ] <- error[t + m - ma]

      # --------------------------------------------------------------------- #
      # Recursive
      # --------------------------------------------------------------------- #
      deta_dalpha  <- rep(0, n)
      deta_dvarphi <- matrix(0, nrow = n, ncol = p1)
      deta_dtheta  <- matrix(0, nrow = n, ncol = q1)

      for (t in (m + 1):n) {
        deta_dalpha[t]    <- 1          - theta %*% deta_dalpha[t - ma]
        deta_dvarphi[t, ] <- P[t - m, ] - theta %*% deta_dvarphi[t - ma, ]
        deta_dtheta[t, ]  <- R[t - m, ] - theta %*% deta_dtheta[t - ma, ]
      }

      s  <- deta_dalpha[(m + 1):n]
      rP <- deta_dvarphi[(m + 1):n, ]
      rR <- deta_dtheta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      ystar_mustar <- ystar - mustar
      mT_ystar_mustar <- crossprod(mT, ystar_mustar)

      U_alpha  <- phi * crossprod(s, mT_ystar_mustar)
      U_varphi <- phi * crossprod(rP, mT_ystar_mustar)
      U_theta  <- phi * crossprod(rR, mT_ystar_mustar)

      U_phi   <- sum(mu1 * ystar_mustar + log(1 - y1)
                     - digamma((1 - mu1) * phi) + digamma(phi))

      escore_vec <- c(U_alpha, U_varphi, U_theta, U_phi)

      # --------------------------------------------------------------------- #

      final <- -1 * escore_vec

      return(final)

    }


    # -------------------------------------------------------------------- #
    # Fisher information matrix
    # -------------------------------------------------------------------- #
    inf_matrix_arma <- function(y, ar, ma,
                                alpha = 0, varphi = 0, theta = 0,
                                phi = 0, link = link) {

      # -------------------------------------------------------------------- #
      P <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]

      # -------------------------------------------------------------------- #
      errorhat <- rep(0, n)
      muhat <- rep(0, n)
      etahat <- rep(NA, n)

      for (t in (m + 1):n) {
        etahat[t] <- alpha + varphi %*% ynew[t - ar] + theta%*%errorhat[t - ma]
        muhat[t] <- linkinv(etahat[t])
        errorhat[t] <- ynew[t] - etahat[t]
      }

      etahat1 <- etahat[(m + 1):n]
      muhat1 <- muhat[(m + 1):n]
      y1 <- y[(m + 1):n]

      # -------------------------------------------------------------------- #
      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m)) R[i, ] <- errorhat[i + m - ma]

      # Recursive ---------------------------------------------------------- #
      deta_dalpha <- rep(0, n)
      deta.dvarphi <- matrix(0, ncol = p1, nrow = n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
        deta.dvarphi[i, ] <- P[(i - m), ] - theta %*% deta.dvarphi[i - ma, ]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
      }

      s <- deta_dalpha[(m + 1):n]
      rP <- deta.dvarphi[(m + 1):n, ]
      rR <- deta_dtheta[(m + 1):n, ]

      # Recursive - end -------------------------------------------------- #

      # Precompute some values
      one_muhat <- 1 - muhat1

      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma(one_muhat * phi)
      vc <- phi * (psi1 * muhat1 - psi2 * one_muhat)
      D <- diag(psi1 * muhat1^2 + psi2 * one_muhat^2 - trigamma(phi))

      # Precompute vector mu_eta and matrix mT
      mu_eta <- mu.eta(eta = etahat1)
      mT <- diag(mu_eta)

      # Precompute W %*% s, W %*% rP, W %*% rR and mT %*% rR
      W <- diag(phi * (psi1 + psi2)) %*% mT^2
      W_diag <- diag(W)

      W_s <- W_diag *  s
      W_rP <- W_diag * rP
      W_rR <- W_diag *  rR
      mT_vc <- mu_eta * vc

      # Compute the intermediate matrices using crossprod()
      K_a_a <- phi * crossprod(s, W_s)
      K_p_a <- phi * crossprod(rP, W_s)
      K_t_a <- phi * crossprod(rR, W_s)

      K_p_p <- phi * crossprod(rP, W_rP)
      K_p_t <- phi * crossprod(rP, W_rR)
      K_t_t <- phi * crossprod(rR, W_rR)

      K_a_phi <- crossprod(s, mT_vc)
      K_p_phi <- crossprod(rP, mT_vc)
      K_t_phi <- crossprod(rR, mT_vc)
      K_phi_phi <- sum(diag(D))

      # Compute the remaining elements
      K_a_p <- t(K_p_a)
      K_t_p <- t(K_p_t)
      K_a_t <- t(K_t_a)
      K_phi_a <- K_a_phi
      K_phi_p <- t(K_p_phi)
      K_phi_t <- t(K_t_phi)

      # Construct the matrix
      # -------------------------------------------------------------------- #
      fisher_info_mat <- rbind(
        cbind(K_a_a, K_a_p, K_a_t, K_a_phi),
        cbind(K_p_a, K_p_p, K_p_t, K_p_phi),
        cbind(K_t_a, K_t_p, K_t_t, K_t_phi),
        cbind(K_phi_a, K_phi_p, K_phi_t, K_phi_phi)
      )

      # -------------------------------------------------------------------- #
      if (any(is.na(ar)) == F) names_varphi <- paste("varphi", ar, sep = "")
      if (any(is.na(ma)) == F) names_theta  <- paste("theta", ma, sep = "")

      names_fisher_info_mat <- c("alpha", names_varphi, names_theta, "phi")
      colnames(fisher_info_mat) <- names_fisher_info_mat
      rownames(fisher_info_mat) <- names_fisher_info_mat

      # -------------------------------------------------------------------- #

      output_list <- list()

      output_list$fisher_info_mat <- fisher_info_mat

      output_list$muhat <- muhat
      output_list$etahat <- etahat
      output_list$errorhat <- errorhat

      # fitted values
      # -------------------------------------------------------------------- #
      fitted_values <- ts(c(rep(NA, m), muhat[(m + 1):n]),
                          start = start(y),
                          frequency = frequency(y))

      # -------------------------------------------------------------------- #
      output_list$fitted <- fitted_values

      return(output_list)

    }


    # -------------------------------------------------------------------- #
    names_par_arma <- c("alpha", names_varphi, names_theta, "phi")

    # -------------------------------------------------------------------- #
    # optimization
    # -------------------------------------------------------------------- #

    # opt_arma <- lbfgs::lbfgs(
    #   call_eval = loglik_arma,
    #   call_grad = score_arma,
    #   vars = start_values,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt_arma <- stats::optim(
      par = start_values,
      fn = loglik_arma,
      gr = score_arma,
      method = "BFGS"
    )

    if (opt_arma$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # convergece status
    z$conv <- opt_arma$convergence
    z$counts <- as.numeric(opt_arma$counts[1])

    # estimates
    coef_arma <- opt_arma$par[1:(p1 + q1 + 2)]
    names(coef_arma) <- c("alpha", names_varphi, names_theta, "phi")

    # output estimates
    z$coef <- coef_arma

    # loglikelihood value
    z$loglik <- -1 * opt_arma$value

    # ------------------------------------------------------------------------
    # Information Fisher Matrix
    # ------------------------------------------------------------------------
    alpha  <- coef_arma[1]
    varphi <- coef_arma[2:(p1 + 1)]
    theta  <- coef_arma[(p1 + 2):(p1 + q1 + 1)]
    phi    <- coef_arma[p1 + q1 + 2]

    # --------------------------------------------------------------------- #
    z$alpha  <- alpha
    z$varphi <- varphi
    z$theta  <- theta
    z$phi    <- phi

    # ------------------------------------------------------------------------
    # Fisher Information Matrix, ARMA
    # ------------------------------------------------------------------------
    output_inf_matrix_arma <- inf_matrix_arma(y = y,
                                              ar = ar,
                                              ma = ma,

                                              alpha = z$alpha,
                                              varphi = z$varphi,
                                              theta = z$theta,
                                              phi = z$phi,

                                              link = link)

    fisher_info_mat <- output_inf_matrix_arma$fisher_info_mat

    # output
    z$fisher_info_mat <- fisher_info_mat
    z$muhat <- output_inf_matrix_arma$muhat
    z$fitted <- output_inf_matrix_arma$fitted
    z$etahat <- output_inf_matrix_arma$etahat

    # ---------------------------------------------------------------------
    # Forecasting
    # ---------------------------------------------------------------------
    errorhat    <- output_inf_matrix_arma$errorhat

    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted

    for (i in 1:h1) {

      ynew_prev[n + i] <- alpha +
        (varphi %*% ynew_prev[n + i - ar]) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      errorhat[n + i] <- 0
    }

  } # End BARMA model


  # =========================================================================
  # BAR model
  # =========================================================================
  if (any(is.na(ar) == FALSE) &&
      any(is.na(ma) == TRUE) &&
      any(is.na(X) == TRUE)) {

    print("\beta AR model")

    z <- list()

    m <- max(p, na.rm = TRUE)
    z$n_obs <- n
    q1 <- 0

    P <- matrix(nrow = n - m, ncol = p1)
    for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

    # initial values
    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    # -------------------------------------------------------------------- #
    # log-likelihood
    # -------------------------------------------------------------------- #
    loglik_ar <- function(z) {

      m <- max(p, na.rm = TRUE)

      # --------------------------------------------------------------------- #
      alpha  <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi   <- z[p1 + 2]

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar]
        mu[t] <- linkinv(eta[t])
      }

      mu1  <- mu[(m + 1):n]
      y1   <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      # log-likelihood - function
      # --------------------------------------------------------------------- #

      ll_terms_ar <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      sum_ll <- -1 * sum(ll_terms_ar)

      return(sum_ll)
    }

    # -------------------------------------------------------------------- #
    # score vector
    # -------------------------------------------------------------------- #
    score_ar <- function(z) {

      m <- max(p, na.rm = TRUE)

      # --------------------------------------------------------------------- #
      alpha  <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi    <- z[p1 + 2]

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar]
        mu[t] <- linkinv(eta[t])
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

      # Precompute vector mu_eta and matrix mT
      mu_eta <- mu.eta(eta = eta1)
      mT <- diag(mu_eta)

      # --------------------------------------------------------------------- #
      ystar_mustar <- ystar - mustar
      mT_ystar_mustar <- mu_eta * ystar_mustar

      U_alpha  <- phi * sum(mu_eta * ystar_mustar)
      U_varphi <- phi * crossprod(P, mT_ystar_mustar)

      U_phi    <- sum(mu1 * ystar_mustar +
                        log(1 - y1) -
                        digamma((1 - mu1) * phi) +
                        digamma(phi))

      score_vec <- c(U_alpha,
                     U_varphi,
                     U_phi)

      # --------------------------------------------------------------------- #

      final <- -1 * score_vec

      return(final)

    }

    names_par_ar <- c("alpha", names_varphi, "phi")

    # opt_ar <- lbfgs::lbfgs(
    #   call_eval = loglik_ar,
    #   call_grad = score_ar,
    #   vars = start_values,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt_ar <- stats::optim(
      par = start_values,
      fn = loglik_ar,
      gr = score_ar,
      method = "BFGS"
    )

    if (opt_ar$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$opt_ar <- opt_ar
    z$conv <- opt_ar$convergence

    # check convergence ----------------------------------------------------- #
    if (z$conv != 0) {

      warning("AR - FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      z$inv_inf_matrix <- 2
      inv_inf_matrix_rest <- z$inv_inf_matrix
      # code 2 implies that don't have matrix to the tests

    } else if (z$conv == 0) {

      coef_ar <- opt_ar$par[1:(p1 + 2)]
      names(coef_ar) <- names_par_ar

      z$coef <- coef_ar
      z$loglik <- -1 * opt_ar$value

      alpha  <- coef_ar[1]
      varphi <- coef_ar[2:(p1 + 1)]
      phi    <- coef_ar[p1 + 2]

      ar_alpha  <- alpha
      ar_varphi <- varphi
      ar_phi    <- phi

      # --------------------------------------------------------------------- #
      z$alpha  <- alpha
      z$varphi <- varphi
      z$phi    <- phi

      # just check
      score_ar_coef <- score_ar(coef_ar)
      z$score_ar_coef <- score_ar_coef

      # ------------------------------------------------------------------------
      # Fisher information matrix - AR
      # ------------------------------------------------------------------------
      etahat   <- rep(NA, n)
      muhat    <- rep(NA, n)

      # E(error)=0
      errorhat <- rep(0, n)

      for (t in (m + 1):n) {
        etahat[t] <- alpha + varphi %*% ynew[t - ar]
        muhat[t] <- linkinv(etahat[t])
        errorhat[t]   <- ynew[t] - etahat[t]
      }
      etahat1 <- etahat[(m + 1):n]
      muhat1  <- muhat[(m + 1):n]
      y1      <- y[(m + 1):n]

      z$etahat <- etahat
      z$muhat <- muhat
      z$fitted <- ts(c(rep(NA, m), muhat[(m + 1):n]),
                     start = start(y),
                     frequency = frequency(y))

      # --------------------------------------------------------------------- #
      # Precompute some values
      one_muhat <- 1 - muhat1

      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma(one_muhat * phi)
      vc <- phi * (psi1 * muhat1 - psi2 * one_muhat)
      D <- diag(psi1 * (muhat1^2) + psi2 * one_muhat^2 - trigamma(phi))

      vI <- as.vector(rep(1, n - m))

      # --------------------------------------------------------------------- #
      # Precompute vector mu_eta and matrix mT
      mu_eta <- mu.eta(eta = etahat1)
      mT <- diag(mu_eta)

      # Precompute W %*% s, W %*% rP
      W <- diag(phi * (psi1 + psi2)) %*% mT^2
      W_diag <- diag(W)

      W_vI <- W_diag * vI
      W_P <- W_diag * P
      mT_vc <- mu_eta * vc

      # Compute the intermediate matrices using crossprod()
      K_a_a <- as.matrix(phi * sum(W_diag))
      K_p_a <- phi * crossprod(P, W_vI)
      K_p_p <- phi * crossprod(P, W_P)

      K_a_phi <- crossprod(vI, mT_vc)
      K_p_phi <- crossprod(P, mT_vc)
      K_phi_phi <- sum(diag(D))

      K_phi_a <- K_a_phi
      K_a_p <- t(K_p_a)
      K_phi_p <- t(K_p_phi)

      fisher_info_mat <- rbind(
        cbind(K_a_a, K_a_p, K_a_phi),
        cbind(K_p_a, K_p_p, K_p_phi),
        cbind(K_phi_a, K_phi_p, K_phi_phi)
      )

      colnames(fisher_info_mat) <- rownames(fisher_info_mat) <- names_par_ar

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

  } # End BAR model


  # =========================================================================
  # BMA model
  # =========================================================================
  if (any(is.na(ar) == TRUE) &&
      any(is.na(ma) == FALSE) &&
      any(is.na(X) == TRUE)) {

    print("\beta MA model")

    z <- list()

    m <- max(q, na.rm = TRUE)
    p1 <- 0

    # initial values
    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    # ---------------------------------------------------------------------- #
    # log-likelihood - MA
    # ---------------------------------------------------------------------- #
    loglik_ma <- function(z) {

      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta   <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t]   <- alpha + theta %*% error[t - ma]
        error[t] <- ynew[t] - eta[t]
      }

      mu1 <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ll_terms_ma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      sum_ll <- -1 * sum(ll_terms_ma)

      return(sum_ll)
    }

    # --------------------------------------------------------------------- #
    # score vector - MA
    # --------------------------------------------------------------------- #
    score_ma <- function(z) {

      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2]

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu  <- rep(NA, n)
      error <- rep(0, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + theta %*% error[t - ma]
        error[t] <- ynew[t] - eta[t]
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- linkinv(eta1)
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m))    R[i, ] <- error[i + m - ma]

      # --------------------------------------------------------------------- #
      # Recursive
      # --------------------------------------------------------------------- #
      deta_dalpha <- rep(0, n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <-              1 - theta %*% deta_dalpha[i - ma]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
      }

      s <- deta_dalpha[(m + 1):n]
      rR <- deta_dtheta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      # Precompute vector mu_eta and matrix mT
      mu_eta <- mu.eta(eta = eta1)
      mT <- diag(mu_eta)

      # Precompute the difference between ystar and mustar
      ystar_mustar <- ystar - mustar

      # Precompute mT %*% (ystar - mustar)
      mT_ystar_mustar <- mu_eta * ystar_mustar

      # Compute U_alpha using crossprod()
      U_alpha <- phi * crossprod(s, mT_ystar_mustar)
      U_theta <- phi * crossprod(rR, mT_ystar_mustar)

      # Compute Uphi
      U_phi <- sum(mu1 * ystar_mustar + log(1 - y1)
                   - digamma((1 - mu1) * phi) + digamma(phi))

      escore_vec <- -1 * c(U_alpha, U_theta, U_phi)

      return(escore_vec)
    }

    names_par_ma <- c("alpha", names_theta, "phi")

    # -------------------------------------------------------------------- #
    # optimization
    # -------------------------------------------------------------------- #
    # opt_ma <- lbfgs::lbfgs(
    #   call_eval = loglik_ma,
    #   call_grad = score_ma,
    #   vars = start_values,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt_ma <- stats::optim(
      par = start_values,
      fn = loglik_ma,
      gr = score_ma,
      method = "BFGS"
    )

    if (opt_ma$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    z$opt_ma <- opt_ma
    z$conv <- opt_ma$convergence


    # check convergence ----------------------------------------------------- #
    if (opt_ma$conv != 0) {

      warning("MA - FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      z$inv_inf_matrix <- 2
      inv_inf_matrix_rest <- z$inv_inf_matrix
      # code 2 implies that don't have matrix to the tests

    } else if (z$conv == 0) {

      coef_ma <- (opt_ma$par)[1:(q1 + 2)]
      names(coef_ma) <- names_par_ma

      z$coef <- coef_ma
      z$loglik <- -1 * opt_ma$value

      alpha <- coef_ma[1]
      theta <- coef_ma[2:(q1 + 1)]
      phi <- coef_ma[q1 + 2]

      z$alpha <- alpha
      z$theta <- theta
      z$phi <- phi

      # -----------------------------------------------------------------------
      # Fisher Information matrix - MA
      # -----------------------------------------------------------------------
      # E(error) = 0

      errorhat <- rep(0, n)
      etahat <- rep(NA, n)
      muhat <- rep(NA, n)

      for (t in (m + 1):n) {
        etahat[t] <- alpha + theta %*% errorhat[t - ma]
        muhat[t] <- linkinv(etahat[t])
        errorhat[t] <- ynew[t] - etahat[t]
      }
      etahat1 <- etahat[(m + 1):n]
      muhat1 <- muhat[(m + 1):n]
      y1 <- y[(m + 1):n]

      z$etahat <- etahat
      z$muhat <- muhat
      z$fitted <- ts(c(rep(NA, m), muhat[(m + 1):n]),
                     start = start(y),
                     frequency = frequency(y))


      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m))  R[i, ] <- errorhat[i + m - ma]

      # --------------------------------------------------------------------- #
      # Recursive
      # --------------------------------------------------------------------- #
      deta_dalpha <- rep(0, n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <-    1         - theta %*% deta_dalpha[i - ma]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
      }

      s <- deta_dalpha[(m + 1):n]
      rR <- deta_dtheta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      # Precompute some values
      one_muhat <- 1 - muhat1

      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma(one_muhat * phi)

      vc <- phi * (psi1 * muhat1 - psi2 * one_muhat)
      D <- diag(c(psi1 * (muhat1^2) + psi2 * one_muhat^2 - trigamma(phi)))

      # Precompute vector mu_eta and matrix mT
      mu_eta <- mu.eta(eta = etahat1)
      mT <- diag(mu_eta)

      # Precompute W %*% s, W %*% rR and mT %*% vc
      W <- diag(phi * (psi1 + psi2)) %*% mT^2
      W_diag <- diag(W)

      W_s <- W_diag * s
      W_rR <- W_diag * rR
      mT_vc <- mu_eta * vc

      # ---------------------------------------- #
      # Compute the intermediate matrices using crossprod()
      K_a_a <- phi * crossprod(s, W_s)
      K_t_a <- phi * crossprod(rR, W_s)
      K_t_t <- phi * crossprod(rR, W_rR)

      K_a_phi <- crossprod(s, mT_vc)
      K_t_phi <-  crossprod(rR, mT_vc)

      # Compute the remaining elements
      K_a_t <- t(K_t_a)
      K_phi_a <- t(K_a_phi)
      K_phi_t <- t(K_t_phi)
      K_phi_phi <- sum(diag(D))

      fisher_info_mat <- rbind(
        cbind(K_a_a, K_a_t, K_a_phi),
        cbind(K_t_a, K_t_t, K_t_phi),
        cbind(K_phi_a, K_phi_t, K_phi_phi)
      )

      # --------------------------------------------------------------------- #
      # Forecasting
      # --------------------------------------------------------------------- #
      ynew_prev <- c(ynew, rep(NA, h1))
      y_prev[1:n] <- z$fitted

      for (i in 1:h1) {

        ynew_prev[n + i] <- alpha + theta %*% errorhat[n + i - ma]
        y_prev[n + i] <- linkinv(ynew_prev[n + i])
        errorhat[n + i] <- 0

      }

    }

  } # End BMA model


  # ====================================================================== #
  # BARMAX model
  # ====================================================================== #
  if (any(is.na(ar) == FALSE) &&
      any(is.na(ma) == FALSE) &&
      any(is.na(X) == FALSE)) {

    print("\beta ARMA X model")

    z <- list()

    m <- max(p, q, na.rm = TRUE)

    X <- as.matrix(X)
    X_hat <- as.matrix(X_hat)
    z$n_obs <- n

    # -------------------------------------------------------------------------
    # Log-likelihood ARMA X
    # -------------------------------------------------------------------------
    loglike_armaX <- function(z) {

      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta <- z[(p1 + 2):(p1 + q1 + 1)]
      phi <- z[p1 + q1 + 2]
      beta <- z[(p1 + q1 + 3):length(z)]

      error <- rep(0, n)
      eta <- rep(NA, n)

      beta_mat <- as.matrix(beta)

      for (i in (m + 1):n) {
        eta[i] <- alpha +
          X[i, ] %*% beta_mat +
          varphi %*% (ynew[i - ar] - X[i - ar, ] %*% beta_mat) +
          theta %*% error[i - ma]

        error[i] <- ynew[i] - eta[i]
      }

      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]

      ll_terms <- dbeta(y1, mu * phi, (1 - mu) * phi, log = TRUE)

      return(sum(-1 * ll_terms))

    }


    # -------------------------------------------------------------------------
    # Score vector ARMA X
    # -------------------------------------------------------------------------
    score_vector_armaX <- function(z) {

      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta <- z[(p1 + 2):(p1 + q1 + 1)]
      phi <- z[p1 + q1 + 2]
      beta <- z[(p1 + q1 + 3):length(z)]

      error <- rep(0, n)
      eta <- rep(NA, n)

      for (i in (m + 1):n) {

        eta[i] <- alpha +
          X[i, ] %*% as.matrix(beta) +

          (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% as.matrix(beta))) +

          (theta %*% error[i - ma])

        error[i] <- ynew[i] - eta[i]

      }

      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

      mT <- diag(mu.eta(eta[(m + 1):n]))

      R <- matrix(rep(NA, (n - m) * q1), ncol = q1)
      for (i in 1:(n - m)) {

        R[i, ] <- error[i + m - ma]

      }

      P <- matrix(rep(NA, (n - m) * p1), ncol = p1)
      for (i in 1:(n - m)) {

        P[i, ] <- ynew[i + m - ar] - X[i + m - ar, ] %*% as.matrix(beta)

      }

      k1 <- length(beta)

      M <- matrix(rep(NA, (n - m) * k1), ncol = k1)
      for (i in 1:(n - m)) {
        for (j in 1:k1) {
          M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
        }
      }


      # Recursive
      deta_dalpha <- rep(0, n)
      deta_dvarphi <- matrix(0, ncol = p1, nrow = n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)
      deta_dbeta <- matrix(0, ncol = k1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
        deta_dvarphi[i, ] <- P[(i - m), ] - theta %*% deta_dvarphi[i - ma, ]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
        deta_dbeta[i, ] <- M[(i - m), ] - theta %*% deta_dbeta[i - ma, ]
      }

      s <- deta_dalpha[(m + 1):n]
      rP <- deta_dvarphi[(m + 1):n, ]
      rR <- deta_dtheta[(m + 1):n, ]
      rM <- deta_dbeta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      # obj. for optimization
      ystar_mustar <- ystar - mustar
      mT_ystar_mustar <- mT %*% ystar_mustar

      U_alpha  <- phi * s     %*% mT_ystar_mustar
      U_varphi <- phi * t(rP) %*% mT_ystar_mustar
      U_theta  <- phi * t(rR) %*% mT_ystar_mustar
      U_beta   <- phi * t(rM) %*% mT_ystar_mustar

      U_phi <- sum(mu * ystar_mustar + log(1 - y1)
                   - digamma((1 - mu) * phi) + digamma(phi))
      # --------------------------------------------------------------------- #

      rval <- c(U_alpha, U_varphi, U_theta, U_phi, U_beta)

      return(-1 * rval)

    }


    # -------------------------------------------------------------------------
    # Fisher information matrix - ARMA X
    # -------------------------------------------------------------------------
    inf_matrix_arma_X <- function(y,
                                  ar, ma, X,
                                  alpha = 0, varphi = 0, theta = 0,
                                  phi = 0, beta = 0) {

      errorhat <- rep(0, n)
      etahat <- rep(NA, n)

      for (i in (m + 1):n) {

        etahat[i] <- alpha +
          X[i, ] %*% as.matrix(beta) +
          (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% as.matrix(beta))) +
          theta %*% errorhat[i - ma]

        errorhat[i] <- ynew[i] - etahat[i]
      }

      muhat1 <- linkinv(etahat[(m + 1):n])
      y1 <- y[(m + 1):n]

      # -------------------------------------------------------------------- #
      R <- matrix(rep(NA, (n - m) * q1), ncol = q1)
      for (i in 1:(n - m)) {

        R[i, ] <- errorhat[i + m - ma]

      }

      P <- matrix(rep(NA, (n - m) * p1), ncol = p1)
      for (i in 1:(n - m)) {

        P[i, ] <- ynew[i + m - ar] - X[i + m - ar, ] %*% as.matrix(beta)

      }

      k1 <- length(beta)

      M <- matrix(rep(NA, (n - m) * k1), ncol = k1)
      for (i in 1:(n - m)) {
        for (j in 1:k1) {
          M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
        }
      }

      # -------------------------------------------------------------------- #
      # recorrencias
      # -------------------------------------------------------------------- #
      deta_dalpha <- rep(0, n)
      deta_dvarphi <- matrix(0, ncol = p1, nrow = n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)
      deta_dbeta <- matrix(0, ncol = k1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
        deta_dvarphi[i, ] <- P[(i - m), ] - theta %*% deta_dvarphi[i - ma, ]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
        deta_dbeta[i, ] <- M[(i - m), ] - theta %*% deta_dbeta[i - ma, ]

      }

      a <- deta_dalpha[(m + 1):n]
      rP <- deta_dvarphi[(m + 1):n, ]
      rR <- deta_dtheta[(m + 1):n, ]
      rM <- deta_dbeta[(m + 1):n, ]

      # -------------------------------------------------------------------- #
      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma((1 - muhat1) * phi)
      mT <- diag(mu.eta(etahat[(m + 1):n]))
      W <- diag(c(phi * (psi1 + psi2))) %*% mT^2
      vc <- phi * (psi1 * muhat1 - psi2 * (1 - muhat1))
      D <- diag(c(psi1 * (muhat1^2) + psi2 * (1 - muhat1)^2 - trigamma(phi)))

      Kaa <- phi * t(a) %*% W %*% a
      Kpa <- phi * t(rP) %*% W %*% a
      Kap <- t(Kpa)
      Kta <- phi * t(rR) %*% W %*% a
      Kat <- t(Kta)
      Kaphi <- t(a) %*% mT %*% vc
      Kphia <- Kaphi
      Kpp <- phi * t(rP) %*% W %*% rP
      Kpt <- phi * t(rP) %*% W %*% rR
      Ktp <- t(Kpt)
      Kpphi <- t(rP) %*% mT %*% vc
      Kphip <- t(Kpphi)
      Ktt <- phi * t(rR) %*% W %*% rR
      Ktphi <- t(rR) %*% mT %*% vc
      Kphit <- t(Ktphi)
      Kphiphi <- sum(diag(D))

      Kba <- phi * t(rM) %*% W %*% a
      Kbb <- phi * t(rM) %*% W %*% rM
      Kbphi <- t(rM) %*% mT %*% vc
      Kbp <- phi * t(rM) %*% W %*% rP
      Kbt <- phi * t(rM) %*% W %*% rR

      Kab <- t(Kba)
      Kphib <- t(Kbphi)
      Kpb <- t(Kbp)
      Ktb <- t(Kbt)

      # -------------------------------------------------------------------- #
      names_varphi <- paste("varphi", ar, sep = "")
      names_theta  <- paste("theta", ma, sep = "")
      names_beta <- colnames(X)

      names_fisher_info_mat <- c("alpha",
                                 names_varphi,
                                 names_theta,
                                 "phi",
                                 names_beta)

      fisher_info_mat <- rbind(
        cbind(Kaa, Kap, Kat, Kaphi, Kab),
        cbind(Kpa, Kpp, Kpt, Kpphi, Kpb),
        cbind(Kta, Ktp, Ktt, Ktphi, Ktb),
        cbind(Kphia, Kphip, Kphit, Kphiphi, Kphib),
        cbind(Kba, Kbp, Kbt, Kbphi, Kbb)
      )

      colnames(fisher_info_mat) <- names_fisher_info_mat
      rownames(fisher_info_mat) <- names_fisher_info_mat

      # -------------------------------------------------------------------- #

      output_list <- list()

      output_list$fisher_info_mat <- fisher_info_mat
      output_list$muhat <- muhat1
      output_list$etahat <- etahat
      output_list$errorhat <- errorhat

      # fitted values
      # -------------------------------------------------------------------- #
      fitted_values <- ts(c(rep(NA, m), muhat1),
                          start = start(y),
                          frequency = frequency(y))

      output_list$fitted <- fitted_values

      return(output_list)

    }

    # -------------------------------------------------------------------------
    # Initial values
    # -------------------------------------------------------------------------
    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    names_par <- c("alpha", names_varphi, names_theta, "phi", names_beta)


    # -------------------------------------------------------------------- #
    # Optimization
    # -------------------------------------------------------------------- #
    # opt <- lbfgs::lbfgs(
    #   vars = start_values,
    #   call_eval = loglike_armaX,
    #   call_grad = score_vector_armaX,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt <- stats::optim(
      par = start_values,
      fn = loglike_armaX,
      gr = score_vector_armaX,
      method = "BFGS"
    )

    if (opt$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    # -------------------------------------------------------------------- #
    # convergece status
    z$conv <- opt$convergence
    # z$counts <- opt$counts
    z$opt <- opt

    # loglikelihood value
    z$loglik <- -1 * opt$value

    # -------------------------------------------------------------------- #
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
    # Fisher Information Matrix, ARMA X
    # ------------------------------------------------------------------------
    output_inf_matrix_arma_X <- inf_matrix_arma_X(
      y = y,
      ar = ar,
      ma = ma,
      X = X,

      alpha = z$alpha,
      varphi = z$varphi,
      theta = z$theta,
      phi = z$phi,
      beta = z$beta
    )

    fisher_info_mat <- output_inf_matrix_arma_X$fisher_info_mat

    # output
    z$muhat <- output_inf_matrix_arma_X$muhat
    z$fitted <- output_inf_matrix_arma_X$fitted
    z$etahat <- output_inf_matrix_arma_X$etahat

    # ------------------------------------------------------------------------
    # Forecasting
    # -----------------------------------------------------------------------
    z$errorhat <- output_inf_matrix_arma_X$errorhat

    errorhat <- z$errorhat
    y_prev[1:n] <- z$fitted
    ynew_prev <- c(ynew, rep(NA, h1))

    X_prev <- rbind(X, X_hat)

    for (i in 1:h1) {

      ynew_prev[n + i] <- alpha +
        X_prev[n + i, ] %*% as.matrix(beta) +
        (varphi %*% (ynew_prev[n + i - ar] -
                       X_prev[n + i - ar, ] %*% as.matrix(beta))) +
        theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])
      errorhat[n + i] <- 0

    }

  } # End BARMAX model


  # ============================================================================
  # BARX model
  # ============================================================================
  if (any(is.na(ar) == FALSE) &&
      any(is.na(ma) == TRUE) &&
      any(is.na(X) == FALSE)) {

    print("\beta AR X Model")

    z <- list()

    q1 <- 0
    m <- max(p, na.rm = TRUE)

    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    # -------------------------------------------------------------------- #
    # Log-likelihood \beta AR X
    # -------------------------------------------------------------------- #
    loglike_arX <- function(z) {

      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]
      beta <- z[(p1 + 3):length(z)]

      # --------------------------------------------------------------------- #
      # E(error)=0
      error <- rep(0, n)
      eta <- rep(NA, n)
      beta_mat <- as.matrix(beta)

      for (i in (m + 1):n) {
        eta[i] <- alpha + X[i, ] %*% beta_mat +
          (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% beta_mat))

        error[i] <- ynew[i] - eta[i]
      }
      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]

      ll_terms <- dbeta(y1, mu * phi, (1 - mu) * phi, log = TRUE)
      return(-1 * sum(ll_terms))

    }

    # -------------------------------------------------------------------- #
    # Score vector \beta AR X
    # -------------------------------------------------------------------- #
    score_vector_arX <- function(z) {

      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]
      beta <- z[(p1 + 3):length(z)]

      # --------------------------------------------------------------------- #
      # E(error)=0
      error <- rep(0, n)
      eta <- rep(NA, n)

      for (i in (m + 1):n) {
        eta[i] <- alpha + X[i, ] %*%
          as.matrix(beta) +
          (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% as.matrix(beta)))
        error[i] <- ynew[i] - eta[i]
      }

      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

      mT <- diag(mu.eta(eta = eta[(m + 1):n]))

      P <- matrix(rep(NA, (n - m) * p1), ncol = p1)
      for (i in 1:(n - m)) {
        P[i, ] <- ynew[i + m - ar] - X[i + m - ar, ] %*% as.matrix(beta)
      }

      M <- matrix(rep(NA, (n - m) * length(beta)), ncol = length(beta))
      for (i in 1:(n - m)) {
        for (j in 1:length(beta)) {
          M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
        }
      }

      # --------------------------------------------------------------------- #
      # obj. for optimization
      ystar_mustar <- ystar - mustar
      mT_ystar_mustar <- mT %*% ystar_mustar

      U_alpha <- phi * sum((mu.eta(eta = eta[(m + 1):n])) * ystar_mustar)
      U_varphi <- phi * t(P) %*% mT_ystar_mustar
      U_beta <- phi * t(M) %*% mT_ystar_mustar

      U_phi <- sum(mu * ystar_mustar +
                     log(1 - y1) -
                     digamma((1 - mu) * phi) +
                     digamma(phi))

      rval <- c(U_alpha, U_varphi, U_phi, U_beta)

      return(-1 * rval)

    }

    # -------------------------------------------------------------------- #
    # Optimization \beta AR X
    # -------------------------------------------------------------------- #
    # opt <- lbfgs::lbfgs(
    #   call_eval = loglike_arX,
    #   call_grad = score_vector_arX,
    #   vars = start_values,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt <- stats::optim(
      par = start_values,
      fn = loglike_arX,
      gr = score_vector_arX,
      method = "BFGS"
    )

    if (opt$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }


    # output
    n_obs <- length(y)
    z$n_obs <- n_obs

    z$conv <- opt$conv

    # -------------------------------------------------------------------- #
    # estimates
    coef <- (opt$par)[1:(p1 + 2 + ncol(X))]
    names_par <- c("alpha", names_varphi, "phi", names_beta)

    names(coef) <- names_par

    # output
    z$coef <- coef

    alpha <- coef[1]
    varphi <- coef[2:(p1 + 1)]
    phi <- coef[p1 + 2]
    beta <- coef[(p1 + 3):length(coef)]

    z$alpha <- alpha
    z$varphi <- varphi
    z$phi <- phi
    z$beta <- beta

    # -------------------------------------------------------------------- #
    # E(error)=0
    errorhat <- rep(0, n)
    etahat <- rep(NA, n)

    for (i in (m + 1):n) {
      etahat[i] <- alpha +
        X[i, ] %*% as.matrix(beta) +
        (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% as.matrix(beta)))

      errorhat[i] <- ynew[i] - etahat[i]

    }
    muhat <- linkinv(etahat[(m + 1):n])

    # -------------------------------------------------------------------- #
    y1 <- y[(m + 1):n]

    z$fitted <- ts(c(rep(NA, m), muhat),
                   start = start(y),
                   frequency = frequency(y))

    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * phi) - digamma((1 - muhat) * phi)


    P <- matrix(rep(NA, (n - m) * p1), ncol = p1)
    for (i in 1:(n - m)) {
      P[i, ] <- ynew[i + m - ar] - X[i + m - ar, ] %*% as.matrix(beta)
    }

    M <- matrix(rep(NA, (n - m) * length(beta)), ncol = length(beta))
    for (i in 1:(n - m)) {
      for (j in 1:length(beta)) {
        M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
      }
    }

    vI <- as.vector(rep(1, n - m))

    psi1 <- trigamma(muhat * phi)
    psi2 <- trigamma((1 - muhat) * phi)
    mT <- diag(mu.eta(etahat[(m + 1):n]))
    W <- diag(c(phi * (psi1 + psi2))) %*% mT^2
    vc <- phi * (psi1 * muhat - psi2 * (1 - muhat))
    D <- diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(phi)))

    Kaa <- as.matrix(phi * sum(diag(W)))
    Kpa <- phi * t(P) %*% W %*% vI
    Kap <- t(Kpa)
    Kaphi <- vI %*% mT %*% vc
    Kphia <- Kaphi
    Kpp <- phi * t(P) %*% W %*% P
    Kpphi <- t(P) %*% mT %*% vc
    Kphip <- t(Kpphi)
    Kphiphi <- sum(diag(D))

    Kba <- phi * t(M) %*% W %*% vI
    Kbb <- phi * t(M) %*% W %*% M
    Kbphi <- t(M) %*% mT %*% vc
    Kbp <- phi * t(M) %*% W %*% P

    Kab <- t(Kba)
    Kphib <- t(Kbphi)
    Kpb <- t(Kbp)

    fisher_info_mat <- rbind(
      cbind(Kaa, Kap, Kaphi, Kab),
      cbind(Kpa, Kpp, Kpphi, Kpb),
      cbind(Kphia, Kphip, Kphiphi, Kphib),
      cbind(Kba, Kbp, Kbphi, Kbb)
    )

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

      errorhat[n + i] <- 0

    }

  } # End BAR model

  # ============================================================================
  # BMAX model
  # ============================================================================
  if (any(is.na(ar) == T) &&
      any(is.na(ma) == F) &&
      any(is.na(X) == F)) {

    print("\beta MA X model")

    z <- list()
    p1 <- 0
    m <- max(q, na.rm = TRUE)

    start_values <- startValues(y, link = link, ar = ar, ma = ma, X = X)

    # -------------------------------------------------------------------- #
    # Log-likelihood \beta MA X
    # -------------------------------------------------------------------- #
    loglike_maX <- function(z) {

      alpha <- z[1]
      theta <- z[(2):(q1 + 1)]
      phi <- z[q1 + 2]
      beta <- z[(q1 + 3):length(z)]

      # E(error)=0
      error <- rep(0, n)
      eta <- rep(NA, n)

      for (i in (m + 1):n) {
        eta[i] <- alpha + X[i, ] %*% as.matrix(beta) + (theta %*% error[i - ma])
        error[i] <- ynew[i] - eta[i]
      }
      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]

      ll <- dbeta(y1, mu * phi, (1 - mu) * phi, log = TRUE)

      return(-1 * sum(ll))

    }

    # -------------------------------------------------------------------- #
    # Score vector \beta MA X
    # -------------------------------------------------------------------- #
    score_vector_maX <- function(z) {

      alpha <- z[1]
      theta <- z[(2):(q1 + 1)]
      phi <- z[q1 + 2]
      beta <- z[(q1 + 3):length(z)]

      # E(error)=0
      error <- rep(0, n)
      eta <- rep(NA, n)

      for (i in (m + 1):n) {
        eta[i] <- alpha + X[i, ] %*% as.matrix(beta) + (theta %*% error[i - ma])
        error[i] <- ynew[i] - eta[i]
      }

      mu <- linkinv(eta[(m + 1):n])
      y1 <- y[(m + 1):n]
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

      mT <- diag(mu.eta(eta[(m + 1):n]))

      R <- matrix(rep(NA, (n - m) * q1), ncol = q1)
      for (i in 1:(n - m)) {
        R[i, ] <- error[i + m - ma]
      }


      k1 <- length(beta)
      M <- matrix(nrow = n - m, ncol = k1)
      for (i in 1:(n - m)) {
        for (j in 1:k1) {
          M[i, j] <- X[i + m, j]
        }
      }

      # Recursive
      deta_dalpha <- rep(0, n)
      deta_dtheta <- matrix(0, ncol = q1, nrow = n)
      deta.dbeta <- matrix(0, ncol = k1, nrow = n)

      for (i in (m + 1):n) {
        deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
        deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
        deta.dbeta[i, ] <- M[(i - m), ] - theta %*% deta.dbeta[i - ma, ]
      }

      s <- deta_dalpha[(m + 1):n]
      rR <- deta_dtheta[(m + 1):n, ]
      rM <- deta.dbeta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      # obj. for optimization
      ystar_mustar <- ystar - mustar
      mT_ystar_mustar <- mT %*% ystar_mustar

      U_alpha <- phi * s %*% mT_ystar_mustar
      U_theta <- phi * t(rR) %*% mT_ystar_mustar
      U_beta <- phi * t(rM) %*% mT_ystar_mustar

      U_phi <- sum(mu * ystar_mustar + log(1 - y1)
                   - digamma((1 - mu) * phi) + digamma(phi))

      rval <- c(U_alpha, U_theta, U_phi, U_beta)

      return(-1 * rval)
    }

    # -------------------------------------------------------------------- #
    # Optimization \beta MA X
    # -------------------------------------------------------------------- #
    # opt <- lbfgs::lbfgs(
    #   call_eval = loglike_maX,
    #   call_grad = score_vector_maX,
    #   vars = start_values,
    #   invisible = 1,
    #   linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
    #   max_iterations = maxit1
    # )
    
    opt <- stats::optim(
      par = start_values,
      fn = loglike_maX,
      gr = score_vector_maX,
      method = "BFGS"
    )

    if (opt$conv != 0) {
      warning("FUNCTION DID NOT CONVERGE!")
    }

    n_obs <- length(y)
    z$n_obs <- n_obs

    names_par <- c("alpha", names_theta, "phi", names_beta)

    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1 + 2 + ncol(X))]
    names(coef) <- names_par
    z$coef <- coef

    alpha <- coef[1]
    theta <- coef[(2):(q1 + 1)]
    phi <- coef[q1 + 2]
    beta <- coef[(q1 + 3):length(coef)]

    z$alpha <- alpha
    z$theta <- theta
    z$phi <- phi

    # -------------------------------------------------------------------- #
    # E(error)=0
    errorhat <- rep(0, n)
    etahat <- rep(NA, n)

    for (i in (m + 1):n) {
      etahat[i] <- alpha + X[i, ] %*%
        as.matrix(beta) +
        theta %*% errorhat[i - ma]
      errorhat[i] <- ynew[i] - etahat[i]
    }
    muhat <- linkinv(etahat[(m + 1):n])

    # -------------------------------------------------------------------- #
    y1 <- y[(m + 1):n]

    z$fitted <- ts(c(rep(NA, m), muhat),
                   start = start(y),
                   frequency = frequency(y))

    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * phi) - digamma((1 - muhat) * phi)

    R <- matrix(rep(NA, (n - m) * q1), ncol = q1)
    for (i in 1:(n - m)) {
      R[i, ] <- errorhat[i + m - ma]
    }

    k1 <- length(beta)
    M <- matrix(rep(NA, (n - m) * k1), ncol = k1)
    for (i in 1:(n - m)) {
      for (j in 1:k1) {
        M[i, j] <- X[i + m, j]
      }
    }

    # Recursive
    deta_dalpha <- rep(0, n)
    deta_dtheta <- matrix(0, ncol = q1, nrow = n)
    deta.dbeta <- matrix(0, ncol = k1, nrow = n)

    for (i in (m + 1):n) {
      deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
      deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
      deta.dbeta[i, ] <- M[(i - m), ] - theta %*% deta.dbeta[i - ma, ]
    }

    a <- deta_dalpha[(m + 1):n]
    rR <- deta_dtheta[(m + 1):n, ]
    rM <- deta.dbeta[(m + 1):n, ]

    psi1 <- trigamma(muhat * phi)
    psi2 <- trigamma((1 - muhat) * phi)
    mT <- diag(mu.eta(etahat[(m + 1):n]))
    W <- diag(c(phi * (psi1 + psi2))) %*% mT^2
    vc <- phi * (psi1 * muhat - psi2 * (1 - muhat))
    D <- diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(phi)))

    Kaa <- phi * t(a) %*% W %*% a
    Kta <- phi * t(rR) %*% W %*% a
    Kat <- t(Kta)
    Kaphi <- t(a) %*% mT %*% vc
    Kphia <- Kaphi
    Ktt <- phi * t(rR) %*% W %*% rR
    Ktphi <- t(rR) %*% mT %*% vc
    Kphit <- t(Ktphi)
    Kphiphi <- sum(diag(D))

    Kba <- phi * t(rM) %*% W %*% a
    Kbb <- phi * t(rM) %*% W %*% rM
    Kbphi <- t(rM) %*% mT %*% vc
    Kbt <- phi * t(rM) %*% W %*% rR

    Kab <- t(Kba)
    Kphib <- t(Kbphi)
    Ktb <- t(Kbt)

    fisher_info_mat <- rbind(
      cbind(Kaa, Kat, Kaphi, Kab),
      cbind(Kta, Ktt, Ktphi, Ktb),
      cbind(Kphia, Kphit, Kphiphi, Kphib),
      cbind(Kba, Kbt, Kbphi, Kbb)
    )

    # ----------------------------------------------------------------------- #
    # Forecasting
    # ----------------------------------------------------------------------- #
    ynew_prev <- c(ynew, rep(NA, h1))
    y_prev[1:n] <- z$fitted

    X_prev <- rbind(X, X_hat)
    for (i in 1:h1) {

      ynew_prev[n + i] <- alpha +
        X_prev[n + i, ] %*% as.matrix(beta) + theta %*% errorhat[n + i - ma]

      y_prev[n + i] <- linkinv(ynew_prev[n + i])

      errorhat[n + i] <- 0
    }

  }  # End BMAX model


  # ============================================================================
  # Final Values
  # ============================================================================

  # Inverse of Fisher information matrix
  vcov <- try(solve(fisher_info_mat, tol = 1e-20), silent = TRUE)


  # Check if the inverse matrix is possible ----------------------------- #
  if (!(typeof(vcov) == "double")) {

    warning("FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
    # Set inv_inf_matrix to 1
    # inv_inf_matrix indicates if the Fisher's information matrix is
    # invertible:
    #             - 0 indicates a invertible matrix
    #             - 1 indicates a non-invertible matrix
    z$inv_inf_matrix <- 1

  } else {

    # output
    z$start_values <- start_values
    z$fisher_info_mat <- fisher_info_mat
    z$forecast <- y_prev[(n + 1):(n + h1)]

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

    model_presentation <- cbind(round(z$coef, 4),
                                round(z_stderror, 4),
                                round(z_zstat, 4),
                                round(z_pvalues, 4))

    colnames(model_presentation) <- c("Estimate",
                                      "Std. Error",
                                      "z value",
                                      "Pr(>|z|)")

    z$model <- model_presentation

    z$link <- link

    # ---------------------------------------------------------------------
    # information criteria
    # ---------------------------------------------------------------------
    if (any(is.na(X) == FALSE)) {

      aux_info1 <- -2 * z$loglik
      aux_info2 <- p1 + q1 + 2 + length(beta)
      log_n <- log(n)

      z$aic <- aux_info1 + 2 * aux_info2
      z$bic <- aux_info1 + log_n * aux_info2
      z$hq  <- aux_info1 + log(log_n) * 2 * aux_info2

    } else {

      aux_info1 <- -2 * z$loglik
      aux_info2 <- p1 + q1 + 2
      log_n <- log(n)

      z$aic <- aux_info1 + 2 * aux_info2
      z$bic <- aux_info1 + log_n * aux_info2
      z$hq  <- aux_info1 + log(log_n) * 2 * aux_info2

    }

    # ---------------------------------------------------------------------
    # Error in the predictor scale
    # ---------------------------------------------------------------------

    # define objects necessary for computing residuals
    fitted_res <- z$fitted[(m + 1):n]
    etahat_res <- z$etahat[(m + 1):n]
    ynew_res <- ynew[(m + 1):n]
    phi_res <- z$phi

    # computing residuals

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

    # final
    resid2_final <- resid2_res_eta / sqrt(resid2_var_res_eta)

    z$resid2 <- resid2_final

  }
  # check inv matrix - end -------------------------------------------- #

  return(z)

}
