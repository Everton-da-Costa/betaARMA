# =========================================================================== #
#      FUNCTIONS FOR THE BETA AUTOREGRESSIVE (BARX) MODEL                     #
# =========================================================================== #
# This file contains the core functions for the BARX model:
#   - loglik_arx: Unpenalized log-likelihood
#   - loglik_arx_ridge: L2-penalized log-likelihood
#   - score_vector_arx: Unpenalized score vector (gradient)
#   - score_vector_arx_ridge: L2-penalized score vector
#   - inf_matrix_arx: Unpenalized Fisher Information Matrix
#   - inf_matrix_arx_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #

#' @title Log-Likelihood for a BARX Model
#' @description This function computes the log-likelihood of the Beta
#'   Autoregressive model with exogenous regressors (BARX).
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) component.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param beta A vector of regression coefficients for exogenous variables.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @importFrom stats dbeta
#'
#' @return The log-likelihood of the BARX model estimators.
#' @keywords internal
loglik_arx <- function(y,
                       ar,
                       X,
                       alpha = 0,
                       varphi = 0,
                       beta = 0,
                       phi = 0,
                       link) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  # Transform the response variable using the link function
  ynew <- linkfun(y)

  # Model components
  p <- max(ar)
  n <- length(y)
  m <- max(p, na.rm = TRUE)

  # Initialize vectors
  eta <- rep(NA, n)
  beta_mat <- as.matrix(beta)

  # Compute eta for t > m
  for (i in (m + 1):n) {
    eta[i] <- alpha + X[i, ] %*% beta_mat +
      (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% beta_mat))
  }

  # Subset for likelihood computation
  mu <- linkinv(eta[(m + 1):n])
  y1 <- y[(m + 1):n]

  # Compute the log-likelihood using the Beta distribution density
  ll_terms <- dbeta(y1, mu * phi, (1 - mu) * phi, log = TRUE)

  # Sum the individual log-likelihood terms to get the total
  sum_ll <- sum(ll_terms)

  return(sum_ll)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BARX Model
#' @description Computes the L2-penalized (Ridge) log-likelihood for a BARX
#'   model.
#'
#' @param y A time series object of data in (0,1).
#' @param ar A numeric vector of positive integers specifying the AR lags.
#' @param X A numeric matrix of external regressors.
#' @param link A character string specifying the link function.
#' @param alpha The numeric intercept term.
#' @param varphi A numeric vector of AR parameters.
#' @param phi The positive precision parameter (not penalized).
#' @param beta A vector of regression coefficients for X (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A single numeric value: the penalized log-likelihood.
#' @keywords internal
loglik_arx_ridge <- function(y,
                             ar,
                             X,
                             link,
                             alpha,
                             varphi,
                             phi,
                             beta,
                             penalty) {
  # Compute the standard (unpenalized) log-likelihood.
  loglik <- loglik_arx(
    y = y,
    ar = ar,
    X = X,
    link = link,
    alpha = alpha,
    varphi = varphi,
    phi = phi,
    beta = beta
  )

  # Determine the number of terms used in BARX maximization
  p_max <- max(ar, na.rm = TRUE)
  a_max <- p_max
  n <- length(y)

  # Define the squared L2 norm of the penalized coefficients (alpha, varphi)
  l2_norm_sq <- alpha^2 + sum(varphi^2)

  # Scale the penalty by the effective number of observations
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Return the penalized log-likelihood
  penalized_loglik <- loglik - penalty * scaled_penalty

  return(penalized_loglik)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BARX Model
#' @description This function computes the score vector of the Beta
#'   Autoregressive model with exogenous regressors (BARX).
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) component.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param beta A vector of regression coefficients for exogenous variables.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @return The score vector of the BARX model estimators.
#' @keywords internal
score_vector_arx <- function(y,
                             ar,
                             X,
                             alpha = 0,
                             varphi = 0,
                             phi = 0,
                             beta = 0,
                             link) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta
  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar)
  p1 <- length(ar)
  k <- length(beta)
  n <- length(y)
  m <- max(p, na.rm = TRUE)

  # Initialize vectors
  eta <- rep(NA, n)
  beta_mat <- as.matrix(beta)

  # Compute linear predictor
  for (i in (m + 1):n) {
    eta[i] <- alpha + X[i, ] %*% beta_mat +
      (varphi %*% (ynew[i - ar] - X[i - ar, , drop = FALSE] %*% beta_mat))
  }

  # Subset for score computation
  eta1 <- eta[(m + 1):n]
  mu1 <- linkinv(eta1)
  y1 <- y[(m + 1):n]
  X1 <- X[(m + 1):n, , drop = FALSE]

  # Intermediate score components
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)
  ystar_mustar <- ystar - mustar

  # Derivatives and design matrices
  mu_eta_val <- mu.eta(eta = eta1)

  # Design matrix P (for autoregressive parameters, varphi)
  P <- matrix(nrow = n - m, ncol = p1)
  for (i in 1:(n - m)) {
    P[i, ] <- ynew[i + m - ar] - X[i + m - ar, , drop = FALSE] %*% beta_mat
  }

  # Design matrix M (for regression coefficients, beta)
  M <- matrix(nrow = n - m, ncol = k)
  for (i in 1:(n - m)) {
    for (j in 1:k) {
      M[i, j] <- X1[i, j] - sum(varphi * X[i + m - ar, j, drop = FALSE])
    }
  }

  # Score vector components
  mT_ystar_mustar <- mu_eta_val * ystar_mustar
  U_alpha <- phi * sum(mT_ystar_mustar)
  U_varphi <- phi * crossprod(P, mT_ystar_mustar)
  U_beta <- phi * crossprod(M, mT_ystar_mustar)
  U_phi <- sum(mu1 * ystar_mustar + log(1 - y1) -
                 digamma((1 - mu1) * phi) + digamma(phi))

  # Return complete score vector in order: alpha, varphi, phi, beta
  score_vec <- c(U_alpha, U_varphi, U_phi, U_beta)
  return(score_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BARX Model
#' @description Computes the L2-penalized score vector for a BARX model.
#'
#' @param y Time series data in (0,1).
#' @param ar Vector of AR lags.
#' @param X Matrix of external regressors.
#' @param link Link function string.
#' @param alpha Intercept term.
#' @param varphi Vector of AR parameters.
#' @param phi Precision parameter (not penalized).
#' @param beta Vector of regression coefficients (not penalized).
#' @param penalty L2 penalty parameter (lambda).
#'
#' @return A numeric vector for the L2-penalized score. The order is:
#'   (alpha, varphi, phi, beta).
#' @keywords internal
score_vector_arx_ridge <- function(y, ar, X, link,
                                   alpha, varphi, phi, beta,
                                   penalty) {
  # Unpenalized score vector from the base function.
  unpenalized_score <- score_vector_arx(
    y = y, ar = ar, X = X, link = link, alpha = alpha,
    varphi = varphi, phi = phi, beta = beta
  )

  # Determine maximum lag order and sample size
  p_max <- max(ar, na.rm = TRUE)
  a_max <- p_max
  n <- length(y)

  # Gradient of the penalty term, matching score vector order:
  # (alpha, varphi, phi, beta)
  beta_gradient <- rep(0, length(beta))
  l2_norm_sq_gradient <- c(2 * alpha, 2 * varphi, 0, beta_gradient)

  # Scale and subtract the penalty gradient
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BARX Model
#' @description This function computes the Fisher information matrix of the
#'   Beta Autoregressive model with exogenous regressors (BetaARX).
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#' @param X A matrix of external regressors.
#' @param alpha The intercept term of the linear predictor.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#' @param phi The precision parameter of the Beta distribution.
#' @param beta A numeric vector of parameters for the external regressors.
#' @param link The link function to be used.
#'
#' @importFrom stats ts start frequency
#'
#' @return A list with the following elements:
#'   - 'fisher_info_mat': The Fisher information matrix.
#'   - 'fitted': The fitted values (conditional mean).
#'   - 'muhat': The truncated vector of fitted conditional means.
#'   - 'etahat': The full vector of the fitted linear predictor.
#'   - 'errorhat': The full vector of errors on the predictor scale.
#' @keywords internal
inf_matrix_arx <- function(y,
                           ar,
                           X,
                           alpha = 0,
                           varphi = 0,
                           phi = 0,
                           beta = 0,
                           link) {
  # Link functions
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta
  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar)
  p1 <- length(ar)
  k1 <- length(beta)
  n <- length(y)
  m <- max(p, na.rm = TRUE)

  # Initialize vectors for fitted values and errors
  errorhat <- rep(0, n)
  etahat <- rep(NA, n)

  # Compute linear predictor (etahat) and errors (errorhat) recursively
  for (i in (m + 1):n) {
    etahat[i] <- alpha +
      X[i, , drop = FALSE] %*% as.matrix(beta) +
      (varphi %*% (ynew[i - ar] -
                     X[i - ar, , drop = FALSE] %*% as.matrix(beta)))
    errorhat[i] <- ynew[i] - etahat[i]
  }

  # Truncate etahat and muhat to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # P: Matrix of lagged (ynew - X*beta) values
  P <- matrix(rep(NA, (n - m) * p1), ncol = p1)
  for (i in 1:(n - m)) {
    P[i, ] <- ynew[i + m - ar] -
      X[i + m - ar, , drop = FALSE] %*% as.matrix(beta)
  }

  # M: Matrix for derivatives with respect to beta
  M <- matrix(rep(NA, (n - m) * k1), ncol = k1)
  for (i in 1:(n - m)) {
    for (j in 1:k1) {
      M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
    }
  }

  # Vector of ones for derivative w.r.t. alpha
  vI <- as.vector(rep(1, n - m))

  # Components for Fisher Information Matrix calculation
  one_minus_muhat <- 1 - muhat
  psi1 <- trigamma(muhat * phi)
  psi2 <- trigamma(one_minus_muhat * phi)
  mu_eta_deriv <- mu.eta(etahat_obs)
  mT <- diag(mu_eta_deriv)
  W_diag_elements <- phi * (psi1 + psi2)
  W <- diag(W_diag_elements) %*% mT^2
  vc <- phi * (psi1 * muhat - psi2 * one_minus_muhat)
  D <- diag(psi1 * (muhat^2) + psi2 * (one_minus_muhat^2) - trigamma(phi))

  # Pre-calculate weighted derivatives
  W_vI <- diag(W) * vI
  W_P <- diag(W) * P
  W_M <- diag(W) * M
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  Kaa <- as.matrix(phi * sum(diag(W)))
  Kpa <- phi * crossprod(P, W_vI)
  Kap <- t(Kpa)
  Kaphi <- crossprod(vI, mT_vc)
  Kphia <- t(Kaphi)
  Kpp <- phi * crossprod(P, W_P)
  Kpphi <- crossprod(P, mT_vc)
  Kphip <- t(Kpphi)
  Kphiphi <- sum(diag(D))
  Kba <- phi * crossprod(M, W_vI)
  Kbb <- phi * crossprod(M, W_M)
  Kbphi <- crossprod(M, mT_vc)
  Kbp <- phi * crossprod(M, W_P)
  Kab <- t(Kba)
  Kphib <- t(Kbphi)
  Kpb <- t(Kbp)

  # Assemble the full Fisher Information Matrix
  # Order: alpha, varphi, phi, beta
  fisher_info_mat <- rbind(
    cbind(Kaa, Kap, Kaphi, Kab),
    cbind(Kpa, Kpp, Kpphi, Kpb),
    cbind(Kphia, Kphip, Kphiphi, Kphib),
    cbind(Kba, Kbp, Kbphi, Kbb)
  )

  # Set matrix row and column names
  names_varphi <- paste0("varphi", ar)
  names_beta <- colnames(X)
  if (is.null(names_beta)) {
    names_beta <- paste0("beta", 1:k1)
  }
  all_names <- c("alpha", names_varphi, "phi", names_beta)
  colnames(fisher_info_mat) <- all_names
  rownames(fisher_info_mat) <- all_names

  # Prepare fitted values as a time series object
  fitted_values <- ts(c(rep(NA, m), muhat),
                      start = start(y), frequency = frequency(y))

  # Prepare output list
  output_list <- list(
    fisher_info_mat = fisher_info_mat,
    fitted = fitted_values,
    muhat = muhat,
    etahat = etahat,
    errorhat = errorhat
  )

  return(output_list)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Fisher Information Matrix for a BARX Model
#' @description This function computes the L2-penalized Fisher Information
#'   Matrix (FIM) for a Beta Autoregressive model with external regressors.
#'
#' @param y A time series of data in the (0,1) interval.
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param X A matrix of external regressors.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param phi The precision parameter (not penalized).
#' @param beta A vector of regression coefficients (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_arx_ridge <- function(y,
                                 ar,
                                 X,
                                 link,
                                 alpha,
                                 varphi,
                                 phi,
                                 beta,
                                 penalty) {
  # Compute the unpenalized Fisher Information Matrix.
  fim_list <- inf_matrix_arx(
    y = y, ar = ar, X = X, link = link, alpha = alpha,
    varphi = varphi, phi = phi, beta = beta
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # Determine the number of terms used in maximization
  p_max <- max(ar, na.rm = TRUE)
  a_max <- p_max
  n <- length(y)

  # Construct the penalty matrix (Hessian of the penalty term).
  # Parameter order: (alpha, varphi, phi, beta)
  p_len <- length(varphi)
  beta_len <- length(beta)

  # Penalty applies to alpha and varphi. No penalty for phi and beta.
  penalty_diag <- c(1, rep(1, p_len), 0, rep(0, beta_len))
  penalty_matrix <- 2 * diag(penalty_diag)

  # Compute the penalized Fisher Information Matrix.
  penalized_fim <- unpenalized_fim + penalty * (n - a_max) * penalty_matrix

  # Prepare the output list.
  output_list <- list(
    fisher_info_mat_penalized = penalized_fim,
    fitted = fim_list$fitted,
    muhat = fim_list$muhat,
    etahat = fim_list$etahat,
    errorhat = fim_list$errorhat
  )

  return(output_list)
}
