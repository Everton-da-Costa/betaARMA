# =========================================================================== #
#      FUNCTIONS FOR THE BETA MOVING AVERAGE (BMAX) MODEL                     #
# =========================================================================== #
# This file contains the core functions for the BMAX model:
#   - loglik_max: Unpenalized log-likelihood
#   - loglik_max_ridge: L2-penalized log-likelihood
#   - score_vector_max: Unpenalized score vector (gradient)
#   - score_vector_max_ridge: L2-penalized score vector
#   - inf_matrix_max: Unpenalized Fisher Information Matrix
#   - inf_matrix_max_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #

#' @title Log-Likelihood for a BMAX Model
#' @description This function computes the conditional log-likelihood of a
#'   Beta regression model with Moving Average (MA) errors and exogenous
#'   regressors (a BMAX model).
#'
#' @param y A time series object of data in (0,1).
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param X A numeric matrix of external regressors.
#' @param link A character string specifying the link function. Valid options
#'   are "logit", "probit", "loglog", or "cloglog".
#' @param alpha The numeric intercept term in the linear predictor.
#' @param theta A numeric vector of MA parameters. Its length must match the
#'   length of the 'ma' vector.
#' @param phi The numeric positive precision parameter.
#' @param beta A vector of regression coefficients for exogenous variables.
#'
#' @importFrom stats dbeta
#'
#' @return A numeric value representing the total conditional log-likelihood.
#' @keywords internal
loglik_max <- function(y,
                       ma,
                       X,
                       link,
                       alpha = 0,
                       theta = 0,
                       phi = 0,
                       beta = 0) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  ynew <- linkfun(y)

  # Model dimensions
  q <- max(ma)
  n <- length(y)
  m <- max(q, na.rm = TRUE)

  # Initialize vectors
  error <- rep(0, n)
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  # Compute eta and mu for t > m
  for (t in (m + 1):n) {
    eta[t] <- alpha + X[t, ] %*% as.matrix(beta) + (theta %*% error[t - ma])
    error[t] <- ynew[t] - eta[t]
  }

  # Truncate for conditional likelihood
  mu1 <- linkinv(eta[(m + 1):n])
  y1 <- y[(m + 1):n]

  # Compute the log-likelihood using the Beta distribution density
  ll_terms <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  # Sum the individual log-likelihood terms
  sum_ll <- sum(ll_terms)

  return(sum_ll)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BMAX Model
#' @description This function computes the L2 penalized log-likelihood of
#'   the Beta Moving Average Regression (BMAR) model.
#'
#' @param y A time series of numbers in the interval (0,1).
#' @param ma A vector specifying the moving average (MA) component.
#' @param X A matrix of external regressors.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term of the linear predictor.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter of the beta distribution.
#' @param beta A vector of regression coefficients for external regressors.
#' @param penalty A non-negative scalar for the L2 penalty.
#'
#' @return The value of the L2 penalized log-likelihood.
#' @keywords internal
loglik_max_ridge <- function(y,
                             ma,
                             X,
                             link,
                             alpha,
                             theta,
                             phi,
                             beta,
                             penalty) {
  # Calculate the standard log-likelihood
  loglik <- loglik_max(
    y,
    ma = ma,
    X = X,
    alpha = alpha,
    theta = theta,
    phi = phi,
    link = link,
    beta = beta
  )

  # Determine number of terms for BMAR maximization
  q_max <- max(ma, na.rm = TRUE)
  a_max <- q_max
  n <- length(y)

  # L2 penalty calculation (applied to alpha and theta)
  l2_norm_sq <- alpha^2 + sum(theta^2)
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Apply the penalty to the log-likelihood
  penalized_loglik <- loglik - penalty * scaled_penalty

  return(penalized_loglik)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BMAX Model
#' @description This function computes the score vector of the Beta Moving
#'   Average model with exogenous regressors (BMAX).
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ma A vector specifying the moving average (MA) component.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term.
#' @param theta A vector of moving average (MA) parameters.
#' @param beta A vector of regression coefficients for exogenous variables.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @return The score vector of the BMAX model estimators.
#' @keywords internal
score_vector_max <- function(y,
                             ma,
                             alpha = 0,
                             theta = 0,
                             phi = 0,
                             beta = 0,
                             X,
                             link) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  q <- max(ma)
  q1 <- length(ma)
  k <- length(beta)
  n <- length(y)
  m <- max(q, na.rm = TRUE)

  # Initialize vectors
  error <- rep(0, n)
  eta <- rep(NA, n)
  beta_mat <- as.matrix(beta)

  # Compute linear predictor and error
  for (i in (m + 1):n) {
    eta[i] <- alpha + X[i, ] %*% beta_mat + (theta %*% error[i - ma])
    error[i] <- ynew[i] - eta[i]
  }

  # Subset for score computation
  eta1 <- eta[(m + 1):n]
  mu1 <- linkinv(eta1)
  y1 <- y[(m + 1):n]

  # Intermediate score components
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)
  ystar_mustar <- ystar - mustar
  mu_eta_val <- mu.eta(eta = eta1)

  # Design matrix R (lagged errors for derivative w.r.t. theta)
  R <- matrix(nrow = n - m, ncol = q1)
  for (i in 1:(n - m)) {
    R[i, ] <- error[i + m - ma]
  }

  # Design matrix M (regressors for derivative w.r.t. beta)
  M <- matrix(nrow = n - m, ncol = k)
  for (i in 1:(n - m)) {
    for (j in 1:k) {
      M[i, j] <- X[i + m, j]
    }
  }

  # Recursive calculation of eta derivatives
  deta_dalpha <- rep(0, n)
  deta_dtheta <- matrix(0, ncol = q1, nrow = n)
  deta_dbeta <- matrix(0, ncol = k, nrow = n)

  for (i in (m + 1):n) {
    deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
    deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
    deta_dbeta[i, ] <- M[(i - m), ] - theta %*% deta_dbeta[i - ma, ]
  }

  # Subset derivatives to the effective sample size
  s <- deta_dalpha[(m + 1):n]
  rR <- deta_dtheta[(m + 1):n, ]
  rM <- deta_dbeta[(m + 1):n, ]

  # Score vector components
  mT_ystar_mustar <- mu_eta_val * ystar_mustar
  U_alpha <- phi * sum(s * mT_ystar_mustar)
  U_theta <- phi * crossprod(rR, mT_ystar_mustar)
  U_beta <- phi * crossprod(rM, mT_ystar_mustar)
  U_phi <- sum(mu1 * ystar_mustar + log(1 - y1) -
                 digamma((1 - mu1) * phi) + digamma(phi))

  # Return complete score vector
  score_vec <- c(U_alpha, U_theta, U_phi, U_beta)
  return(score_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BMAX Model
#' @description Computes the L2-penalized score vector for a BMAX model.
#'
#' @param y Time series data in (0,1).
#' @param ma Vector of MA lags.
#' @param X Matrix of external regressors.
#' @param link Link function string.
#' @param alpha Intercept term.
#' @param theta Vector of MA parameters.
#' @param phi Precision parameter (not penalized).
#' @param beta Vector of regression coefficients (not penalized).
#' @param penalty L2 penalty parameter (lambda).
#'
#' @return A numeric vector for the L2-penalized score. The order is:
#'   (alpha, theta, phi, beta).
#' @keywords internal
score_vector_max_ridge <- function(y, ma, X, link,
                                   alpha, theta, phi, beta,
                                   penalty) {
  # Unpenalized score vector
  unpenalized_score <- score_vector_max(
    y = y, ma = ma, X = X, link = link, alpha = alpha,
    theta = theta, phi = phi, beta = beta
  )

  # Determine max lag order and sample size
  q_max <- max(ma, na.rm = TRUE)
  a_max <- q_max
  n <- length(y)

  # Gradient of the penalty term (matches score vector order)
  beta_gradient <- rep(0, length(beta))
  l2_norm_sq_gradient <- c(2 * alpha, 2 * theta, 0, beta_gradient)

  # Scale and subtract the penalty gradient
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BMAX Model
#' @description This function computes the Fisher information matrix for the
#'   Beta Moving Average model with exogenous regressors (BMAX).
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ma A numeric vector specifying the moving average (MA) lags.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term of the linear predictor.
#' @param theta A numeric vector of moving average (MA) parameters.
#' @param beta A numeric vector of regression coefficients.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function to be used ("logit", "probit", etc.).
#'
#' @importFrom stats ts start frequency
#'
#' @return A list with the following elements:
#'   - 'fisher_info_mat': The Fisher information matrix.
#'   - 'fitted': A time series object of the fitted conditional mean values.
#'   - 'muhat': The truncated vector of fitted conditional mean values.
#'   - 'etahat': The full vector of the linear predictor values.
#'   - 'errorhat': The full vector of errors on the predictor scale.
#' @keywords internal
inf_matrix_max <- function(y,
                           ma,
                           X,
                           alpha = 0,
                           theta = 0,
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
  q <- max(ma)
  q1 <- length(ma)
  k1 <- length(beta)
  n <- length(y)
  m <- max(q, na.rm = TRUE)

  # Initialize vectors for fitted values and errors
  errorhat <- rep(0, n)
  etahat <- rep(NA, n)

  # Compute linear predictor (etahat) and errors (errorhat) recursively
  for (i in (m + 1):n) {
    etahat[i] <- alpha +
      X[i, , drop = FALSE] %*% as.matrix(beta) +
      theta %*% errorhat[i - ma]
    errorhat[i] <- ynew[i] - etahat[i]
  }

  # Truncate to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # R: Matrix of lagged errorhat values
  R <- matrix(rep(NA, (n - m) * q1), ncol = q1)
  for (i in 1:(n - m)) {
    R[i, ] <- errorhat[i + m - ma]
  }

  # M: Matrix for derivatives with respect to beta
  M <- matrix(rep(NA, (n - m) * k1), ncol = k1)
  for (i in 1:(n - m)) {
    for (j in 1:k1) {
      M[i, j] <- X[i + m, j]
    }
  }

  # Recursive calculation of eta derivatives
  deta_dalpha <- rep(0, n)
  deta_dtheta <- matrix(0, ncol = q1, nrow = n)
  deta.dbeta <- matrix(0, ncol = k1, nrow = n)

  for (i in (m + 1):n) {
    deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
    deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
    deta.dbeta[i, ] <- M[(i - m), ] - theta %*% deta.dbeta[i - ma, ]
  }

  # Truncate derivatives to the observed range
  a <- deta_dalpha[(m + 1):n]
  rR <- deta_dtheta[(m + 1):n, ]
  rM <- deta.dbeta[(m + 1):n, ]

  # Components for Fisher Information Matrix
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
  W_a <- diag(W) * a
  W_rR <- diag(W) * rR
  W_rM <- diag(W) * rM
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  Kaa <- phi * crossprod(a, W_a)
  Kta <- phi * crossprod(rR, W_a)
  Kat <- t(Kta)
  Kaphi <- crossprod(a, mT_vc)
  Kphia <- t(Kaphi)
  Ktt <- phi * crossprod(rR, W_rR)
  Ktphi <- crossprod(rR, mT_vc)
  Kphit <- t(Ktphi)
  Kphiphi <- sum(diag(D))
  Kba <- phi * crossprod(rM, W_a)
  Kbb <- phi * crossprod(rM, W_rM)
  Kbphi <- crossprod(rM, mT_vc)
  Kbt <- phi * crossprod(rM, W_rR)
  Kab <- t(Kba)
  Kphib <- t(Kbphi)
  Ktb <- t(Kbt)

  # Assemble the full Fisher Information Matrix
  fisher_info_mat <- rbind(
    cbind(Kaa, Kat, Kaphi, Kab),
    cbind(Kta, Ktt, Ktphi, Ktb),
    cbind(Kphia, Kphit, Kphiphi, Kphib),
    cbind(Kba, Kbt, Kbphi, Kbb)
  )

  # Set matrix row and column names
  names_theta  <- paste("theta", ma, sep = "")
  names_beta <- colnames(X)
  if (is.null(names_beta)) {
    names_beta <- paste0("beta", 1:k1)
  }
  names_fisher <- c("alpha", names_theta, "phi", names_beta)
  colnames(fisher_info_mat) <- names_fisher
  rownames(fisher_info_mat) <- names_fisher

  # Prepare output
  fitted_values <- ts(c(rep(NA, m), muhat),
                      start = start(y), frequency = frequency(y))
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

#' @title Penalized Fisher Information Matrix for a BMAX Model
#' @description This function computes the L2-penalized Fisher Information
#'   Matrix (FIM) for a BMAX model. The penalty is applied to the intercept
#'   (alpha) and the MA parameters (theta).
#'
#' @param y A time series of data in the (0,1) interval.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param X A matrix of external regressors.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter (not penalized).
#' @param beta A vector of regression coefficients (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_max_ridge <- function(y,
                                 ma,
                                 X,
                                 link,
                                 alpha,
                                 theta,
                                 phi,
                                 beta,
                                 penalty) {
  # Compute the unpenalized Fisher Information Matrix
  fim_list <- inf_matrix_max(
    y = y, ma = ma, X = X, link = link, alpha = alpha,
    theta = theta, phi = phi, beta = beta
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # Determine number of terms for BMAX maximization
  q_max <- max(ma, na.rm = TRUE)
  a_max <- q_max
  n <- length(y)

  # Construct the penalty matrix (Hessian of the penalty term)
  q_len <- length(theta)
  beta_len <- length(beta)
  penalty_diag <- c(1, rep(1, q_len), 0, rep(0, beta_len))
  penalty_matrix <- 2 * diag(penalty_diag)

  # Compute the penalized Fisher Information Matrix
  penalized_fim <- unpenalized_fim + penalty * (n - a_max) * penalty_matrix

  # Prepare the output list
  output_list <- list()
  output_list$fisher_info_mat_penalized <- penalized_fim
  output_list$fitted <- fim_list$fitted
  output_list$muhat <- fim_list$muhat
  output_list$etahat <- fim_list$etahat
  output_list$errorhat <- fim_list$errorhat

  return(output_list)
}
