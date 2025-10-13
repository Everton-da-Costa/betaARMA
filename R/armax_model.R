# =========================================================================== #
#     FUNCTIONS FOR THE BETA AUTOREGRESSIVE MOVING AVERAGE (BARMAX) MODEL     #
# =========================================================================== #
# This file contains the core functions for the BARMAX model:
#   - loglik_armax: Unpenalized log-likelihood
#   - loglik_armax_ridge: L2-penalized log-likelihood
#   - score_vector_armax: Unpenalized score vector (gradient)
#   - score_vector_armax_ridge: L2-penalized score vector
#   - inf_matrix_armax: Unpenalized Fisher Information Matrix
#   - inf_matrix_armax_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #

#' @title Log-Likelihood for a BARMAX Model
#' @description This function computes the log-likelihood of the Beta
#'   Autoregressive Moving Average model with exogenous regressors (BARMAX).
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter of the Beta distribution.
#' @param beta A vector of regression coefficients for exogenous variables.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @importFrom stats dbeta
#'
#' @return The log-likelihood of the BARMAX model estimators.
#' @keywords internal
loglik_armax <- function(y,
                         ar, ma,
                         X, link,
                         alpha = 0, varphi = 0, theta = 0,
                         phi = 0, beta = 0) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv

  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # Initialize vectors
  error <- rep(0, n)
  eta <- rep(NA, n)
  beta_mat <- as.matrix(beta)

  # Compute eta and errors recursively
  for (i in (m + 1):n) {
    eta[i] <- alpha +
      X[i, ] %*% beta_mat +
      varphi %*% (ynew[i - ar] - X[i - ar, ] %*% beta_mat) +
      theta %*% error[i - ma]

    error[i] <- ynew[i] - eta[i]
  }

  # Subset for likelihood computation
  mu <- linkinv(eta[(m + 1):n])
  y1 <- y[(m + 1):n]

  # Compute log-likelihood (Note: warnings are suppressed)
  ll_terms <- suppressWarnings(dbeta(y1, mu * phi, (1 - mu) * phi, log = TRUE))

  return(sum(ll_terms))
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BARMAX Model
#' @description Computes the L2-penalized (Ridge) log-likelihood for a BARMAX
#'   model.
#'
#' @details The penalty is applied to the intercept (`alpha`), AR (`varphi`),
#'   and MA (`theta`) coefficients. The regressor (`beta`) and precision
#'   (`phi`) parameters are not penalized.
#'
#' @param y A time series of numbers in the interval (0,1).
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param X A matrix of external regressors.
#' @param link The link function to be used ("logit", "probit", "cloglog").
#' @param alpha The intercept term of the linear predictor.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter of the beta distribution.
#' @param beta A vector of regression coefficients for external regressors.
#' @param penalty A non-negative scalar for the L2 penalty.
#'
#' @return The value of the L2 penalized log-likelihood.
#' @keywords internal
loglik_armax_ridge <- function(y,
                               ar, ma,
                               X, link,
                               alpha,
                               varphi,
                               theta,
                               phi,
                               beta,
                               penalty) {
  # Compute the standard (unpenalized) log-likelihood.
  loglik <- loglik_armax(y,
                         ar, ma,
                         X = X,
                         alpha = alpha,
                         varphi = varphi,
                         theta = theta,
                         phi = phi,
                         link = link,
                         beta = beta)

  # Determine the number of terms used in BARMAX maximization
  p_max <- max(ar, na.rm = TRUE)
  q_max <- max(ma, na.rm = TRUE)
  a_max <- max(p_max, q_max)
  n <- length(y)

  # Calculate the penalty term.
  l2_norm_sq <- alpha^2 + sum(varphi^2) + sum(theta^2)
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Return the penalized log-likelihood
  ell2_penalty <- loglik - penalty * scaled_penalty

  return(ell2_penalty)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BARMAX Model
#' @description This function computes the score vector of the Beta
#'   Autoregressive Moving Average model with exogenous regressors (BARMAX).
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter of the Beta distribution.
#' @param beta A vector of regression coefficients for exogenous variables.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @return The score vector of the BARMAX estimators.
#' @keywords internal
score_vector_armax <- function(y,
                               ar, ma,
                               X, link,
                               alpha = 0,
                               varphi = 0,
                               theta = 0,
                               phi = 0,
                               beta = 0) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar)
  q <- max(ma)
  p1 <- length(ar)
  q1 <- length(ma)
  k1 <- length(beta)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # Initialize vectors
  error <- rep(0, n)
  eta <- rep(NA, n)

  # Compute eta and errors recursively
  for (i in (m + 1):n) {
    eta[i] <- alpha +
      X[i, ] %*% as.matrix(beta) +
      (varphi %*% (ynew[i - ar] - X[i - ar, ] %*% as.matrix(beta))) +
      (theta %*% error[i - ma])
    error[i] <- ynew[i] - eta[i]
  }

  # Subset for score computation
  eta1 <- eta[(m + 1):n]
  mu <- linkinv(eta1)
  y1 <- y[(m + 1):n]
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  # Derivative and lagged matrices
  mT <- diag(mu.eta(eta1))
  R <- matrix(NA, nrow = n - m, ncol = q1)
  for (i in 1:(n - m)) R[i, ] <- error[i + m - ma]

  P <- matrix(NA, nrow = n - m, ncol = p1)
  for (i in 1:(n - m)) {
    P[i, ] <- ynew[i + m - ar] - X[i + m - ar, ] %*% as.matrix(beta)
  }

  M <- matrix(NA, nrow = n - m, ncol = k1)
  for (i in 1:(n - m)) {
    for (j in 1:k1) {
      M[i, j] <- X[i + m, j] - sum(varphi * X[i + m - ar, j])
    }
  }

  # Recursive calculation of eta derivatives
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

  # Score vector components
  ystar_mustar <- ystar - mustar
  mT_ystar_mustar <- mT %*% ystar_mustar

  U_alpha  <- phi * s     %*% mT_ystar_mustar
  U_varphi <- phi * t(rP) %*% mT_ystar_mustar
  U_theta  <- phi * t(rR) %*% mT_ystar_mustar
  U_beta   <- phi * t(rM) %*% mT_ystar_mustar

  U_phi <- sum(mu * ystar_mustar + log(1 - y1)
               - digamma((1 - mu) * phi) + digamma(phi))

  # Return complete score vector
  score_vec <- c(U_alpha, U_varphi, U_theta, U_phi, U_beta)

  return(score_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BARMAX Model
#' @description Computes the L2-penalized score vector for a BARMAX model.
#'
#' @param y Time series data in (0,1).
#' @param ar Vector of AR lags.
#' @param ma Vector of MA lags.
#' @param X Matrix of external regressors.
#' @param link Link function string.
#' @param alpha Intercept term.
#' @param varphi Vector of AR parameters.
#' @param theta Vector of MA parameters.
#' @param phi Precision parameter (not penalized).
#' @param beta Vector of regression coefficients (not penalized).
#' @param penalty L2 penalty parameter (lambda).
#'
#' @return A numeric vector for the L2-penalized score. The order is:
#'   (alpha, varphi, theta, phi, beta).
#' @keywords internal
score_vector_armax_ridge <- function(y,
                                     ar,
                                     ma,
                                     X,
                                     link,
                                     alpha,
                                     varphi,
                                     theta,
                                     phi,
                                     beta,
                                     penalty) {
  # Compute the standard (unpenalized) score vector.
  unpenalized_score <- score_vector_armax(
    y = y, ar = ar, ma = ma, X = X, link = link, alpha = alpha,
    varphi = varphi, theta = theta, phi = phi, beta = beta
  )

  # Determine maximum lag order and sample size.
  p_max <- max(ar, na.rm = TRUE)
  q_max <- max(ma, na.rm = TRUE)
  a_max <- max(p_max, q_max)
  n <- length(y)

  # Gradient of the penalty term.
  # Order must match score vector: (alpha, varphi, theta, phi, beta).
  beta_gradient <- rep(0, length(beta))
  l2_norm_sq_gradient <- c(2 * alpha, 2 * varphi, 2 * theta, 0, beta_gradient)

  # Scale and subtract the penalty gradient.
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BARMAX Model
#' @description This function computes the Fisher information matrix for the
#'   BARMAX model.
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#' @param X A matrix of exogenous variables (regressors).
#' @param alpha The intercept term of the linear predictor.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#' @param theta A numeric vector of moving average (MA) parameters.
#' @param beta A numeric vector of regression coefficients for X.
#' @param phi The precision parameter of the Beta distribution.
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
inf_matrix_armax <- function(y,
                             ar,
                             ma,
                             X,
                             alpha = 0,
                             varphi = 0,
                             theta = 0,
                             phi = 0,
                             beta = 0,
                             link = "logit") {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar); q <- max(ma)
  p1 <- length(ar); q1 <- length(ma); k1 <- length(beta)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # Initialize vectors
  errorhat <- rep(0, n)
  etahat <- rep(NA, n)

  # Compute linear predictor (etahat) and errors (errorhat) recursively
  for (i in (m + 1):n) {
    etahat[i] <- alpha +
      X[i, , drop = FALSE] %*% as.matrix(beta) +
      (varphi %*% (ynew[i - ar] -
                     X[i - ar, , drop = FALSE] %*% as.matrix(beta))) +
      theta %*% errorhat[i - ma]
    errorhat[i] <- ynew[i] - etahat[i]
  }

  # Truncate to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # Prepare lagged matrices for derivative calculations
  R <- matrix(NA, (n - m), q1); for (i in 1:(n-m)) R[i,] <- errorhat[i+m-ma]
  P <- matrix(NA, (n-m), p1); for (i in 1:(n-m)) P[i,] <- ynew[i+m-ar] -
    X[i+m-ar, , drop = FALSE] %*% as.matrix(beta)
  M <- matrix(NA, (n-m), k1); for (i in 1:(n-m)) for (j in 1:k1) {
    M[i,j] <- X[i+m,j] - sum(varphi * X[i+m-ar,j])}

  # Recursive calculation of eta derivatives
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

  # Components for Fisher Information Matrix
  one_minus_muhat <- 1 - muhat
  psi1 <- trigamma(muhat * phi)
  psi2 <- trigamma(one_minus_muhat * phi)
  mu_eta_deriv <- mu.eta(etahat_obs)
  mT <- diag(mu_eta_deriv)

  W_diag_elements <- phi * (psi1 + psi2)
  W <- diag(W_diag_elements) %*% mT^2
  vc <- phi * (psi1 * muhat - psi2 * one_minus_muhat)
  D <- diag(psi1*(muhat^2) + psi2*(one_minus_muhat^2) - trigamma(phi))

  # Pre-calculate weighted derivatives
  W_a <- diag(W) * a; W_rP <- diag(W) * rP
  W_rR <- diag(W) * rR; W_rM <- diag(W) * rM
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  K_a_a <- phi * crossprod(a, W_a); K_p_a <- phi * crossprod(rP, W_a)
  K_t_a <- phi * crossprod(rR, W_a); K_phi_a <- crossprod(a, mT_vc)
  K_p_p <- phi * crossprod(rP, W_rP); K_p_t <- phi * crossprod(rP, W_rR)
  K_p_phi <- crossprod(rP, mT_vc); K_t_t <- phi * crossprod(rR, W_rR)
  K_t_phi <- crossprod(rR, mT_vc); K_phi_phi <- sum(diag(D))
  K_b_a <- phi * crossprod(rM, W_a); K_b_b <- phi * crossprod(rM, W_rM)
  K_b_phi <- crossprod(rM, mT_vc); K_b_p <- phi * crossprod(rM, W_rP)
  K_b_t <- phi * crossprod(rM, W_rR)

  # Assemble the full Fisher Information Matrix
  fisher_info_mat <- rbind(
    cbind(K_a_a, t(K_p_a), t(K_t_a), t(K_phi_a), t(K_b_a)),
    cbind(K_p_a, K_p_p, K_p_t, K_p_phi, t(K_b_p)),
    cbind(K_t_a, t(K_p_t), K_t_t, K_t_phi, t(K_b_t)),
    cbind(K_phi_a, t(K_p_phi), t(K_t_phi), K_phi_phi, t(K_b_phi)),
    cbind(K_b_a, K_b_p, K_b_t, K_b_phi, K_b_b)
  )

  # Set matrix row and column names
  names_varphi <- paste0("varphi", ar)
  names_theta  <- paste0("theta", ma)
  names_beta <- colnames(X)
  if (is.null(names_beta)) names_beta <- paste0("beta", 1:k1)
  param_names <- c("alpha", names_varphi, names_theta, "phi", names_beta)
  colnames(fisher_info_mat) <- param_names
  rownames(fisher_info_mat) <- param_names

  # Prepare output list
  fitted_values <- ts(c(rep(NA, m), muhat),
                      start = start(y), frequency = frequency(y))
  output_list <- list(
    fisher_info_mat = fisher_info_mat,
    fitted = fitted_values,
    muhat = muhat,
    etahat = etahat, errorhat = errorhat
    )

  return(output_list)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Fisher Information Matrix for a BARMAX Model
#' @description This function computes the L2-penalized FIM for a BARMAX model.
#'
#' @param y A time series of numbers in the interval (0,1).
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param X A matrix of external regressors.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter (not penalized).
#' @param beta A vector of regression coefficients (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_armax_ridge <- function(y,
                                   ar,
                                   ma,
                                   X,
                                   link,
                                   alpha,
                                   varphi,
                                   theta,
                                   phi,
                                   beta,
                                   penalty) {
  # Compute the unpenalized Fisher Information Matrix.
  fim_list <- inf_matrix_armax(
    y = y, ar = ar, ma = ma, X = X, link = link, alpha = alpha,
    varphi = varphi, theta = theta, phi = phi, beta = beta
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # Determine maximum lag and sample size
  p_max <- max(ar, na.rm = TRUE); q_max <- max(ma, na.rm = TRUE)
  a_max <- max(p_max, q_max); n <- length(y)

  # Construct the penalty matrix (Hessian of the penalty term).
  # Parameter order: (alpha, varphi, theta, phi, beta)
  p_len <- length(varphi); q_len <- length(theta); beta_len <- length(beta)

  # Penalty applies to alpha, varphi, and theta.
  penalty_diag <- c(1, rep(1, p_len), rep(1, q_len), 0, rep(0, beta_len))
  penalty_matrix <- 2 * diag(penalty_diag)

  # Compute the penalized Fisher Information Matrix.
  penalized_fim <- unpenalized_fim + penalty * (n - a_max) * penalty_matrix

  # Prepare the output list (Note: contains a bug).
  output_list <- list()
  output_list$fisher_info_mat_penalized <- penalized_fim
  output_list$fitted <- fim_list$fitted
  output_list$muhat <- fim_list$muhat
  output_list$etahat <- fim_list$etahat
  output_list$errorhat <- fim_list$errorhat

  return(output_list)
}
