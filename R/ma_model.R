# =========================================================================== #
#        FUNCTIONS FOR THE BETA MOVING AVERAGE (BMA) MODEL                    #
# =========================================================================== #
# This file contains the core functions for the BMA model:
#   - loglik_ma: Unpenalized log-likelihood
#   - loglik_ma_ridge: L2-penalized log-likelihood
#   - score_vector_ma: Unpenalized score vector (gradient)
#   - score_vector_ma_ridge: L2-penalized score vector
#   - inf_matrix_ma: Unpenalized Fisher Information Matrix
#   - inf_matrix_ma_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #

#' @title Log-Likelihood for a BMA Model
#' @description This function computes the conditional log-likelihood of a
#'   Beta Moving Average (BMA) model.
#'
#' @param y A time series object of data in (0,1).
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param alpha The numeric intercept term in the linear predictor.
#' @param theta A numeric vector of MA parameters. Its length must match the
#'   length of the 'ma' vector.
#' @param phi The numeric positive precision parameter.
#' @param link A character string specifying the link function ("logit",
#'   "probit", "loglog", or "cloglog").
#'
#' @importFrom stats dbeta
#'
#' @return A numeric value representing the total conditional log-likelihood
#'   of the BMA model for the given data and parameters.
#' @keywords internal
loglik_ma <- function(y,
                       ma,
                       alpha = 0,
                       theta = 0,
                       phi = 0,
                       link) {
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
    eta[t] <- alpha + theta %*% error[t - ma]
    mu[t] <- linkinv(eta[t])
    error[t] <- ynew[t] - eta[t]
  }

  # Subset for likelihood computation
  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # Compute the log-likelihood using the Beta distribution density
  ll_terms_ma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  # Sum the individual log-likelihood terms to get the total
  sum_ll <- sum(ll_terms_ma)

  return(sum_ll)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BMA Model
#' @description Computes the L2-penalized (Ridge) log-likelihood for a BMA
#'   model.
#'
#' @details The penalty is applied to the intercept (`alpha`) and MA (`theta`)
#'   coefficients. The precision parameter (`phi`) is not penalized.
#'
#' @param y A numeric vector or ts object of data in (0,1).
#' @param ma A vector of positive integer MA lags.
#' @param alpha The numeric intercept term.
#' @param theta Vector of MA parameters, matching the length of `ma`.
#' @param phi The positive precision parameter.
#' @param link The link function (e.g., "logit", "probit").
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A single numeric value: the penalized log-likelihood.
#' @keywords internal
loglik_ma_ridge <- function(y,
                             ma,
                             alpha,
                             theta,
                             phi,
                             link,
                             penalty) {
  # Compute the standard (unpenalized) log-likelihood.
  loglik <- loglik_ma(
    y,
    ma = ma,
    alpha = alpha,
    theta = theta,
    phi = phi,
    link = link
  )

  # Determine the number of terms used in BMA maximization
  q_max <- max(ma)
  a_max <- max(q_max)
  n <- length(y)

  # Calculate the penalty term.
  l2_norm_sq <- alpha^2 + sum(theta^2)
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Return the penalized log-likelihood
  ell2_penalty <- loglik - penalty * scaled_penalty

  return(ell2_penalty)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BMA Model
#' @description This function computes the score vector of the Beta Moving
#'   Average (BMA) model.
#'
#' @param y A time series object of data in (0,1).
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param alpha The numeric intercept term in the linear predictor.
#' @param theta A numeric vector of MA parameters.
#' @param phi The numeric positive precision parameter.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @return The score vector of the BMA estimators.
#' @keywords internal
score_vector_ma <- function(y,
                            ma,
                            alpha = 0,
                            theta = 0,
                            phi = 0,
                            link) {
  # Link functions setup
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  q <- max(ma)
  q_len <- length(theta)
  n <- length(y)
  m <- max(q, na.rm = TRUE)

  # Initialize vectors
  error <- rep(0, n)
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  # Compute eta, mu, and errors recursively
  for (t in (m + 1):n) {
    eta[t] <- alpha + theta %*% error[t - ma]
    mu[t] <- linkinv(eta[t])
    error[t] <- ynew[t] - eta[t]
  }

  # Subset for score computation
  eta1 <- eta[(m + 1):n]
  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # Score components
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

  # Matrix of lagged errors
  R <- matrix(nrow = n - m, ncol = q_len)
  for (i in 1:(n - m)) R[i, ] <- error[i + m - ma]

  # Recursive calculation of eta derivatives
  deta.dalpha <- rep(0, n)
  deta.dtheta <- matrix(0, ncol = q_len, nrow = n)

  for (i in (m + 1):n) {
    deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i - ma]
    deta.dtheta[i, ] <- R[(i - m), ] - theta %*% deta.dtheta[i - ma, ]
  }

  s <- deta.dalpha[(m + 1):n]
  rR <- deta.dtheta[(m + 1):n, ]

  # Score vector components
  mu_eta <- mu.eta(eta = eta1)
  ystar_mustar <- ystar - mustar
  mT_ystar_mustar <- mu_eta * ystar_mustar

  U_alpha <- phi * crossprod(s, mT_ystar_mustar)
  U_theta <- phi * crossprod(rR, mT_ystar_mustar)
  U_phi <- sum(mu1 * ystar_mustar + log(1 - y1) -
                 digamma((1 - mu1) * phi) + digamma(phi))

  # Return complete score vector
  escore_vec <- c(U_alpha, U_theta, U_phi)

  return(escore_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BMA Model
#' @description This function computes the L2-penalized (Ridge) score vector
#'   for a Beta Moving Average (BMA) model.
#'
#' @param y A time series object of data in (0,1).
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param alpha The numeric intercept term in the linear predictor.
#' @param theta A numeric vector of MA parameters.
#' @param phi The numeric positive precision parameter (not penalized).
#' @param link A character string specifying the link function.
#' @param penalty The numeric regularization parameter (lambda).
#'
#' @return A numeric vector representing the L2-penalized score.
#' @keywords internal
score_vector_ma_ridge <- function(y,
                                  ma,
                                  alpha,
                                  theta,
                                  phi,
                                  link,
                                  penalty) {
  # Compute the standard (unpenalized) score vector.
  unpenalized_score <- score_vector_ma(
    y = y,
    ma = ma,
    alpha = alpha,
    theta = theta,
    phi = phi,
    link = link
  )

  # Determine the number of terms used in BMA maximization
  q_max <- max(ma, na.rm = TRUE)
  a_max <- q_max
  n <- length(y)

  # Compute the gradient of the penalty term.
  l2_norm_sq_gradient <- c(2 * alpha, 2 * theta, 0)

  # Scale the penalty gradient.
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient

  # Compute the final penalized score.
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BMA Model
#' @description This function computes the Fisher information matrix of the
#'   Beta Moving Average (BMA) model.
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ma A numeric vector specifying the moving average (MA) lags.
#' @param alpha The intercept term.
#' @param theta A numeric vector of moving average (MA) parameters.
#' @param phi The precision parameter of the BMA model.
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
inf_matrix_ma <- function(y,
                          ma,
                          alpha = 0,
                          theta = 0,
                          phi = 0,
                          link) {
  # Link functions
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  q <- max(ma)
  q1 <- length(ma)
  n <- length(y)
  m <- max(q, na.rm = TRUE)

  # Initialize vectors for fitted values and errors
  errorhat <- rep(0, n)
  etahat <- rep(NA, n)

  # Compute fitted values and errors (errorhat) recursively
  for (t in (m + 1):n) {
    etahat[t] <- alpha + theta %*% errorhat[t - ma]
    errorhat[t] <- ynew[t] - etahat[t]
  }

  # Truncate etahat and muhat to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # R: Matrix of lagged errorhat values corresponding to MA lags
  R <- matrix(nrow = n - m, ncol = q1)
  for (i in 1:(n - m)) R[i, ] <- errorhat[i + m - ma]

  # Recursive calculation of eta derivatives
  deta_dalpha <- rep(0, n)
  deta_dtheta <- matrix(0, ncol = q1, nrow = n)

  for (i in (m + 1):n) {
    deta_dalpha[i] <- 1 - theta %*% deta_dalpha[i - ma]
    deta_dtheta[i, ] <- R[(i - m), ] - theta %*% deta_dtheta[i - ma, ]
  }

  s <- deta_dalpha[(m + 1):n]
  rR <- deta_dtheta[(m + 1):n, ]

  # Components for Fisher Information Matrix calculation
  one_minus_muhat <- 1 - muhat
  psi1 <- trigamma(muhat * phi)
  psi2 <- trigamma(one_minus_muhat * phi)

  mu_eta_deriv <- mu.eta(eta = etahat_obs)
  mT <- diag(mu_eta_deriv)

  W_diag_elements <- phi * (psi1 + psi2)
  W <- diag(W_diag_elements) %*% mT^2

  vc <- phi * (psi1 * muhat - psi2 * one_minus_muhat)
  D <- diag(psi1 * (muhat^2) + psi2 * (one_minus_muhat^2) - trigamma(phi))

  # Pre-calculate weighted derivatives
  W_s <- diag(W) * s
  W_rR <- diag(W) * rR
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  K_a_a <- phi * crossprod(s, W_s)
  K_t_a <- phi * crossprod(rR, W_s)
  K_t_t <- phi * crossprod(rR, W_rR)
  K_a_phi <- crossprod(s, mT_vc)
  K_t_phi <- crossprod(rR, mT_vc)
  K_phi_phi <- sum(diag(D))

  # Assemble the symmetric matrix
  K_a_t <- t(K_t_a)
  K_phi_a <- t(K_a_phi)
  K_phi_t <- t(K_t_phi)

  # Construct the Fisher information matrix
  fisher_info_mat <- rbind(
    cbind(K_a_a, K_a_t, K_a_phi),
    cbind(K_t_a, K_t_t, K_t_phi),
    cbind(K_phi_a, K_phi_t, K_phi_phi)
  )

  # Name the rows and columns of the matrix
  names_theta <- paste0("theta", ma)
  param_names <- c("alpha", names_theta, "phi")
  colnames(fisher_info_mat) <- param_names
  rownames(fisher_info_mat) <- param_names

  # Prepare fitted values as a time series object
  fitted_values <- ts(c(rep(NA, m), muhat),
                      start = start(y),
                      frequency = frequency(y))

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

#' @title Penalized Fisher Information Matrix for a BMA Model
#' @description This function computes the L2-penalized Fisher Information
#'   Matrix (FIM) for a Beta Moving Average (BMA) model.
#'
#' @param y A time series of data in the (0,1) interval.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_ma_ridge <- function(y,
                                ma,
                                link,
                                alpha,
                                theta,
                                phi,
                                penalty) {
  # Compute the unpenalized Fisher Information Matrix.
  fim_list <- inf_matrix_ma(
    y = y,
    ma = ma,
    link = link,
    alpha = alpha,
    theta = theta,
    phi = phi
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # Determine the number of terms used in BMA maximization
  q_max <- max(ma, na.rm = TRUE)
  a_max <- q_max
  n <- length(y)

  # Construct the penalty matrix.
  q_len <- length(theta)

  # Penalty applies to alpha and theta. No penalty for phi.
  penalty_diag <- c(1, rep(1, q_len), 0)
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
