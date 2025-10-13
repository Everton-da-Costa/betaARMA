# =========================================================================== #
#           FUNCTIONS FOR THE BETA AUTOREGRESSIVE (BAR) MODEL                 #
# =========================================================================== #
# This file contains the core functions for the BAR model:
#   - loglik_ar: Unpenalized log-likelihood
#   - loglik_ar_ridge: L2-penalized log-likelihood
#   - score_vector_ar: Unpenalized score vector (gradient)
#   - score_vector_ar_ridge: L2-penalized score vector
#   - inf_matrix_ar: Unpenalized Fisher Information Matrix
#   - inf_matrix_ar_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #

#' @title Log-Likelihood for a BAR Model
#' @description This function computes the log-likelihood of the Beta
#'   Autoregressive (BAR) model.
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) component.
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function to be used ("logit", "probit", "loglog",
#'   "cloglog").
#'
#' @importFrom stats dbeta
#'
#' @return A numeric value representing the total conditional log-likelihood
#'   of the BAR model for the given data and parameters.
#' @keywords internal
loglik_ar <- function(y,
                      ar,
                      alpha = 0,
                      varphi = 0,
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
  mu <- rep(NA, n)

  # Compute eta and mu for t > m
  for (t in (m + 1):n) {
    eta[t] <- alpha + varphi %*% ynew[t - ar]
    mu[t] <- linkinv(eta[t])
  }

  # Subset for likelihood computation
  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # Compute the log-likelihood using the Beta distribution density
  ll_terms_ar <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  # Sum the individual log-likelihood terms to get the total
  sum_ll <- sum(ll_terms_ar)

  return(sum_ll)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BAR Model
#' @description Computes the L2-penalized (Ridge) log-likelihood for a BAR
#'   model.
#'
#' @details The penalty is applied to the intercept (`alpha`) and AR (`varphi`)
#'   coefficients. The precision parameter (`phi`) is not penalized.
#'
#' @param y A numeric vector or ts object of data in (0,1).
#' @param ar A vector of positive integer AR lags.
#' @param alpha The numeric intercept term.
#' @param varphi Vector of AR parameters, matching the length of `ar`.
#' @param phi The positive precision parameter.
#' @param link The link function (e.g., "logit", "probit").
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A single numeric value: the penalized log-likelihood.
#' @keywords internal
loglik_ar_ridge <- function(y,
                            ar,
                            alpha,
                            varphi,
                            phi,
                            link,
                            penalty) {
  # Compute the standard (unpenalized) log-likelihood.
  loglik <- loglik_ar(
    y,
    ar = ar,
    alpha = alpha,
    varphi = varphi,
    phi = phi,
    link = link
  )

  # Determine the number of terms used in BAR maximization
  p_max <- max(ar)
  a_max <- max(p_max)
  n <- length(y)

  # Calculate the penalty term.
  l2_norm_sq <- alpha^2 + sum(varphi^2)
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Return the penalized log-likelihood
  ell2_penalty <- loglik - penalty * scaled_penalty

  return(ell2_penalty)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BAR Model
#' @description This function computes the score vector of the Beta
#'   Autoregressive (BAR) model.
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) component.
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param phi The precision parameter of the Beta distribution.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @return The score vector of the BAR estimators.
#' @keywords internal
score_vector_ar <- function(y,
                            ar,
                            alpha = 0,
                            varphi = 0,
                            phi = 0,
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
  n <- length(y)
  m <- max(p, na.rm = TRUE)

  # Initialize vectors
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  # Compute eta and mu for t > m
  for (t in (m + 1):n) {
    eta[t] <- alpha + varphi %*% ynew[t - ar]
    mu[t] <- linkinv(eta[t])
  }

  # Subset for score computation
  eta1 <- eta[(m + 1):n]
  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # Score components
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

  # Derivative matrix and vector
  mu_eta <- mu.eta(eta = eta1)
  mT <- diag(mu_eta)

  # Design matrix for AR terms
  P <- matrix(nrow = n - m, ncol = p1)
  for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

  # Intermediate score calculations
  ystar_mustar <- ystar - mustar
  mT_ystar_mustar <- mu_eta * ystar_mustar

  # Score vector components
  U_alpha <- phi * sum(mu_eta * ystar_mustar)
  U_varphi <- phi * crossprod(P, mT_ystar_mustar)
  U_phi <- sum(mu1 * ystar_mustar +
                 log(1 - y1) -
                 digamma((1 - mu1) * phi) +
                 digamma(phi))

  # Return complete score vector
  score_vec <- c(U_alpha, U_varphi, U_phi)

  return(score_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BAR Model
#' @description This function computes the L2-penalized (Ridge) score vector
#'   for a Beta Autoregressive (BAR) model.
#'
#' @param y A time series object of data in (0,1).
#' @param ar A numeric vector of positive integers specifying the AR lags.
#' @param alpha The numeric intercept term in the linear predictor.
#' @param varphi A numeric vector of AR parameters.
#' @param phi The numeric positive precision parameter (not penalized).
#' @param link A character string specifying the link function.
#' @param penalty The numeric regularization parameter (lambda).
#'
#' @return A numeric vector representing the L2-penalized score.
#' @keywords internal
score_vector_ar_ridge <- function(y,
                                  ar,
                                  alpha,
                                  varphi,
                                  phi,
                                  link,
                                  penalty) {
  # Compute the standard (unpenalized) score vector.
  unpenalized_score <- score_vector_ar(
    y = y,
    ar = ar,
    alpha = alpha,
    varphi = varphi,
    phi = phi,
    link = link
  )

  # Determine the number of terms used in BAR maximization
  p_max <- max(ar, na.rm = TRUE)
  a_max <- p_max
  n <- length(y)

  # Compute the gradient of the penalty term.
  # Derivative of penalty w.r.t. (alpha, varphi, phi)
  l2_norm_sq_gradient <- c(2 * alpha, 2 * varphi, 0)

  # Scale the penalty gradient.
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient

  # Compute the final penalized score.
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BAR Model
#' @description This function computes the Fisher information matrix of the
#'   Beta Autoregressive (BAR) model.
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#' @param alpha The intercept term.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#' @param phi The precision parameter of the BAR model.
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
inf_matrix_ar <- function(y,
                          ar,
                          alpha = 0,
                          varphi = 0,
                          phi = 0,
                          link) {

  # Link functions
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  p <- max(ar)
  p1 <- length(ar)
  n <- length(y)
  m <- max(p, na.rm = TRUE)

  # Initialize vectors for fitted values and errors
  etahat <- rep(NA, n)
  errorhat <- rep(0, n)

  # Compute fitted values and errors (errorhat) recursively
  for (t in (m + 1):n) {
    etahat[t] <- alpha + varphi %*% ynew[t - ar]
    errorhat[t] <- ynew[t] - etahat[t]
  }

  # Truncate etahat and muhat to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # P: Matrix of lagged ynew values corresponding to AR lags
  P <- matrix(nrow = n - m, ncol = p1)
  for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

  # Components for Fisher Information Matrix calculation
  one_minus_muhat <- 1 - muhat
  psi1 <- trigamma(muhat * phi)
  psi2 <- trigamma(one_minus_muhat * phi)

  mu_eta_deriv <- mu.eta(eta = etahat_obs)
  mT <- diag(mu_eta_deriv)

  W_diag_elements <- phi * (psi1 + psi2)
  W <- diag(W_diag_elements) %*% mT^2

  vc <- phi * (psi1 * muhat - psi2 * one_minus_muhat)
  D <- diag(psi1 * (muhat^2) + psi2 * one_minus_muhat^2 - trigamma(phi))

  vI <- as.vector(rep(1, n - m))

  # Pre-calculate weighted derivatives
  W_vI <- diag(W) * vI
  W_P <- diag(W) * P
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  K_a_a <- as.matrix(phi * sum(diag(W)))
  K_p_a <- phi * crossprod(P, W_vI)
  K_p_p <- phi * crossprod(P, W_P)

  K_a_phi <- crossprod(vI, mT_vc)
  K_p_phi <- crossprod(P, mT_vc)
  K_phi_phi <- sum(diag(D))

  # Assemble the symmetric matrix
  K_phi_a <- K_a_phi
  K_a_p <- t(K_p_a)
  K_phi_p <- t(K_p_phi)

  # Construct the Fisher information matrix
  fisher_info_mat <- rbind(
    cbind(K_a_a, K_a_p, K_a_phi),
    cbind(K_p_a, K_p_p, K_p_phi),
    cbind(K_phi_a, K_phi_p, K_phi_phi)
  )

  # Name the rows and columns of the matrix for clarity
  names_varphi <- paste0("varphi", ar)
  names_fisher_info_mat <- c("alpha", names_varphi, "phi")
  colnames(fisher_info_mat) <- names_fisher_info_mat
  rownames(fisher_info_mat) <- names_fisher_info_mat

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

#' @title Penalized Fisher Information Matrix for a BAR Model
#' @description This function computes the L2-penalized Fisher Information
#'   Matrix (FIM) for a Beta Autoregressive (BAR) model.
#'
#' @param y A time series of data in the (0,1) interval.
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param phi The precision parameter (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_ar_ridge <- function(y,
                                ar,
                                link,
                                alpha,
                                varphi,
                                phi,
                                penalty) {
  # Compute the unpenalized Fisher Information Matrix.
  fim_list <- inf_matrix_ar(
    y = y,
    ar = ar,
    link = link,
    alpha = alpha,
    varphi = varphi,
    phi = phi
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # Determine the number of terms used in BAR maximization
  p_max <- max(ar, na.rm = TRUE)
  a_max <- p_max
  n <- length(y)

  # Construct the penalty matrix.
  p_len <- length(varphi)

  # Penalty applies to alpha and varphi. No penalty for phi.
  penalty_diag <- c(1, rep(1, p_len), 0)
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
