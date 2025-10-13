# =========================================================================== #
#           FUNCTIONS FOR THE BETA AUTOREGRESSIVE MOVING AVERAGE              #
#                           (BARMA) MODEL                                     #
# =========================================================================== #
# This file contains the core functions for the BARMA model:
#   - loglik_arma: Unpenalized log-likelihood
#   - loglik_arma_ridge: L2-penalized log-likelihood
#   - score_vector_arma: Unpenalized score vector (gradient)
#   - score_vector_arma_ridge: L2-penalized score vector
#   - inf_matrix_arma: Unpenalized Fisher Information Matrix
#   - inf_matrix_arma_ridge: L2-penalized Fisher Information Matrix
# --------------------------------------------------------------------------- #
#'
#' @param y Data, a time series of numbers in (0,1).
#' @param ar A vector specifying the autoregressive (AR) component.

#' @title Log-Likelihood for a BARMA Model
#' @description This function computes the log-likelihood of the BARMA model.
#' @param ma A vector specifying the moving average (MA) component.
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter of the BARMA model.
#' @param link The link function ("logit", "probit", "loglog", "cloglog").
#'
#' @importFrom stats dbeta
#'
#' @return The likelihood of the BARMA estimators.
#' @keywords internal
loglik_arma <- function(y, ar, ma,
                        alpha = 0, varphi = 0, theta = 0,
                        phi = 0, link) {

  # Link functions
  # ----------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)

  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta

  ynew <- linkfun(y)

  # ----------------------------------------------------------------------- #
  p <-  max(ar)
  q <-  max(ma)

  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # --------------------------------------------------------------------- #
  error <- rep(0, n)
  eta   <- rep(NA, n)
  mu    <- rep(NA, n)

  for (t in (m + 1):n) {

    eta[t] <- (alpha +
                 crossprod(varphi, ynew[t - ar]) +
                 crossprod(theta, error[t - ma]))

    error[t]  <- ynew[t] - eta[t]
  }

  mu1  <- linkinv(eta = eta[(m + 1):n])
  y1   <- y[(m + 1):n]

  ll_terms_arma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  final <- sum(ll_terms_arma)

  return(final)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Log-Likelihood for a BARMA Model
#' @description Computes the L2-penalized (Ridge) log-likelihood for a BARMA
#'   model.
#'
#' @details The penalty is applied to the intercept (`alpha`), AR (`varphi`),
#'   and MA (`theta`) coefficients. The precision parameter (`phi`) is not
#'   penalized.
#'
#' @param y A numeric vector or ts object of data in (0,1).
#' @param ar A vector of positive integer AR lags.
#' @param ma A vector of positive integer MA lags.
#' @param alpha The numeric intercept term.
#' @param varphi Vector of AR parameters, matching the length of `ar`.
#' @param theta Vector of MA parameters, matching the length of `ma`.
#' @param phi The positive precision parameter.
#' @param link The link function (e.g., "logit", "probit").
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A single numeric value: the penalized log-likelihood.
#' @keywords internal
loglik_arma_ridge <- function(y,
                              ar,
                              ma,
                              alpha,
                              varphi,
                              theta,
                              phi,
                              link,
                              penalty) {
  # Compute the standard (unpenalized) log-likelihood.
  # --------------------------------------------------------------------- #
  loglik <- loglik_arma(
    y,
    ar = ar,
    ma = ma,
    alpha = alpha,
    varphi = varphi,
    theta = theta,
    phi = phi,
    link = link
  )

  # To determine the number of terms used in BARMA maximization
  # --------------------------------------------------------------------- #
  p_max <- max(ar)
  q_max <- max(ma)
  a_max <- max(p_max, q_max)

  # sample size
  n <- length(y)

  # Calculate the penalty term.
  # --------------------------------------------------------------------- #
  # The L2-norm of penalized coefficients
  l2_norm_sq <- alpha^2 + sum(varphi^2) + sum(theta^2)
  scaled_penalty <- (n - a_max) * l2_norm_sq

  # Return the penalized log-likelihood
  ell2_penalty <- loglik - penalty * scaled_penalty

  return(ell2_penalty)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Score Vector for a BARMA Model
#' @description This function computes the score vector for a Beta
#'   Autoregressive Moving Average (BARMA) model. The score vector is the
#'   gradient of the log-likelihood function with respect to the model
#'   parameters.
#'
#' @param y A time series of data, with values in the interval (0, 1).
#' @param ar A numeric vector of positive integers specifying the AR lags.
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param alpha The numeric intercept term in the linear predictor.
#' @param varphi A numeric vector of AR parameters.
#' @param theta A numeric vector of MA parameters.
#' @param phi The numeric positive precision parameter.
#' @param link A character string for the link function, e.g., "logit".
#'
#' @return A numeric vector containing the score values, ordered as
#'   (alpha, varphi, theta, phi).
#' @keywords internal
score_vector_arma <- function(y, ar, ma,
                              alpha = 0, varphi = 0, theta = 0,
                              phi = 0, link) {

  # Link functions setup
  # ----------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  # ------------------------------------------------------------------------- #
  p <-  max(ar)
  q <-  max(ma)
  p1 <- length(ar)
  q1 <- length(ma)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # Initialize vectors for the recursive calculations
  # ------------------------------------------------------------------------- #
  error <- rep(0, n)
  eta <- rep(NA, n)

  # Recursively compute the linear predictor (eta) and errors
  for (t in (m + 1):n) {
    eta[t] <- alpha +
      crossprod(varphi, ynew[t - ar]) +
      crossprod(theta, error[t - ma])
    error[t]  <- ynew[t] - eta[t]
  }

  # Subset all series to the effective sample size (n - m)
  eta1 <- eta[(m + 1):n]
  mu1  <- linkinv(eta = eta1)
  y1   <- y[(m + 1):n]

  # Derivatives of the log-likelihood function
  # ------------------------------------------------------------------------- #
  # Design matrices for AR (P) and MA (R) components
  P <- matrix(nrow = n - m, ncol = p1)
  for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

  R <- matrix(nrow = n - m, ncol = q1)
  for (t in 1:(n - m)) R[t, ] <- error[t + m - ma]

  # Recursively compute the derivatives of eta w.r.t. each parameter
  deta_dalpha  <- rep(0, n)
  deta_dvarphi <- matrix(0, nrow = n, ncol = p1)
  deta_dtheta  <- matrix(0, nrow = n, ncol = q1)

  for (t in (m + 1):n) {
    deta_dalpha[t]    <- 1 - crossprod(theta, deta_dalpha[t - ma])
    deta_dvarphi[t, ] <- P[t - m, ] - crossprod(theta, deta_dvarphi[t - ma, ])
    deta_dtheta[t, ]  <- R[t - m, ] - crossprod(theta, deta_dtheta[t - ma, ])
  }

  # Subset derivatives to the effective sample size
  s  <- deta_dalpha[(m + 1):n]
  rP <- deta_dvarphi[(m + 1):n, ]
  rR <- deta_dtheta[(m + 1):n, ]

  # Compute score vector components
  # ------------------------------------------------------------------------- #
  mu_eta <- mu.eta(eta = eta1)
  ystar  <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

  ystar_mustar <- ystar - mustar
  mT_ystar_mustar <- mu_eta * ystar_mustar

  # Score component for each parameter
  U_alpha  <- phi * crossprod(s, mT_ystar_mustar)
  U_varphi <- phi * crossprod(rP, mT_ystar_mustar)
  U_theta  <- phi * crossprod(rR, mT_ystar_mustar)
  U_phi   <- sum(
    mu1 * ystar_mustar + log(1 - y1) -
      digamma((1 - mu1) * phi) +
      digamma(phi)
  )

  # Combine into a single score vector
  escore_vec <- c(U_alpha, U_varphi, U_theta, U_phi)

  return(escore_vec)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Penalized Score Vector for a BARMA Model
#' @description This function computes the L2-penalized score vector for a
#'   BARMA model.
#'
#' @param y A time series object of data in (0,1).
#' @param ar A numeric vector of positive integers specifying the AR lags.
#' @param ma A numeric vector of positive integers specifying the MA lags.
#' @param alpha The numeric intercept term.
#' @param varphi A numeric vector of AR parameters.
#' @param theta A numeric vector of MA parameters.
#' @param phi The numeric positive precision parameter (not penalized).
#' @param link A character string specifying the link function.
#' @param penalty The numeric regularization parameter (lambda).
#'
#' @return A numeric vector for the L2-penalized score. The order is:
#'   (alpha, varphi, theta, phi).
#' @keywords internal
score_vector_arma_ridge <- function(y,
                                    ar, ma,
                                    alpha, varphi, theta, phi,
                                    link, penalty) {

  # Compute the standard (unpenalized) score vector.
  unpenalized_score <- score_vector_arma(
    y = y,
    ar = ar,
    ma = ma,
    alpha = alpha,
    varphi = varphi,
    theta = theta,
    phi = phi,
    link = link
  )

  # Determine the maximum lag order.
  # --------------------------------------------------------------------- #
  p_max <- max(ar, na.rm = TRUE)
  q_max <- max(ma, na.rm = TRUE)
  a_max <- max(p_max, q_max)

  # Sample size
  n <- length(y)

  # Compute the gradient of the penalty term.
  # --------------------------------------------------------------------- #
  # Derivative of penalty w.r.t (alpha, varphi, theta, phi)
  l2_norm_sq_gradient <- c(2 * alpha, 2 * varphi, 2 * theta, 0)

  # Scale the penalty gradient.
  penalty_gradient <- penalty * (n - a_max) * l2_norm_sq_gradient

  # Compute the final penalized score.
  penalized_score <- unpenalized_score - penalty_gradient

  return(penalized_score)
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#' @title Fisher Information Matrix for a BARMA Model
#' @description This function computes the Fisher information matrix of the
#'   Beta Autoregressive Moving Average (BARMA) model. It also computes
#'   auxiliary values like fitted values and residuals.
#'
#' @param y Data, a time series of numbers in the open interval (0,1).
#' @param ar A numeric vector specifying the autoregressive (AR) lags.
#' @param ma A numeric vector specifying the moving average (MA) lags.
#' @param alpha The intercept term.
#' @param varphi A numeric vector of autoregressive (AR) parameters.
#' @param theta A numeric vector of moving average (MA) parameters.
#' @param phi The precision parameter of the BARMA model.
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
inf_matrix_arma <- function(y,
                            ar,
                            ma,
                            alpha = 0,
                            varphi = 0,
                            theta = 0,
                            phi = 0,
                            link) {

  # Link functions
  # ----------------------------------------------------------------------- #
  link_structure <- make_link_structure(link)
  linkfun <- link_structure$linkfun
  linkinv <- link_structure$linkinv
  mu.eta  <- link_structure$mu.eta

  ynew <- linkfun(y)

  # Model dimensions
  # ----------------------------------------------------------------------- #
  p <- max(ar)
  q <- max(ma)
  p1 <- length(ar)
  q1 <- length(ma)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE)

  # Initialize vectors for fitted values and errors
  # ------------------------------------------------------------------------- #
  errorhat <- rep(0, n)
  etahat <- rep(NA, n)

  # Compute fitted values and errors (errorhat) recursively
  # ------------------------------------------------------------------------- #
  for (t in (m + 1):n) {
    etahat[t]    <- (alpha +
                       crossprod(varphi, ynew[t - ar]) +
                       crossprod(theta, errorhat[t - ma]))

    errorhat[t] <- ynew[t] - etahat[t]
  }

  # Truncate etahat and muhat to observed range
  etahat_obs <- etahat[(m + 1):n]
  muhat <- linkinv(etahat_obs)

  # Prepare lagged AR and MA components for derivative calculations
  # ------------------------------------------------------------------------- #
  P <- matrix(nrow = n - m, ncol = p1)
  for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

  R <- matrix(nrow = n - m, ncol = q1)
  for (t in 1:(n - m)) R[t, ] <- errorhat[t + m - ma]

  # Recursive calculation of eta derivatives w.r.t. parameters
  # ------------------------------------------------------------------------- #
  deta_dalpha  <- rep(0, n)
  deta_dvarphi <- matrix(0, nrow = n, ncol = p1)
  deta_dtheta  <- matrix(0, nrow = n, ncol = q1)

  for (t in (m + 1):n) {
    deta_dalpha[t]    <- 1 - crossprod(theta, deta_dalpha[t - ma])
    deta_dvarphi[t, ] <- P[t - m, ] - crossprod(theta, deta_dvarphi[t - ma, ])
    deta_dtheta[t, ]  <- R[t - m, ] - crossprod(theta, deta_dtheta[t - ma, ])
  }

  # Truncate derivatives to the observed range
  s  <- deta_dalpha[(m + 1):n]
  rP <- deta_dvarphi[(m + 1):n, ]
  rR <- deta_dtheta[(m + 1):n, ]

  # Components for Fisher Information Matrix calculation
  # ------------------------------------------------------------------------- #
  one_minus_muhat <- 1 - muhat
  psi1 <- trigamma(muhat * phi)
  psi2 <- trigamma(one_minus_muhat * phi)

  mu_eta_deriv <- mu.eta(eta = etahat_obs)
  mT <- diag(mu_eta_deriv)

  W_diag_elements <- phi * (psi1 + psi2)
  W <- diag(W_diag_elements) %*% mT^2

  vc <- phi * (psi1 * muhat - psi2 * one_minus_muhat)
  D <- diag(psi1 * muhat^2 + psi2 * one_minus_muhat^2 - trigamma(phi))

  W_s <- diag(W) *  s
  W_rP <- diag(W) * rP
  W_rR <- diag(W) * rR
  mT_vc <- mu_eta_deriv * vc

  # Compute blocks of the Fisher Information Matrix
  # ------------------------------------------------------------------------- #
  K_a_a <- phi * crossprod(s, W_s)
  K_p_p <- phi * crossprod(rP, W_rP)
  K_t_t <- phi * crossprod(rR, W_rR)

  K_p_a <- phi * crossprod(rP, W_s)
  K_t_a <- phi * crossprod(rR, W_s)
  K_p_t <- phi * crossprod(rP, W_rR)

  K_a_phi <- crossprod(s, mT_vc)
  K_p_phi <- crossprod(rP, mT_vc)
  K_t_phi <- crossprod(rR, mT_vc)

  K_phi_phi <- sum(diag(D))

  # Assemble the symmetric matrix
  K_a_p <- t(K_p_a)
  K_a_t <- t(K_t_a)
  K_t_p <- t(K_p_t)
  K_phi_a <- K_a_phi
  K_phi_p <- t(K_p_phi)
  K_phi_t <- t(K_t_phi)

  # Construct the Fisher information matrix
  # Order: alpha, varphi (AR), theta (MA), phi
  # ------------------------------------------------------------------------- #
  fisher_info_mat <- rbind(
    cbind(K_a_a, K_a_p, K_a_t, K_a_phi),
    cbind(K_p_a, K_p_p, K_p_t, K_p_phi),
    cbind(K_t_a, K_t_p, K_t_t, K_t_phi),
    cbind(K_phi_a, K_phi_p, K_phi_t, K_phi_phi)
  )

  # Name the rows and columns of the matrix for clarity
  # ------------------------------------------------------------------------- #
  names_varphi <- paste0("varphi", ar)
  names_theta  <- paste0("theta", ma)

  names_fisher_info_mat <- c("alpha", names_varphi, names_theta, "phi")
  colnames(fisher_info_mat) <- names_fisher_info_mat
  rownames(fisher_info_mat) <- names_fisher_info_mat

  # Prepare fitted values as a time series object
  # ------------------------------------------------------------------------- #
  fitted_values <- ts(c(rep(NA, m), muhat),
                      start = start(y),
                      frequency = frequency(y))

  # Prepare output list
  # ------------------------------------------------------------------------- #
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

#' @title Penalized Fisher Information Matrix for a BARMA Model
#' @description This function computes the L2-penalized Fisher Information
#'   Matrix (FIM) for a BARMA model.
#'
#' @param y A time series of data in the (0,1) interval.
#' @param ar A vector specifying the autoregressive (AR) lags.
#' @param ma A vector specifying the moving average (MA) lags.
#' @param link The link function ("logit", "probit", "cloglog").
#' @param alpha The intercept term.
#' @param varphi A vector of autoregressive (AR) parameters.
#' @param theta A vector of moving average (MA) parameters.
#' @param phi The precision parameter (not penalized).
#' @param penalty The non-negative L2 penalty parameter (lambda).
#'
#' @return A list containing the penalized FIM and other model outputs.
#' @keywords internal
inf_matrix_arma_ridge <- function(y,
                                  ar,
                                  ma,
                                  link,
                                  alpha,
                                  varphi,
                                  theta,
                                  phi,
                                  penalty) {
  # Compute the unpenalized Fisher Information Matrix.
  fim_list <- inf_matrix_arma(
    y = y,
    ar = ar,
    ma = ma,
    link = link,
    alpha = alpha,
    varphi = varphi,
    theta = theta,
    phi = phi
  )
  unpenalized_fim <- fim_list$fisher_info_mat

  # To determine the number of terms used in BARMA maximization
  # --------------------------------------------------------------------- #
  p_max <- max(ar, na.rm = TRUE)
  q_max <- max(ma, na.rm = TRUE)
  a_max <- max(p_max, q_max)

  # Sample size
  n <- length(y)

  # Construct the penalty matrix.
  # --------------------------------------------------------------------- #
  p_len <- length(varphi)
  q_len <- length(theta)

  # Penalty applies to alpha, varphi, and theta. No penalty for phi.
  penalty_diag <- c(1, rep(1, p_len), rep(1, q_len), 0)
  penalty_matrix <- 2 * diag(penalty_diag)

  # Compute the penalized Fisher Information Matrix.
  penalized_fim <- unpenalized_fim + penalty * (n - a_max) * penalty_matrix

  # Prepare the output list.
  # ------------------------------------------------------------------------- #
  output_list <- list(
    fisher_info_mat_penalized = penalized_fim,
    fitted = fim_list$fitted,
    muhat = fim_list$muhat,
    etahat = fim_list$etahat,
    errorhat = fim_list$errorhat
  )

  return(output_list)
}
