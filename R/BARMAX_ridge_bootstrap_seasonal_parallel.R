#' Parallelized Seasonal Block Bootstrap for BARMA(X) Models
#'
#' @description
#' Implements a robust, parallelized seasonal block bootstrap for BARMA(X) models.
#' This function is designed to provide stable parameter estimates by resampling
#' the data and refitting a ridge-penalized model, as described in
#' Cribari-Neto, Costa, and Fonseca (2025). It includes a retry mechanism to
#' handle model non-convergence.
#'
#' @param n_cores Number of CPU cores for parallel processing. Defaults to 1.
#' @param y The time series data, typically a 'ts' object.
#' @param seed Integer for the random number generator seed for reproducibility.
#' @param ar_vec Numeric vector for autoregressive (AR) lags.
#' @param ma_vec Numeric vector for moving average (MA) lags.
#' @param fit_ridge An object from the `barma` function for the original data.
#' @param penalty The penalty term (lambda) for the ridge regression.
#' @param regressors_list A list containing the time series (`y_sample_ts`) and
#'   all required regressor vectors (`hs`, `hc`, `hs_hat`, `hc_hat`).
#' @param block_length Integer specifying the block size for the bootstrap.
#' @param n_boot_rep The number of bootstrap replications. Defaults to 200.
#'
#' @return
#' A list containing the bootstrap results, including a matrix of coefficient
#' estimates (`mat_boot_estimates`) and their means
#' (`mean_bootstrap_estimates`), along with convergence statistics.
#'
#' @return
#' A list containing the bootstrap results, including a matrix of coefficient
#' estimates (`mat_boot_estimates`) and their means
#' (`mean_bootstrap_estimates`), along with convergence statistics.
#'
#' @references
#' Cribari-Neto, F., Costa, E., & Fonseca, R. V. (2025). Numerical stability
#' enhancements in beta autoregressive moving average model estimation.
#' *Manuscript submitted for publication*.
#'
#' Dudek, A. E., Leskow, J., Paparoditis, E., & Politis, D. N. (2014).
#' A GENERALIZED BLOCK BOOTSTRAP FOR SEASONAL TIME SERIES.
#' *Journal of Time Series Analysis*, 35(2), 89–114.
#'
#' @import doMC
#' @import foreach
#' @importFrom doRNG %dorng%
#' @export
BARMAX_ridge_bootstrap_seasonal_parallel <- function(n_cores = 1,
                                                     y = NULL,
                                                     seed = 3,
                                                     ar_vec = NA,
                                                     ma_vec = NA,
                                                     fit_ridge = NULL,
                                                     penalty = NULL,
                                                     regressors_list = NULL,
                                                     block_length = NULL,
                                                     n_boot_rep = 200) {
  # --------------------------------------------------------------------------
  # Input Validation
  # --------------------------------------------------------------------------
  stopifnot(
    "fit_ridge cannot be NULL" = !is.null(fit_ridge),
    "regressors_list cannot be NULL" = !is.null(regressors_list)
  )

  n_coef <- length(fit_ridge$coef)

  # --------------------------------------------------------------------------
  # Parallel settings
  # --------------------------------------------------------------------------
  registerDoMC(n_cores)

  comb <- function(...) {
    mapply("rbind", ..., SIMPLIFY = FALSE)
  }

  # --------------------------------------------------------------------------
  # Bootstrap with foreach loop
  # --------------------------------------------------------------------------
  res <- foreach(
    i = 1:n_boot_rep,
    .combine = "comb",
    .multicombine = TRUE,
    .options.RNG = seed
  ) %dorng% {
    flag <- TRUE
    count_fail_fit_boot <- c()
    count_conv_fail_boot <- c()
    count_conv_ok_boot <- c()

    while (flag) {
      block_bootstrap_seasonal <- bootstrap_seasonal(
        y = y,
        block_length = block_length,
        seasonal_period = 12
      )

      ids_y_bootstrap <- block_bootstrap_seasonal$ids_y_bootstrap

      # Fit BARMA ridge to the bootstrap sample with corrected list access
      fit_ridge_boot <- barma(
        y = ts(y[ids_y_bootstrap], frequency = 12),
        ar = ar_vec,
        ma = ma_vec,
        X = cbind(
          hs = regressors_list$hs,
          hc = regressors_list$hc
        ),
        X_hat = cbind(
          hs_hat = regressors_list$hs_hat,
          hc_hat = regressors_list$hc_hat
        ),
        penalty = penalty
      )

      # Check convergence status
      if (!is.list(fit_ridge_boot)) {
        flag <- TRUE
        count_fail_fit_boot <- rbind(count_fail_fit_boot, 1)
      } else if (fit_ridge_boot$conv != 0) {
        flag <- TRUE
        count_conv_fail_boot <- rbind(count_conv_fail_boot, 1)
      } else if (fit_ridge_boot$conv == 0) {
        flag <- FALSE
        count_conv_ok_boot <- rbind(count_conv_ok_boot, 1)
        mat_boot_estimates <- fit_ridge_boot$coef
      }
    } # end while

    # Return results for this replication
    list(
      mat_boot_estimates = mat_boot_estimates,
      count_fail_fit_boot = count_fail_fit_boot,
      count_conv_ok_boot = count_conv_ok_boot,
      count_conv_fail_boot = count_conv_fail_boot
    )
  }

  # --------------------------------------------------------------------------
  # Post-processing results
  # --------------------------------------------------------------------------
  mat_boot_estimates <- data.frame(res$mat_boot_estimates)
  mean_bootstrap_estimates <- as.numeric(colMeans(mat_boot_estimates))

  # Final output list
  output_list <- list(
    fit_ridge = fit_ridge,
    opt_ridge = list(conv = fit_ridge$conv, coef = fit_ridge$coef),
    mat_boot_estimates = mat_boot_estimates,
    mean_bootstrap_estimates = mean_bootstrap_estimates,
    fail_fit_boot = sum(unlist(res$count_fail_fit_boot)),
    conv_ok_boot = sum(unlist(res$count_conv_ok_boot)),
    conv_fail_boot = sum(unlist(res$count_conv_fail_boot))
  )

  return(output_list)
}
