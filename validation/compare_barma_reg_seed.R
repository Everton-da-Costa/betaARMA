# ============================================================================
# DESCRIPTION:
# Function to compare "new" vs "old" BARMA functions for a
# single, specified seed. Useful for debugging specific failures.
# ============================================================================

# The new code is now validated to be numerically identical to the original
# implementation for $\beta$ARMA(p,q) models.

# Key Improvements:

# S3 Method Extraction: 
# All post-estimation logic (forecasts, residuals, summary
# statistics) has been extracted from the main barma function into separate S3
# methods (forecast.barma, residuals.barma, summary.barma). This makes the code
# cleaner and aligns with R package standards.

# Unified Model Interface: 
# The new barma function now serves as a single interface, handling $\beta$ARMA,
# $\beta$AR, and $\beta$MA models automatically based on the ar and ma 
# arguments. This replaces the large, complex if/else structure of the original
# barma_old.R.

# Computational Efficiency (Vectorization): 
# The loglik_barma and score_vector_barma functions now use crossprod() for 
# vector multiplication inside the loop. This is significantly more efficient 
# than the previous approach in barma_old.R.

# Single-Pass Loop: 
# The new score_vector_barma function calculates the eta value and its 
# recursive derivatives in a single for loop, whereas the old code required two
# separate passes (one for eta/error, one for derivatives).This vignette runs a
# single-seed comparison to prove these changes produce identical numerical
# results.
# ============================================================================

rm(list = ls())

library(forecast)
library(dplyr)
library(here)
library(stats)

# ============================================================================
# Helper function
# ============================================================================
reorder_values <- function(fit) {
  # 1. Setup indices
  n_total <- length(fit$coef)
  n_ar    <- length(fit$varphi)
  n_ma    <- length(fit$theta)
  
  # The fixed ARMA part (AR + MA + Intercept)
  indx_arma <- seq_len(n_ar + n_ma + 1) 
  
  # Locate phi, after ARMA
  indx_phi <- n_ar + n_ma + 2
  
  # Locate Regressors: Everything else
  all_indices <- seq_len(n_total)
  indx_reg <- setdiff(all_indices, c(indx_arma, indx_phi))
  
  # 2. Define the target order: [ARMA] -> [Regressors] -> [phi]
  new_order <- c(indx_arma, indx_reg, indx_phi)
  
  # 3. Reorder
  # Using drop=FALSE ensures we keep matrix dimensions even if it's 1x1
  coef_ord   <- fit$coef[new_order]
  vcov_ord   <- fit$vcov[new_order, new_order, drop = FALSE]
  
  return(list(
    coef_old_ord = coef_ord,
    fit_vcov_ord = vcov_ord
  ))
}



# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# After fitting barma_old
# fit_old <- barma_old(y = y, ar = ar, ma = ma, X = X, X_hat = X_hat, ...)

# Reorder
# reorder_list <- reorder_values(fit_old)

# Extract reordered values
# coef_old <- reorder_list$coef_old_ord
# fisher_info_mat_old <- reorder_list$fisher_info_mat_old_ord
# fit_vcov <- reorder_list$fit_vcov_ord

# Verify parameter counts
# cat("n_ar:", reorder_list$n_ar, "\n")
# cat("n_ma:", reorder_list$n_ma, "\n")
# cat("n_beta:", reorder_list$n_beta, "\n")


# ============================================================================
# MAIN FUNCTION: Compare Single Seed
# ============================================================================
compare_barma_seed <- function(
    ar = NULL,
    ma = NULL,
    seed = 1,
    n = 250,
    burn = 1000,
    link = "logit",
    varphi_true = 0.4,
    theta_true = NA,
    alpha_true = 0,
    phi_true = 20,
    X = NA,
    X_hat = NA) {
  
  # ------------------------------------------------------------------------- #
  # --- 1. Load All Functions ---
  # ------------------------------------------------------------------------- #
  
  cat("============================================================\n")
  cat("BARMA Single Seed Comparison (Seed:", seed, ")\n")
  cat("============================================================\n")
  cat("Loading functions...\n")
  
  # --- Load OLD functions ---
  source(here("R", "barma_old.R"))
  
  # --- Load NEW functions ---
  source(here("R", "barma.R"))
  barma_new <- barma
  
  source(here("R", "loglik_barma.R"))
  source(here("R", "score_vector_barma.R"))
  source(here("R", "fim_barma.R"))
  
  source(here("R", "residuals.barma.R"))
  source(here("R", "fitted.barma.R"))
  source(here("R", "coef.barma.R"))
  source(here("R", "print.barma.R"))
  source(here("R", "summary.barma.R"))
  source(here("R", "print.summary.barma.R"))
  source(here("R", "forecast.barma.R"))
  
  # --- Load shared helpers ---
  source(here("R", "start_values.R"))
  source(here("R", "simu_barma.R"))
  source(here("R", "make_link_structure.R"))
  
  # ------------------------------------------------------------------------- #
  # --- 2. Generate Data for the Specific Seed ---
  # ------------------------------------------------------------------------- #
  
  cat("\nGenerating data for seed:", seed, "\n")
  
  nburn <- n + burn
  set.seed(seed)
  y_burn <- simu_barma(
    n = nburn,
    varphi = varphi_true,
    theta = theta_true,
    alpha = alpha_true,
    phi = phi_true
  )
  y <- ts(y_burn[(burn + 1):nburn], frequency = 12)
  
  # ------------------------------------------------------------------------- #
  # --- 3. Fit OLD Model ---
  # ------------------------------------------------------------------------- #
  
  cat("Fitting barma_old()...\n")
  time_old_run <- system.time({
    fit_old <- barma_old(
      y = y,
      ar = ar,
      ma = ma,
      link = link,
      h1 = 6,
      X = X,
      X_hat = X_hat
    )
  })
  cat("  Time (old):", round(time_old_run["elapsed"], 4), "s\n")
  
  # ------------------------------------------------------------------------- #
  # --- 4. Fit NEW Model ---
  # ------------------------------------------------------------------------- #
  if (any(is.na(X))) X = NULL
  
  cat("Fitting barma() (new)...\n")
  time_new_run <- system.time({
    fit_new <- barma_new(
      y = y,
      ar = ar,
      ma = ma,
      link = link,
      xreg = X
    )
  })
  cat("  Time (new):", round(time_new_run["elapsed"], 4), "s\n")
  
  # ------------------------------------------------------------------------- #
  # --- 5. Reorder values ---
  # ------------------------------------------------------------------------- #
  
  if (!is.null(fit_old$beta)) {
    reorder_values_list <- reorder_values(fit_old)
    
    coef_old = reorder_values_list$coef_old_ord
    vcov_old = reorder_values_list$fit_vcov_ord
    
  } else {
    
    coef_old = fit_old$coef
    vcov_old = fit_old$vcov
    
  }
  
  # ------------------------------------------------------------------------- #
  # Forecast
  # ------------------------------------------------------------------------- #
  
  if (!is.null(fit_old$beta)) {
    
    forecast_new <- forecast.barma(object = fit_new, xreg = X_hat,  h = 6)
    
  } else {
    
    forecast_new <- forecast.barma(object = fit_new,  h = 6)
    
  }
  
  
  
  # ------------------------------------------------------------------------- #
  # --- 6. Run All Comparisons ---
  # ------------------------------------------------------------------------- #
  cat("--------------------------------------------", "\n")
  print(
    rbind(
      barma_new = c(fit_new$conv, fit_new$coef),
      barma_old = c(fit_old$conv, coef_old),
      difference = c(NA, fit_new$coef - coef_old)
    )
  )
  cat("--------------------------------------------", "\n")
  
  cat("\n--- Comparison Results ---\n")
  
  # --- S3 Method outputs from NEW function ---
  start_values_new <- fit_new$start_values
  sum_new <- summary(fit_new)
  res_new <- residuals(fit_new)
  fit_new_fitted <- fitted(fit_new)
  coef_new <- coef(fit_new)
  
  
  # --- List outputs from OLD function ---
  start_values_old <- fit_old$start_values
  res_old <- fit_old$resid2
  fit_old_fitted <- fit_old$fitted
  # coef_old <- fit_old$coef
  # vcov_old <- fit_old$vcov
  sum_old_table <- fit_old$model
  loglik_old <- fit_old$loglik
  aic_old <- fit_old$aic
  bic_old <- fit_old$bic
  hq_old <- fit_old$hq
  
  forecast_old <- fit_old$forecast
  
  # --- Index for effective observations ---
  max_lag <- fit_new$max_lag
  idx_effective <- (max_lag + 1):length(y)
  
  # --- Helper function for printing results ---
  check_equality <- function(name, old_val, new_val, tol = 1e-8) {
    result <- all.equal(old_val, new_val, tolerance = tol)
    if (isTRUE(result)) {
      cat(sprintf("  ✓ %-12s: Identical (tol=%.0e)\n", name, tol))
    } else {
      cat(sprintf("  ✗ %-12s: DIFFERENT\n", name))
      print(result) # Print the details of the difference
    }
  }
  # --- C0. Start values ---
  cat("\n1. Comparing Start values:\n")
  check_equality("Start values", 
                 as.numeric(start_values_old), 
                 as.numeric(start_values_new))
  
  # --- C1. Coefficients ---
  cat("\n1. Comparing Coefficients:\n")
  check_equality("Coefficients", coef_old, coef_new)
  
  # --- C2. Log-Likelihood & Criteria ---
  cat("\n2. Comparing Log-Likelihood and Information Criteria:\n")
  check_equality("Log-Likelihood", loglik_old, fit_new$loglik)
  check_equality("AIC", aic_old, sum_new$aic)
  check_equality("BIC", bic_old, sum_new$bic)
  check_equality("HQ", hq_old, sum_new$hq)
  
  # --- C3. Model Outputs ---
  cat("\n3. Comparing Model Outputs:\n")
  check_equality("VCov Matrix", 
                 as.numeric(vcov_old), 
                 as.numeric(sum_new$vcov))
  
  # Compare effective values (stripping NAs)
  check_equality("Fitted Values",
                 fit_old_fitted[idx_effective],
                 as.numeric(fit_new_fitted[idx_effective])
  )
  check_equality("Residuals",
                 res_old,
                 as.numeric(res_new[idx_effective])
  )
  
  check_equality("Forecast",
                 as.numeric(forecast_old),
                 as.numeric(forecast_new)
  )
  
  cat("\n============================================================\n")
  cat("Comparison complete.\n")
  cat("============================================================\n\n")
  
  # Return both objects invisibly for inspection
  return(invisible(list(fit_new = fit_new, fit_old = fit_old)))
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# --- Run the comparison for a specific seed ---
# 
seed = 20
ar = 1
ma = 1:4
n = 100
burn = 1000
link = "logit"
varphi_true = 0.4
theta_true = NA
alpha_true = 0
phi_true = 20
X = NA
X_hat = NA

compare_barma_seed(
  ar = ar,
  ma = ma,
  seed = seed,
  n = n,
  burn = burn,
  link = link,
  varphi_true = varphi_true,
  theta_true = theta_true,
  alpha_true = alpha_true,
  phi_true = phi_true,
  X = X,
  X_hat = X_hat
)


vec_train <- 1:n
vec_test <- (max(vec_train) + 1):(max(vec_train) + 6)

X = cbind(
  hs =  sin(2 * pi * vec_train / 12),
  hc =  cos(2 * pi * vec_train / 12)
)

X_hat = cbind(
  hs_test =  sin(2 * pi * vec_test / 12),
  hc_test =  cos(2 * pi * vec_test / 12)
)

compare_barma_seed(
  ar = ar,
  ma = ma,
  seed = seed,
  n = n,
  burn = burn,
  link = link,
  varphi_true = varphi_true,
  theta_true = theta_true,
  alpha_true = alpha_true,
  phi_true = phi_true,
  X = X,
  X_hat = X_hat
)
