#' Block Bootstrap for Seasonal Time Series (Dudek et al., 2013)
#'
#' Implements the generalized block bootstrap method for seasonal time series
#' data, as described in Dudek et al. (2013). This method is designed to
#' resample seasonal time series in a way that preserves its periodic structure.
#'
#' @param y A time series object containing the original time series data.
#' @param block_length An integer specifying the block size (`b`),
#'    which determines the length of contiguous observations sampled together
#'    during resampling.
#' @param seasonal_period An integer specifying the seasonal period length
#'    (`d`), such as 12 (for monthly data with yearly seasonality) or
#'    24, 36, etc.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{y_star_ts}{A time series object (`ts`) containing the bootstrapped
#'   time series data.}
#'   \item{ids_y_bootstrap}{A numeric vector containing the indices of the
#'   sampled blocks for each bootstrap iteration.}
#' }
#'
#' @details
#' This function implements the following steps (from Dudek et al., 2013):
#'
#' - **Step 1**: Define block size (`b`) and calculate the number of blocks
#'    (`l`). Identify the set of starting indices for the blocks.
#'
#' - **Step 2**: For each block starting index (`t`), calculate the set of
#'   possible indices (`S_\{t,n\}`) from which blocks can be sampled.
#'   This ensures that the seasonal structure is preserved.
#'
#' - **Step 3**: Randomly sample blocks from the time series based on the sets
#'   calculated in Step 2, and assemble the bootstrapped series by joining
#'   these blocks.
#'
#' The method ensures that the seasonal dependencies in the data are maintained
#' by carefully defining the set of candidate indices for resampling.
#'
#' @references
#' Dudek, A. E., Leskow, J., Paparoditis, E., & Politis, D. N. (2014).
#' A GENERALIZED BLOCK BOOTSTRAP FOR SEASONAL TIME SERIES.
#' _Journal of Time Series Analysis_, 35(2), 89–114. DOI: 10.1111/jtsa.12053
#'
#' @keywords internal
bootstrap_seasonal <- function(y,
                               block_length,
                               seasonal_period = 12) {

  # Input Validation
  stopifnot(
    "y must be a numeric vector or ts object" =
      is.numeric(y) || is.ts(y),
    "y must not be empty" =
      length(y) > 0,
    "block_length must be a positive integer" =
      block_length > 0 && block_length == round(block_length),
    "seasonal_period must be a positive integer" =
      seasonal_period > 0 && seasonal_period == round(seasonal_period),
    "block_length cannot be greater than the length of the time series" =
      block_length <= length(y)
  )

  # Rename objects to align with the notation in Dudek et al. (2013)
  d <- seasonal_period # Seasonal period length (d), e.g., 12, 24, 36, ...
  b <- block_length # Block length (b)

  # ------------------------------------------------------------------------- #
  # Step 1: Define parameters
  # ------------------------------------------------------------------------- #
  n <- length(y) # Number of observations
  l <- floor(n / b) # Number of blocks
  ind_set_t_end <- l * b + 1 # Definition of the last index in the block set

  # Set of indices: t = 1, b + 1, 2b + 1, ..., lb + 1
  ind_set_t <- seq(1, ind_set_t_end, b)

  # ------------------------------------------------------------------------- #
  # Step 2: Calculate the set S_{t,n}
  # ------------------------------------------------------------------------- #
  # Function to calculate R1,n and R2,n
  calc_R <- function(t, n, b, d) {
    R1 <- floor((t - 1) / d)
    R2 <- floor((n - b - t) / d)
    return(list(R1 = R1, R2 = R2))
  }

  # Function to calculate the set S_{t,n}
  calc_S <- function(t, n, b, d) {
    R <- calc_R(t, n, b, d)
    R1 <- R$R1
    R2 <- R$R2
    S <- seq(t - d * R1, t + d * R2, by = d)
    return(S)
  }

  # ------------------------------------------------------------------------- #
  # Step 3: Join the l + 1 blocks
  # ------------------------------------------------------------------------- #
  # Pre-allocate y_star vector to store the bootstrapped data
  y_star <- vector("numeric", length =  b * (l + 1))

  # Pre-allocate a list to store indices for each bootstrap block
  ids_y_bootstrap_list <- vector("list", length = length(ind_set_t))

  # Initialize a counter to track the position in the list of indices
  counter <- 1

  # Loop through each index in ind_set_t to perform block bootstrap
  for (i in seq_along(ind_set_t)) {
    ind_t <- ind_set_t[i]

    # Randomly sample
    k_t <- sample(calc_S(t = ind_t, n = n, b = b, d = d), 1)

    # Define the range of indices for the current bootstrap block
    start_idx <- k_t
    end_idx <- k_t + b - 1

    # Store the indices of the bootstrap block in the list
    ids_y_bootstrap_list[[counter]] <- start_idx:end_idx

    # Update the y_star vector with the sampled values for this block
    y_star[ind_t:(ind_t + b - 1)] <- y[start_idx:end_idx]

    # Increment the counter for the next block
    counter <- counter + 1
  }

  # ------------------------------------------------------------------------- #
  # Output: Post-processing and storage
  # ------------------------------------------------------------------------- #

  # Trim y_star to match the original length of the data
  y_star <- y_star[1:n]
  y_star_ts <- ts(y_star, frequency = seasonal_period)

  # Flatten the list of bootstrap indices into a single vector
  ids_y_bootstrap_aux <- unlist(ids_y_bootstrap_list)

  # Trim the bootstrap indices to match the original length of the data
  ids_y_bootstrap <- ids_y_bootstrap_aux[1:n]

  output_list <- list()
  output_list$y_star_ts <- y_star_ts
  output_list$ids_y_bootstrap <- ids_y_bootstrap


  return(output_list)
}

# ## Examples
# y <- seq(0.001, 0.306, by = 0.001)
# block_length <- 24
# seed <- 0
# seasonal_period <- 12
#
# ## ------------------------------------------------------------------------ #
# ## Example 1
# ## ------------------------------------------------------------------------ #
# block_bootstrap_seasonal <- bootstrap_seasonal(
#   y = y,
#   block_length = block_length,
#   seasonal_period = seasonal_period
#   )
#
# block_bootstrap_seasonal$ids_y_bootstrap
#
# ids_y_bootstrap <- block_bootstrap_seasonal$ids_y_bootstrap
#
# length(y[ids_y_bootstrap])
# length(y[ids_y_bootstrap])
#
# all(y[ids_y_bootstrap] == block_bootstrap_seasonal$y_star_ts)
#
# y[ids_y_bootstrap] == block_bootstrap_seasonal$y_star_ts
#
#
# ## ------------------------------------------------------------------------ #
# ##  Example 2
# ## ------------------------------------------------------------------------ #
# set.seed(0)
# numobs <- length(y)
# b <- 24
# d <- 12
# l <- floor(numobs / b)
# ind.loop <- seq(from = 1, to = l * b + 1, by = b)
# ind.loop.len <- l + 1
# y.boot.raw <- vector("numeric", length = b * (l + 1))
# for (i in 1:ind.loop.len) {
#   t <- ind.loop[i]
#   R1 <- floor((t - 1) / d)
#   R2 <- floor((numobs - b - t) / d)
#   ind.set <- seq(from = t - d * R1, to = t + d * R2, by = d)
#   ind.star <- sample(ind.set, size = 1)
#   block.random <- y[ind.star:(ind.star + 23)]
#   y.boot.raw[t:(t + (b - 1))] <- block.random
# }
#
# y.boot <- y.boot.raw[1:numobs]
#
# all(y[ids_y_bootstrap] == block_bootstrap_seasonal$y_star_ts)
# y.boot == block_bootstrap_seasonal$y_star_ts
