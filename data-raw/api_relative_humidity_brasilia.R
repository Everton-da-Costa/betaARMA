# =============================================================================
# data-raw/data_brasilia.R
#
# Provenance
# ----------
# Source      : NASA Prediction of Worldwide Energy Resources (POWER) project
#               Temporal/Daily endpoint, parameter RH2M
#               (Relative Humidity at 2 meters above the surface)
#               https://power.larc.nasa.gov/
# Variable    : Monthly relative humidity in Brasília, Brazil,
#               expressed as a proportion (0 to 1).
# Coordinates : latitude = -15.7797, longitude = -47.9257
#               From: https://power.larc.nasa.gov/data-access-viewer/
# Period      : January 1999 -- June 2024 (306 monthly observations)
# Method      : Daily values downloaded via the NASA POWER API and averaged
#               to monthly proportions.
#
# Processing notes
# ----------------
# 1. Daily relative humidity values (%) are downloaded from the NASA POWER
#    API for the coordinates of Brasília.
# 2. Values are converted to proportions by dividing by 100.
# 3. Daily proportions are averaged within each calendar month to obtain
#    monthly series.
# 4. All values are strictly in (0, 1) — no replacements required.
#
# Reference
# ---------
# Cribari-Neto, F., Costa, E., & Fonseca, R. V. (2025). Numerical stability
# enhancements in beta autoregressive moving average model estimation.
# Brazilian Journal of Probability and Statistics, 39(4), 410-437.
# https://doi.org/10.1214/25-BJPS645
# =============================================================================

library(httr)      # For making API requests
library(jsonlite)  # For handling JSON data
library(dplyr)     # For data manipulation
library(lubridate) # For date operations
library(zoo)       # For yearmon time index

# -----------------------------------------------------------------------------
# 1. API REQUEST — Daily relative humidity for Brasília
# -----------------------------------------------------------------------------

# Coordinates for Brasília, Brazil
latitude  <- -15.7797
longitude <- -47.9257

# NASA POWER daily endpoint
base_url   <- "https://power.larc.nasa.gov/api/temporal/daily/point"
parameters <- "RH2M"       # Relative Humidity at 2 meters
start_date <- "19990101"   # January 1, 1999  (format: YYYYMMDD)
end_date   <- "20240601"   # June 1, 2024     (format: YYYYMMDD)

url <- paste0(
  base_url,
  "?parameters=", parameters,
  "&community=AG",
  "&latitude=",   latitude,
  "&longitude=",  longitude,
  "&start=",      start_date,
  "&end=",        end_date,
  "&format=JSON"
)

response <- GET(url)

if (response$status_code != 200) {
  stop("NASA POWER API request failed. Status code: ", response$status_code)
}

data <- fromJSON(content(response, as = "text"))

# -----------------------------------------------------------------------------
# 2. PROCESS DAILY DATA
# -----------------------------------------------------------------------------

# Extract daily RH2M values (named list: date string -> value)
rh2m_daily <- data$properties$parameter$RH2M

# Build daily data frame: date and proportion (divided by 100)
brasilia_daily_df <- data.frame(
  time = as.Date(names(rh2m_daily), format = "%Y%m%d"),
  y    = unlist(rh2m_daily) / 100   # Convert % to proportion
)

rownames(brasilia_daily_df) <- NULL

# -----------------------------------------------------------------------------
# 3. AGGREGATE DAILY TO MONTHLY
# -----------------------------------------------------------------------------

brasilia_monthly_df <- brasilia_daily_df %>%
  mutate(
    year  = year(time),
    month = month(time)
  ) %>%
  group_by(year, month) %>%
  summarise(y = mean(y, na.rm = TRUE), .groups = "drop") %>%
  data.frame()

# -----------------------------------------------------------------------------
# 4. BUILD ts OBJECT AND DATA FRAME
# -----------------------------------------------------------------------------

# ts object: monthly, starting January 1999
brasilia_ts <- ts(
  round(brasilia_monthly_df$y, 4),
  start     = c(1999, 1),
  frequency = 12
)

# data frame with yearmon time column (used in vignette Feature Engineering)
brasilia_df <- data.frame(
  time = as.yearmon(time(brasilia_ts)),
  y    = as.numeric(brasilia_ts)
)

# -----------------------------------------------------------------------------
# 5. VALIDATION CHECKS
# -----------------------------------------------------------------------------
stopifnot(
  "Values outside (0, 1)" =
    all(brasilia_ts > 0 & brasilia_ts < 1),
  "Missing values in series" =
    !anyNA(brasilia_ts),
  "Series does not start in January 1999" =
    start(brasilia_ts)[1] == 1999 && start(brasilia_ts)[2] == 1,
  "Series does not end in June 2024" =
    end(brasilia_ts)[1] == 2024 && end(brasilia_ts)[2] == 6
)

message(
  "brasilia_ts: ", length(brasilia_ts), " observations, ",
  format(min(brasilia_df$time)), " to ", format(max(brasilia_df$time))
)

# -----------------------------------------------------------------------------
# 6. SAVE TO PACKAGE data/
# -----------------------------------------------------------------------------
usethis::use_data(brasilia_ts, brasilia_df, overwrite = TRUE)