## Code to prepare `brasilia_ts` and `brasilia_df` datasets

# Load necessary libraries for this script.
library(here)
library(zoo)
library(usethis)

## Load data

# Set the system's locale to correctly parse month abbreviations if needed.
# This is good practice if your date column uses non-English month names.
original_locale <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "en_US.UTF-8") # Or "pt_BR.UTF-8" if applicable
on.exit(Sys.setlocale("LC_TIME", original_locale))

# Read the raw CSV data from the data-raw directory.
# Ensure you have a 'brasilia_dataset_raw.csv' file in that location.
brasilia_raw_df <- read.csv(
  here::here("data-raw", "dataset_relative_humidity_brasilia_raw.csv")
)

## Create the Data Frame Object

# Create the data frame object that will be used internally in your package.
# It contains the numeric value and a proper year-month column for plotting.
brasilia_df <- data.frame(
  y = brasilia_raw_df$y,
  time = as.yearmon(brasilia_raw_df$time, format = "%b %Y")
)

## Create the Time Series Object

# Create the final `ts` object that will be exported to the user.
brasilia_ts <- ts(
  data = brasilia_df$y,
  start = c(
    as.numeric(format(min(brasilia_df$time), "%Y")), # Start Year
    as.numeric(format(min(brasilia_df$time), "%m"))  # Start Month
  ),
  frequency = 12
)

# Add metadata attributes to the data frame
attr(brasilia_df, "source") <- "NASA POWER Project"
attr(brasilia_df, "description") <- "Monthly relative humidity in Brasília, Brazil"
attr(brasilia_df, "date_processed") <- Sys.Date()

# Save data objects for the package
usethis::use_data(brasilia_df, overwrite = TRUE)
usethis::use_data(brasilia_ts, overwrite = TRUE)
