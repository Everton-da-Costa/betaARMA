#' Monthly Relative Humidity in Brasília, Brazil (Data Frame)
#'
#' @description
#' A time series dataset containing the monthly average relative humidity in
#' Brasília, Brazil, expressed as a proportion (ranging from 0 to 1). [cite: 8] The data
#' covers a period with pronounced seasonal differences due to the city's
#' tropical savanna climate, which features extremely dry winter months. [cite: 30]
#'
#' @format A data frame with 2 variables:
#' \describe{
#'   \item{time}{Month and year of observation (yearmon format)}
#'   \item{y}{Relative humidity as a proportion (between 0 and 1)}
#' }
#' @source The data were obtained from the NASA (National Aeronautics and Space
#' Administration) Prediction of Worldwide Energy Resources (POWER) project. [cite: 38]
"brasilia_df"

#' Brasília Relative Humidity Monthly Time Series
#'
#' Monthly relative humidity data from Brasília, Brazil, as a time series object,
#' starting from January 1999. [cite: 37]
#'
#' @format A time series object (`ts`) with 306 monthly observations from
#' January 1999 to June 2024. [cite: 37]
#'
#' @source The data were obtained from the NASA (National Aeronautics and Space
#' Administration) Prediction of Worldwide Energy Resources (POWER) project. [cite: 38]
"brasilia_ts"
