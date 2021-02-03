#' Dose response curve helper function
#'
#' This function peforms the four-parameter dose response curve fit. It is the helper function for the fit in the \code{fit_drc_4p} function.
#'
#' @param data a data frame that contains at least the dose and response column the model should be fitted to.
#' @param response the name of the column that contains the response values.
#' @param dose the name of the column that contains the dose values.
#' @param log_logarithmic logical indicating if a logarithmic or log-logarithmic model is fitted. 
#' If response values form a symmetric curve for non-log transformed dose values, a logarithmic model instead
#' of a log-logarithmic model should be used. Usually biological dose response data has a log-logarithmic distribution, which is the 
#' reason this is the default. Log-logarithmic models are symmetric if dose values are log transformed. 
#' @param pb progress bar object. This is only necessary if the function is used in an iteration.
#'
#' @return An object of class \code{drc}. If no fit was performed a character vector with content "no_fit".
#'
#' @examples
#' \dontrun{
#' drc_4p(
#'   data,
#'   response = intensity,
#'   dose = concentration
#' )
#' }
drc_4p <- function(data, response, dose, log_logarithmic = TRUE, pb = NULL) {
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop("Package \"drc\" is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!is.null(pb)) pb$tick()
  if (log_logarithmic == TRUE) {
    result <- tryCatch(
      {
        suppressWarnings(drc::drm(
          stats::as.formula(paste(ensym(response), "~", ensym(dose))),
          data = data,
          fct = drc::LL.4(names = c("hill", "min_value", "max_value", "ec_50")),
          control = drc::drmc(otrace = TRUE)
        ))
      },
      error = function(error) {
        c("no_fit")
      }
    )
    return(result)
  }
  if (log_logarithmic == FALSE) {
    result <- tryCatch(
      {
        suppressWarnings(drc::drm(
          stats::as.formula(paste(ensym(response), "~", ensym(dose))),
          data = data,
          fct = drc::L.4(names = c("hill", "min_value", "max_value", "ec_50")),
          control = drc::drmc(otrace = TRUE)
        ))
      },
      error = function(error) {
        c("no_fit")
      }
    )
    return(result)
  }
}