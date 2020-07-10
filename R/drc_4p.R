#' Dose response curve helper function
#'
#' This function peforms the four-parameter dose response curve fit. It is the helper function for the fit in the \code{fit_drc_4p} function.
#'
#' @param data a data frame that contains at least the dose and response column the model should be fitted to.
#' @param response the name of the column that contains the response values.
#' @param dose the name of the column that contains the dose values.
#' @param pb progress bar object. This is only necessary if the function is used in an iteration. 
#'
#' @return An object of class \code{drc}. If no fit was performed a character vector with content "no_fit". 
#' @importFrom drc drm
#' @importFrom drc LL.4
#'
#' @examples
#' \dontrun{
#' drc_4p(
#' data,
#' response = intensity,
#' dose = concentration
#' )
#' }
drc_4p <- function(data, response, dose, pb = NULL) {
  if (!is.null(pb)) pb$tick()
  tryCatch({
    suppressWarnings(drc::drm(
      stats::as.formula(paste(ensym(response), "~", ensym(dose))),
      data = data,
      fct = drc::LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
    ))
  }, error = function(error) {
    c("no_fit")
  })
}