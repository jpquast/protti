#' Intensity normalisation
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `normalise()`.
#' The normalisation method in the new function needs to be provided as an argument.
#'
#' @keywords internal
#' @export
median_normalisation <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "median_normalisation()",
    "normalise()",
    details = "This function has been renamed."
  )

  normalise(...)
}
#' Intensity normalisation
#'
#' Performs normalisation on intensities. For median normalisation the normalised intensity is the
#' original intensity minus the run median plus the global median. This is also the way it is
#' implemented in the Spectronaut search engine.
#'
#' @param data a data frame containing at least sample names and intensity values.
#' @param sample a column in the \code{data} data frame that contains the sample names.
#' @param intensity_log2 a column in the \code{data} data frame that contains the log2 transformed
#' intensity values to be normalised.
#' @param method a character value specifying the method to be used for normalisation. Default
#' is "median".
#'
#' @return A data frame with a column called \code{normalised_intensity_log2} containing the
#' normalised intensity values.
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median
#' @export
#'
#' @examples
#' data <- data.frame(
#'   r_file_name = c("s1", "s2", "s3", "s1", "s2", "s3"),
#'   intensity_log2 = c(18, 19, 17, 20, 21, 19)
#' )
#'
#' normalise(data,
#'   sample = r_file_name,
#'   intensity_log2 = intensity_log2,
#'   method = "median"
#' )
normalise <-
  function(data,
           sample,
           intensity_log2,
           method = "median") {
    if (method == "median") {
      median_normalised <- data %>%
        dplyr::distinct() %>%
        dplyr::mutate(global_median = stats::median({{ intensity_log2 }}, na.rm = TRUE)) %>%
        dplyr::group_by({{ sample }}) %>%
        dplyr::mutate(run_median = stats::median({{ intensity_log2 }}, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(normalised_intensity_log2 = {{ intensity_log2 }} - .data$run_median + .data$global_median) %>%
        dplyr::select(-.data$run_median, -.data$global_median)

      return(median_normalised)
    }
  }
