#' Median normalisation
#'
#' Performs median normalisation on intensities. The normalised intensity is the original intensity minus the run median plus the global median. This is also the way it is implemented in Spectronaut.
#'
#' @param data A data frame containing at least sample names and intensity values.
#' @param sample The name of the column containing the sample names.
#' @param intensity_log2 The name of the column containing the log2 transformed intensity values to be normalised.
#'
#' @return A dataframe with a column called \code{normalised_intensity_log2} containing the normalised intensity values.
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
#' median_normalisation(data,
#'   sample = r_file_name,
#'   intensity_log2 = intensity_log2
#' )
median_normalisation <-
  function(data, sample, intensity_log2) {
    data %>%
      dplyr::distinct() %>%
      dplyr::mutate(global_median = stats::median({{ intensity_log2 }}, na.rm = TRUE)) %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::mutate(run_median = stats::median({{ intensity_log2 }}, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(normalised_intensity_log2 = {{ intensity_log2 }} - .data$run_median + .data$global_median) %>%
      dplyr::select(-.data$run_median, -.data$global_median)
  }
