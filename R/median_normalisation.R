#' Median normalisation
#'
#' Performs median normalisation on intensities. The normalised intensity is the original intensity minus the run median plus the global median. This is also the way it is implemented in Spectronaut.
#'
#' @param data A data frame containing at least sample names, grouping variables and intensity values.
#' @param sample The name of the column containing the sample names.
#' @param grouping The name of the column containing the grouping variable - this can be peptides, precursors or proteins.
#' @param log2intensity The name of the column containing the log2 transformed intensity values to be normalised.
#' @param na.rm Logical indicating whether missing values should be removed. Default is TRUE.
#'
#' @return A new column in the original dataframe called \code{normalised_intensity_log2} containing the normalised intensity values.
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median
#' @export
#'
#' @examples
#' \dontrun{
#' median_normalisation(data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' log2intensity = intensity_log2)
#' }
median_normalisation <-
  function(data, sample, grouping, log2intensity, na.rm = TRUE)
  {
    data %>%
      dplyr::filter(!is.na({{log2intensity}})) %>%
      dplyr::mutate(global_median = stats::median({{log2intensity}}), na.rm = na.rm) %>%
      dplyr::group_by({{sample}}) %>%
      dplyr::mutate(run_median = stats::median({{log2intensity}}), na.rm = na.rm) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(normalised_intensity_log2 = {{log2intensity}} - .data$run_median + .data$global_median) %>%
      dplyr::select(-.data$run_median, -.data$global_median)
  }
