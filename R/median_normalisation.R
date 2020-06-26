#' Median normalisation
#'
#' Median normalise intensity values for all samples of mass spectrometry data. The intensity values are ideally log2 transformed. The output is a new column called normalised_intensity_log2.
#'
#' @param data A data frame containing at least sample names and intensity values.
#' @param file_column_name The name of the column containing the sample names.
#' @param intensity_column_name The name of the column containing the intensity values.
#' @param na.rm A logical indicating wther missing values should be removed.
#'
#' @return A data frame that contains the input data and an additional column with normalised intensity values.
#' @export
#'
#' @examples
#' \dontrun{
#' median_normalisation(data, sample_name, intensity_log2, na.rm = TRUE)
#' }
median_normalisation <-
function(data, file_column_name, intensity_column_name, na.rm = FALSE){
  # define global variables to prevent notes in function check
  median_run_intensity <- median_intensity <- NULL
  # main function
  data %>%
    dplyr::group_by({{file_column_name}})%>%
    dplyr::mutate(median_run_intensity = stats::median({{intensity_column_name}}, na.rm = na.rm))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(median_intensity = stats::median(unique(median_run_intensity), na.rm = na.rm))%>%
    dplyr::mutate(normalised_intensity_log2 = {{intensity_column_name}}/median_run_intensity*median_intensity)%>%
    dplyr::select(- median_run_intensity, -median_intensity)
}
