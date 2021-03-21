#' Imputation of missing values
#'
#' \code{impute} is calculating imputation values for missing data depending on the selected method.
#'
#' @param data a dataframe that is ideally the output from the \code{assign_missingness} function. It should containing at least the
#' input variables. For each "reference_vs_treatment" comparison, there should be the pair of the reference and treatment condition.
#' That means the reference condition should be doublicated once for every treatment.
#' @param sample the name of the column containing the sample names.
#' @param grouping the name of the column containing precursor or peptide identifiers.
#' @param intensity the name of the column containing intensity values. Note: The input intensities should be log2 transformed.
#' @param condition the name of the column containing the conditions.
#' @param comparison the name of the column containing the comparisons of treatment/reference pairs. This is an output of the
#' \code{assign_missingnes} function.
#' @param missingness the name of the column that contains the missingness type of the data determines how values for imputation are sampled.
#' This should at least contain \code{"MAR"} or \code{"MNAR"}.
#' @param noise the name of the column that contains the noise value for the precursor/peptide. Is only required if \code{method = "noise"}.
#' Note: Noise values need to be log2 transformed.
#' @param method the method to be used for imputation. For \code{method = "ludovic"}, MNAR missingness is sampled from a normal distribution
#' around a value that is three lower (log2) than the lowest intensity value recorded for the precursor/peptide and that has a spread of
#' the mean standard deviation for the precursor/peptide. For \code{method = "noise"}, MNAR missingness is sampled from a normal distribution
#' around the mean noise for the precursor/peptide and that has a spread of the mean standard deviation (from each condition) for the precursor/peptide.
#' @param skip_log2_transform_error logical, if FALSE a check is performed to validate that input values are log2 transformed. If input
#' values are > 40 the test is failed and an error is thrown.
#' @param retain_columns A vector indicating if certain columns should be retained from the input data frame. Default is not retaining
#' additional columns \code{retain_columns = NULL}. Specific columns can be retained by providing their names (not in quotations marks,
#' just like other column names, but in a vector).
#'
#' @return A data frame that contains an \code{imputed_intensity} and \code{imputed} column in addition to the required input columns.
#' The \code{imputed} column indicates if a value was imputed. The \code{imputed_intensity} column contains imputed intensity values
#' for previously missing intensities.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' impute(
#'   data,
#'   sample = r_file_name,
#'   grouping = eg_precursor_id,
#'   intensity = intensity_log2,
#'   condition = r_condition,
#'   comparison = comparison,
#'   missingness = missingness,
#'   method = "ludovic"
#' )
#' }
impute <- function(data, sample, grouping, intensity, condition, comparison, missingness, noise = NULL, method, skip_log2_transform_error = FALSE, retain_columns = NULL) {
  noise_missing <- missing(noise) # check if argument noise was provided or not
  result <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}, {{ condition }}, {{ comparison }}, {{ missingness }}, {{ noise }}) %>%
    dplyr::group_by({{ grouping }}, {{ condition }}, {{ comparison }}) %>%
    dplyr::mutate(mean = mean({{ intensity }}, na.rm = TRUE)) %>%
    dplyr::mutate(n_replicates = dplyr::n()) %>%
    dplyr::mutate(sd = stats::sd({{ intensity }}, na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::mutate(min = min({{ intensity }}, na.rm = TRUE)) %>%
    dplyr::mutate(noise_mean = ifelse(noise_missing, NA, mean({{ noise }}, na.rm = TRUE))) %>%
    dplyr::mutate(sd = mean(unique(.data$sd), na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}, {{ comparison }}) %>%
    dplyr::mutate(
      impute = dplyr::case_when(
        missingness == "MNAR" ~
        calculate_imputation(
          n_replicates,
          min = min,
          noise = noise_mean,
          mean = mean,
          sd = sd,
          missingness = "MNAR",
          method = method,
          skip_log2_transform_error = skip_log2_transform_error
        ),
        missingness == "MAR" ~
        calculate_imputation(
          n_replicates,
          min = min,
          noise = noise_mean,
          mean = mean,
          sd = sd,
          missingness = "MAR",
          method = method,
          skip_log2_transform_error = skip_log2_transform_error
        )
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imputed_intensity = ifelse(is.na({{ intensity }}) == TRUE, .data$impute, {{ intensity }})) %>%
    dplyr::mutate(imputed = is.na({{ intensity }}) & !is.na(.data$imputed_intensity)) %>%
    dplyr::select(-.data$impute, -.data$mean, -.data$sd, -.data$min, -.data$n_replicates, -.data$noise_mean) %>%
    dplyr::arrange({{ grouping }})

  if (missing(retain_columns)) {
    return(result)
  } else {
    result <- data %>%
      dplyr::select(!!enquo(retain_columns), colnames(result)[!colnames(result) %in% c("imputed_intensity", "imputed")]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c("imputed_intensity", "imputed")]) %>%
      dplyr::arrange({{ grouping }})
    return(result)
  }
}
