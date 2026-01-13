#' Imputation of missing values
#'
#' \code{impute} is calculating imputation values for missing data depending on the selected
#' method.
#'
#' @param data a data frame that is ideally the output from the \code{assign_missingness} function.
#' It should containing at least the input variables. For each "reference_vs_treatment" comparison,
#' there should be the pair of the reference and treatment condition. That means the reference
#' condition should be doublicated once for every treatment.
#' @param sample a character column in the \code{data} data frame that contains the sample names.
#' @param grouping a character column in the \code{data} data frame that contains the precursor or
#' peptide identifiers.
#' @param intensity_log2 a numeric column in the \code{data} data frame that contains the intensity
#' values.
#' @param condition a character or numeric column in the \code{data} data frame that contains the
#' the conditions.
#' @param comparison a character column in the \code{data} data frame that contains the the
#' comparisons of treatment/reference pairs. This is an output of the \code{assign_missingnes}
#' function.
#' @param missingness a character column in the \code{data} data frame that contains the
#' missingness type of the data determines how values for imputation are sampled. This should at
#' least contain \code{"MAR"} or \code{"MNAR"}. Missingness assigned as `NA` will not be imputed.
#' @param noise a numeric column in the \code{data} data frame that contains the noise value for
#' the precursor/peptide. Is only required if \code{method = "noise"}. Note: Noise values need to
#' be log2 transformed.
#' @param method a character value that specifies the method to be used for imputation. For
#' \code{method = "ludovic"}, MNAR missingness is sampled from a normal distribution around a
#' value that is three lower (log2) than the lowest intensity value recorded for the
#' precursor/peptide and that has a spread of the mean standard deviation for the
#' precursor/peptide. For \code{method = "noise"}, MNAR missingness is sampled from a normal
#' distribution around the mean noise for the precursor/peptide and that has a spread of the
#' mean standard deviation (from each condition) for the precursor/peptide. Both methods impute
#' MAR data using the mean and variance of the condition with the missing data.
#' @param skip_log2_transform_error a logical value that determines if a check is performed to
#' validate that input values are log2 transformed. If input values are > 40 the test is failed
#' and an error is returned.
#' @param retain_columns a vector that indicates columns that should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return A data frame that contains an \code{imputed_intensity} and \code{imputed} column in
#' addition to the required input columns. The \code{imputed} column indicates if a value was
#' imputed. The \code{imputed_intensity} column contains imputed intensity values for previously
#' missing intensities.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stringr str_sort
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 10,
#'   frac_change = 0.5,
#'   n_replicates = 4,
#'   n_conditions = 2,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' head(data, n = 24)
#'
#' # Assign missingness information
#' data_missing <- assign_missingness(
#'   data,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   ref_condition = "all",
#'   retain_columns = c(protein, peptide_intensity)
#' )
#'
#' head(data_missing, n = 24)
#'
#' # Perform imputation
#' data_imputed <- impute(
#'   data_missing,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity_missing,
#'   condition = condition,
#'   comparison = comparison,
#'   missingness = missingness,
#'   method = "ludovic",
#'   retain_columns = c(protein, peptide_intensity)
#' )
#'
#' head(data_imputed, n = 24)
impute <- function(data,
                   sample,
                   grouping,
                   intensity_log2,
                   condition,
                   comparison = comparison,
                   missingness = missingness,
                   noise = NULL,
                   method = "ludovic",
                   skip_log2_transform_error = FALSE,
                   retain_columns = NULL) {
  noise_missing <- missing(noise) # check if argument noise was provided or not

  result <- data %>%
    dplyr::distinct(
      {{ sample }},
      {{ grouping }},
      {{ intensity_log2 }},
      {{ condition }},
      {{ comparison }},
      {{ missingness }},
      {{ noise }}
    ) %>%
    dplyr::group_by({{ grouping }}, {{ condition }}, {{ comparison }}) %>%
    dplyr::mutate(mean = mean({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::mutate(sd = stats::sd({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::mutate(min = min({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::mutate(noise_mean = ifelse(noise_missing, NA, mean({{ noise }}, na.rm = TRUE))) %>%
    dplyr::mutate(sd = mean(unique(.data$sd), na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}, {{ sample }}, {{ comparison }}) %>%
    dplyr::mutate(
      impute =
        calculate_imputation(
          min = .data$min,
          noise = .data$noise_mean,
          mean = mean,
          sd = .data$sd,
          missingness = {{ missingness }},
          method = method,
          skip_log2_transform_error = skip_log2_transform_error
        )
    ) %>%
    dplyr::group_by({{ grouping }}, {{ sample }}, {{ missingness }}) %>%
    dplyr::mutate(impute = .data$impute[1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imputed_intensity = ifelse(is.na({{ intensity_log2 }}) == TRUE,
      .data$impute,
      {{ intensity_log2 }}
    )) %>%
    dplyr::mutate(imputed = is.na({{ intensity_log2 }}) & !is.na(.data$imputed_intensity)) %>%
    dplyr::select(-c("impute", "mean", "sd", "min", "noise_mean")) %>%
    dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))

  if (missing(retain_columns)) {
    return(result)
  } else {
    result <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(result)[!colnames(result) %in% c(
          "imputed_intensity",
          "imputed"
        )]
      ) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c(
        "imputed_intensity",
        "imputed"
      )]) %>%
      dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))
    return(result)
  }
}
