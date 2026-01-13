#' Imputation of Missing Values Using Random Forest Imputation
#'
#' `impute_randomforest` performs imputation for missing values in the data using the random
#' forest-based method implemented in the `missForest` package.
#'
#' The function imputes missing values by building random forests, where missing values are
#' predicted based on other available values within the dataset. For each variable with missing
#' data, the function trains a random forest model using the available (non-missing) data in
#' that variable, and subsequently predicts the missing values.
#'
#' In addition to the imputed values, users can choose to retain additional columns from the
#' original input data frame that were not part of the imputation process.
#'
#' This function allows passing additional parameters to the underlying `missForest` function,
#' such as controlling the number of trees used in the random forest models or specifying the
#' stopping criteria. For a full list of parameters, refer to the `missForest` documentation.
#'
#' To enable parallelisation, ensure that the `doParallel` package is installed and loaded:
#' ```
#' install.packages("doParallel")
#' library(doParallel)
#' ```
#' Then register the desired number of cores for parallel processing:
#' ```
#' registerDoParallel(cores = 6)
#' ```
#' To leverage parallelisation during the imputation, pass `parallelize = "variables"`
#' as an argument to the `missForest` function.
#'
#' @param data A data frame that contains the input variables. This should include columns for
#' the sample names, precursor or peptide identifiers, and intensity values.
#' @param sample A character column in the `data` data frame that contains the sample names.
#' @param grouping A character column in the `data` data frame that contains the precursor or
#' peptide identifiers.
#' @param intensity_log2 A numeric column in the `data` data frame that contains the intensity
#' values.
#' @param retain_columns A character vector indicating which columns should be retained from the
#' input data frame. These columns will be preserved in the output alongside the imputed values.
#' By default, no additional columns are retained (`retain_columns = NULL`), but specific
#' columns can be retained by providing their names as a vector.
#' @param ... Additional parameters to pass to the `missForest` function. These parameters
#' can control aspects such as the number of trees (`ntree`) and the stopping criteria
#' (`maxiter`).
#'
#' @return A data frame that contains an `imputed_intensity` column with the imputed values
#' and an `imputed` column indicating whether each value was imputed (`TRUE`) or not
#' (`FALSE`), in addition to any columns retained via `retain_columns`.
#'
#' @author Elena Krismer
#'
#' @references
#' Stekhoven, D.J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation
#' for mixed-type data. Bioinformatics, 28(1), 112-118. https://doi.org/10.1093/bioinformatics/btr597
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stringr str_sort
#' @importFrom tidyr pivot_wider pivot_longer
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
#' # Perform imputation
#' data_imputed <- impute_randomforest(
#'   data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity_missing
#' )
#'
#' head(data_imputed, n = 24)
impute_randomforest <- function(
  data,
  sample,
  grouping,
  intensity_log2,
  retain_columns = NULL,
  ...
) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    message(strwrap("Package \"missForest\" is needed for this function to work. Please install it.",
      prefix = "\n", initial = ""
    ))
    return(invisible(NULL))
  }

  input <- data %>%
    distinct({{ sample }}, {{ grouping }}, {{ intensity_log2 }})

  # Pivot to wide format and remove the sample column
  data_wide_pre <- input %>%
    tidyr::pivot_wider(
      id_cols = {{ sample }},
      names_from = {{ grouping }},
      values_from = {{ intensity_log2 }},
      names_repair = "minimal"
    )

  # get samples names in the correct order to add back later
  sample_names <- data_wide_pre %>%
    dplyr::pull({{ sample }})

  data_wide <- data_wide_pre %>%
    dplyr::select(-{{ sample }}) %>% # Exclude the sample column for imputation
    as.data.frame()

  # Perform the random forest imputation
  data_imputed_rf <- missForest::missForest(data_wide, ...)

  # Get the imputed values
  data_imputed <- data_imputed_rf$ximp

  # Add the sample column back to the imputed data
  data_imputed_complete <- data_imputed %>%
    mutate({{ sample }} := sample_names)

  # Convert back to long format
  data_imputed_long <- data_imputed_complete %>%
    tidyr::pivot_longer(
      cols = -{{ sample }},
      names_to = rlang::as_name(rlang::enquo(grouping)),
      values_to = "imputed_intensity"
    ) %>%
    dplyr::left_join(input, by = c(rlang::as_name(rlang::enquo(grouping)), rlang::as_name(rlang::enquo(sample)))) %>%
    dplyr::mutate(imputed = is.na({{ intensity_log2 }}) & !is.na(.data$imputed_intensity))

  # Join the retained columns if specified
  if (!missing(retain_columns)) {
    result <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(data_imputed_long)[!colnames(data_imputed_long) %in% c(
          "imputed_intensity",
          "imputed"
        )]
      ) %>%
      dplyr::distinct() %>%
      dplyr::right_join(data_imputed_long, by = colnames(data_imputed_long)[!colnames(data_imputed_long) %in% c(
        "imputed_intensity",
        "imputed"
      )]) %>%
      dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))
  } else {
    result <- data_imputed_long
  }

  return(result)
}
