#' Imputation of Missing Values Using Random Forest Imputation
#'
#' \code{impute_randomforest} performs imputation for missing values in the data using the random
#' forest-based method implemented in the \code{missForest} package.
#'
#' The function imputes missing values by building random forests, where missing values are
#' predicted based on other available values within the dataset. For each variable with missing
#' data, the function trains a random forest model using the available (non-missing) data in
#' that variable, and subsequently predicts the missing values.
#'
#' In addition to the imputed values, users can choose to retain additional columns from the
#' original input data frame that were not part of the imputation process.
#'
#' This function allows passing additional parameters to the underlying \code{missForest} function,
#' such as controlling the number of trees used in the random forest models or specifying the
#' stopping criteria. For a full list of parameters, refer to the \code{missForest} documentation.
#'
#' To enable parallelization, ensure that the `doParallel` package is installed and loaded:
#' ```
#' install.packages("doParallel")
#' library(doParallel)
#' ```
#' Then register the desired number of cores for parallel processing:
#' ```
#' registerDoParallel(cores = 6)
#' ```
#' To leverage parallelization during the imputation, pass `parallelize = "variables"`
#' as an argument to the `missForest` function.
# `
#'
#' Stekhoven, D.J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation
#' for mixed-type data. Bioinformatics, 28(1), 112-118. https://doi.org/10.1093/bioinformatics/btr597
#'
#' @param data A data frame that contains the input variables. This should include columns for
#' the sample names, precursor or peptide identifiers, and intensity values.
#' @param sample A character column in the \code{data} data frame that contains the sample names.
#' @param grouping A character column in the \code{data} data frame that contains the precursor or
#' peptide identifiers.
#' @param intensity_log2 A numeric column in the \code{data} data frame that contains the intensity
#' values.
#' @param retain_columns A character vector indicating which columns should be retained from the
#' input data frame. These columns will be preserved in the output alongside the imputed values.
#' By default, no additional columns are retained (\code{retain_columns = NULL}), but specific
#' columns can be retained by providing their names as a vector.
#' @param ... Additional parameters to pass to the \code{missForest} function. These parameters
#' can control aspects such as the number of trees (\code{ntree}) and the stopping criteria
#' (\code{maxiter}).
#'
#' @return A data frame that contains an \code{imputed_intensity} column with the imputed values
#' and an \code{imputed} column indicating whether each value was imputed (\code{TRUE}) or not
#' (\code{FALSE}), in addition to any columns retained via \code{retain_columns}.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stringr str_sort
#' @import missForest
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
    ...) {
  # Convert inputs to symbols
  sample_sym <- rlang::ensym(sample)
  grouping_sym <- rlang::ensym(grouping)
  intensity_sym <- rlang::ensym(intensity_log2)

  # Pivot to wide format and remove the sample column
  data_wide <- data %>%
    tidyr::pivot_wider(
      id_cols = !!sample_sym,
      names_from = !!grouping_sym,
      values_from = !!intensity_sym
    ) %>%
    dplyr::select(-!!sample_sym) # Exclude the sample column for imputation

  # Convert to numeric and suppress warnings for coercion
  data_wide <- suppressWarnings(data.frame(lapply(data_wide, function(x) as.numeric(as.character(x)))))

  # Perform the random forest imputation
  data_imputed_rf <- missForest::missForest(data_wide, ...)

  # Get the imputed values
  data_imputed <- data_imputed_rf$ximp

  # Add the sample column back to the imputed data
  data_imputed[[as.character(sample_sym)]] <- data %>%
    tidyr::pivot_wider(
      id_cols = !!sample_sym,
      names_from = !!grouping_sym,
      values_from = !!intensity_sym
    ) %>%
    dplyr::pull(!!sample_sym)

  # Convert back to long format
  data_imputed <- data_imputed %>%
    as.data.frame() %>%
    tidyr::pivot_longer(
      cols = -as.character(sample_sym),
      names_to = as.character(grouping_sym),
      values_to = as.character(intensity_sym)
    )

  # Join the retained columns if specified
  if (!is.null(retain_columns)) {
    result <- data_imputed %>%
      dplyr::left_join(
        data %>%
          dplyr::select(retain_columns, !!sample_sym),
        by = as.character(sample_sym)
      )
  } else {
    result <- data_imputed
  }

  return(result)
}

#' @references
#' Stekhoven, D.J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation
#' for mixed-type data. Bioinformatics, 28(1), 112-118. https://doi.org/10.1093/bioinformatics/btr597
