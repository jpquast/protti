#' Sampling of values for imputation
#'
#' \code{calculate_imputation} is a helper function that is used in the \code{impute} function. Depending on the type of missingness and method, it samples values from a normal distribution that can be used for the imputation. Note: The input intensities should be log2 transformed.
#'
#' @param n_replicates number of replicates for each condition
#' @param min minimal intensity value of the precursor/peptide. Is only required if \code{method = "ludovic"} and \code{missingness = "MNAR"}.
#' @param noise noise value for the precursor/peptide. Is only required if \code{method = "noise"} and \code{missingness = "MNAR"}.
#' @param mean mean intensity value of the condition with missing values for a given precursor/peptide. Is only required if \code{missingness = "MAR"}.
#' @param sd mean of the standard deviation of all conditions for a given precursor/peptide. 
#' @param missingness the missingness type of the data determines how values for imputation are sampled. This can be \code{"MAR"} or \code{"MNAR"}. 
#' @param method the method to be used for imputation. For \code{method = "ludovic"}, MNAR missingness is sampled around a value that is three lower (log2) than the lowest intensity value recorded for the precursor/peptide. For \code{method = "noise"}, MNAR missingness is sampled around the noise value for the precursor/peptide.
#' @param skip_log2_transform_error logical, if FALSE a check is performed to validate that input values are log2 transformed. If input values are > 40 the test is failed and an error is thrown. 
#' 
#' @return A vector of values for the imputation of missing data. The length of the vector depends on the number of replicates.
#'
#' @examples
#' \dontrun{
#' calculate_imputation(
#' n_replicates = 3,
#' min = 25.4,
#' sd = 0.25,
#' missingness = "MNAR",
#' method = "ludovic"
#' )
#' }
calculate_imputation <-
  function(n_replicates,
           min = NULL,
           noise = NULL,
           mean = NULL,
           sd,
           missingness = c("MNAR", "MAR"),
           method = c("ludovic", "noise"),
           skip_log2_transform_error = FALSE)
  {
    set.seed(123)
    if ((ifelse(is.null(min), FALSE, min) > 40 | ifelse(is.null(mean), FALSE, mean) > 40 | ifelse(is.null(noise), FALSE, noise) > 40) & skip_log2_transform_error == FALSE) {
      stop("Input intensities seem not to be log2 transformed. If they are and you want to proceed set the skip_log2_transform_error argument to TRUE. Notice that this function does not give correct results for non-log2 transformed data.")
    }
    if (method == "ludovic") {
      if (missingness == "MNAR") {
        result <- stats::rnorm(n_replicates, mean = min - 3, sd = sd)
      }
      if (missingness == "MAR") {
        result <- stats::rnorm(n_replicates, mean = mean, sd = sd)
      }
    }
    if (method == "noise") {
      if (missingness == "MNAR") {
        result <- stats::rnorm(n_replicates, mean = noise, sd = sd)
      }
      if (missingness == "MAR") {
        result <- stats::rnorm(n_replicates, mean = mean, sd = sd)
      }
    }
    result
  }