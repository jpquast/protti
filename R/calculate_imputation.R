#' Sampling of values for imputation
#'
#' \code{calculate_imputation} is a helper function that is used in the \code{impute} function.
#' Depending on the type of missingness and method, it samples values from a normal distribution
#' that can be used for the imputation. Note: The input intensities should be log2 transformed.
#'
#' @param min a numeric value specifying the minimal intensity value of the precursor/peptide.
#' Is only required if \code{method = "ludovic"} and \code{missingness = "MNAR"}.
#' @param noise a numeric value specifying a noise value for the precursor/peptide. Is only
#' required if \code{method = "noise"} and \code{missingness = "MNAR"}.
#' @param mean a numeric value specifying the mean intensity value of the condition with missing
#' values for a given precursor/peptide. Is only required if \code{missingness = "MAR"}.
#' @param sd a numeric value specifying the mean of the standard deviation of all conditions for
#' a given precursor/peptide.
#' @param missingness a character value specifying the missingness type of the data determines
#' how values for imputation are sampled. This can be \code{"MAR"} or \code{"MNAR"}.
#' @param method a character value specifying the method to be used for imputation. For
#' \code{method = "ludovic"}, MNAR missingness is sampled around a value that is three lower
#' (log2) than the lowest intensity value recorded for the precursor/peptide. For
#' \code{method = "noise"}, MNAR missingness is sampled around the noise value for the
#' precursor/peptide.
#' @param skip_log2_transform_error a logical value, if FALSE a check is performed to validate that
#' input values are log2 transformed. If input values are > 40 the test is failed and an error is
#' returned.
#'
#' @return A value sampled from a normal distribution with the input parameters. Method specifics
#' are applied to input parameters prior to sampling.
calculate_imputation <-
  function(min = NULL,
           noise = NULL,
           mean = NULL,
           sd,
           missingness = c("MNAR", "MAR"),
           method = c("ludovic", "noise"),
           skip_log2_transform_error = FALSE) {
    if ((ifelse(is.na(ifelse(is.null(min), 0, min) > 40),
      FALSE,
      ifelse(is.null(min), 0, min) > 40
    ) |
      ifelse(is.na(ifelse(is.null(mean), 0, mean) > 40),
        FALSE,
        ifelse(is.null(mean), 0, mean) > 40
      ) |
      ifelse(is.na(ifelse(is.null(noise), 0, noise) > 40),
        FALSE,
        ifelse(is.null(noise), 0, noise) > 40
      )) &
      skip_log2_transform_error == FALSE) {
      stop(strwrap("Input intensities seem not to be log2 transformed. If they are and you want
                   to proceed set the skip_log2_transform_error argument to TRUE. Notice that
                   this function does not give correct results for non-log2 transformed data.",
        prefix = "\n", initial = ""
      ))
    }
    if (!(missingness %in% c("MNAR", "MAR"))) {
      return(NA)
    }
    if (method == "ludovic") {
      if (missingness == "MNAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = min - 3, sd = sd))
      }
      if (missingness == "MAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = mean, sd = sd))
      }
    }
    if (method == "noise") {
      if (missingness == "MNAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = noise, sd = sd))
      }
      if (missingness == "MAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = mean, sd = sd))
      }
    }
    result
  }
