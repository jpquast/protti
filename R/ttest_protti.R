#' Perform Welch's t-test
#'
#' Performs a Welch's t-test and calculates p-values between two groups.
#'
#' @param mean1 a numeric vector that contains the means of group1.
#' @param mean2 a numeric vector that contains the means of group2.
#' @param sd1 a numeric vector that contains the standard deviations of group1.
#' @param sd2 a numeric vector that contains the standard deviations of group2.
#' @param n1 a numeric vector that contains the number of replicates used for the calculation of
#' each mean and standard deviation of group1.
#' @param n2 a numeric vector that contains the number of replicates used for the calculation of
#' each mean and standard deviation of group2.
#' @param log_values a logical value that indicates if values are log transformed. This determines
#' how fold changes are calculated. Default is \code{log_values = TRUE}.
#'
#' @return A data frame that contains the calculated differences of means, standard error, t
#' statistic and p-values.
#' @importFrom stats pt
#' @export
#'
#' @examples
#' ttest_protti(
#'   mean1 = 10,
#'   mean2 = 15.5,
#'   sd1 = 1,
#'   sd2 = 0.5,
#'   n1 = 3,
#'   n2 = 3
#' )
ttest_protti <- function(mean1, mean2, sd1, sd2, n1, n2, log_values = TRUE) {
  std_error <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
  # Welch-Satterwhite equation to estimate the degrees of freedom
  df <- ((sd1^2 / n1) + (sd2^2 / n2))^2 / (sd1^4 / (n1^2 * (n1 - 1)) + sd2^4 / (n2^2 * (n2 - 1)))
  # fold change calculation
  if (log_values == TRUE) {
    diff <- mean1 - mean2
  } else {
    diff <- mean1 / mean2
  }
  # t statistic calculation
  t <- (diff) / std_error
  result <- data.frame(cbind(diff, std_error, t, 2 * pt(-abs(t), df)))
  colnames(result) <- c("diff", "std_error", "t_statistic", "pval")
  return(result)
}
