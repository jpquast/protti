#' Perform Welch's t-test
#'
#' Performs a Welch's t-test and calculates p-values between two groups.
#'
#' @param mean1 Vector containing the means of group1.
#' @param mean2 Vector containing the means of group2.
#' @param sd1 Vector containing the standard deviations of group1.
#' @param sd2 Vector containing the standard deviations of group2.
#' @param n1 Vector containing the number of replicates used for the calculation of each mean and standard deviation of group1.
#' @param n2 Vector containing the number of replicates used for the calculation of each mean and standard deviation of group2.
#'
#' @return A data frame that contains the calculated differences of means, standard error, t statistic and p-values.
#' @importFrom stats pt
#' @export
#'
#' @examples
#' \dontrun{
#' ttest_protti(
#' mean1 = 10,
#' mean2 = 15.5,
#' sd1 = 1,
#' sd2 = 0.5,
#' n1 = 3,
#' n2 = 3)
#' }
ttest_protti <- function(mean1, mean2, sd1, sd2, n1, n2) {
  std_error <- sqrt( (sd1^2/n1) + (sd2^2/n2))
  #Welch-Satterwhite equation to estimate the degrees of freedom
  df <- ((sd1^2/n1) + (sd2^2/n2))^2 / (sd1^4/(n1^2 * (n1-1)) + sd2^4/(n2^2 * (n2-1)))
  #t statistic calculation
  t <- (mean1 - mean2)/std_error
  result <- data.frame(cbind(mean1 - mean2, std_error, t, 2*pt(-abs(t),df)))
  colnames(result) <- c("difference_of_means", "std_error", "t_statistic", "p_value")
  return(result)
}