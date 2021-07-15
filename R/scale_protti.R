#' Scaling a vector
#'
#' \code{scale_protti} is used to scale a numeric vector either between 0 and 1 or around a centered value using the standard deviation.
#' If a vector containing only one value or repeatedly the same value is provided, 1 is returned as the scaled value for
#' \code{method = "01"} and 0 is returned for \code{metod = "center"}.
#'
#' @param x a numeric vector
#' @param method the method to be used for scaling. "01" scales the vector between 0 and 1. "center" scales the vector equal to \code{base::scale} around a center. This is done by subtracting the mean from every value and then deviding them by the standard deviation.
#'
#' @return A scaled numeric vector.
#' @export
#'
#' @examples
#' scale_protti(c(1, 2, 1, 4, 6, 8), method = "01")
scale_protti <- function(x, method) {
  if (is.numeric(x) == FALSE) {
    stop("x is a ", typeof(x), " vector but needs to be a numeric vector!")
  }
  if (method == "01") {
    result <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

    if ((max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) == 0) {
      result <- rep(1, length(x))
    }
  }
  if (method == "center") {
    result <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

    if (stats::sd(x, na.rm = TRUE) == 0) {
      result <- rep(0, length(x))
    }
  }
  result
}
