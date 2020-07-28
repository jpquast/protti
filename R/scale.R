#' Scaling a vector
#'
#' \code{scale} is used to scale a numeric vector either between 0 and 1 or around a centered value using the standard deviation.
#'
#' @param x a numeric vector
#' @param method the method to be used for scaling. "01" scales the vector between 0 and 1. "center" scales the vector equal to \code{base::scale} around a center. This is done by subtracting the mean from every value and then deviding them by the standard deviation.
#'
#' @return A scaled numeric vector.
#' @export
#'
#' @examples
#' \dontrun{
#' scale(c(1, 2, 1, 4, 6, 8), method = "01")
#' }
scale <- function(x, method)
{
  if (is.numeric(x) == FALSE) {
    stop("x is a ", typeof(x), " vector but needs to be a numeric vector!")
  }
  if (method == "01") {
    result <- (x - min(x)) / (max(x) - min(x))
  }
  if (method == "center") {
    result <- (x - mean(x)) / stats::sd(x)
  }
  result
}