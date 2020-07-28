#' Read, clean and convert
#'
#' The function uses the very fast \code{fread} function form the \code{data.table} package. The column names of the resulting data table are made more r-friendly using \code{clean_names} from the \code{janitor} package. It replaces "." and " " with "_" and converts names to lower case which is also known as snake_case. In the end the data table is converted to a tibble.
#'
#' @param filename the path to the file.
#' @param ... additional arguments for the fread function.
#'
#' @importFrom data.table fread
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#'
#' @return A data frame (with class tibble) that contains the content of the specified file.
#' @export
#'
#' @examples
#' \dontrun{
#' read_protti("folder\filename")
#' }
read_protti <-
function(filename, ...){
  data.table::fread(filename, ...)%>%
  janitor::clean_names()%>%
  tibble::as_tibble()
}
