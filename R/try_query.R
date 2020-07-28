#' Query from url
#'
#' Downloads data from url using the \code{fread} function from the \code{data.table} package. If an error occures during the query (for example due to no connection) the function waits 3 seconds and tries again.  
#'
#' @param url Url of website that contains the table that should be downloaded.
#' @param max_tries Number of times the function tries to download the data in case an error occures.
#' @param silent logical, if TRUE no individual messages are printed after each try that failed.
#'
#' @return A data frame that contains the table from the url.
#'
#' @examples
#' \dontrun{
#' try_query("http://www.uniprot.org/uniprot/?query=id:P36578&format=tab")
#' }
try_query <-
  function (url, max_tries = 5, silent = TRUE)
  {
    for (i in 1:max_tries) {
      result <- tryCatch({
        utils::read.delim(url, quote = "", stringsAsFactors = FALSE, fill = FALSE)
      },
      error = function(error) {
        NULL
      })
      if (!is.null(result)) {
        return(result)
      } else{
        if (!silent) {message(paste0(i, ". try failed!"))}
        Sys.sleep(3)
      }
    }
    stop(paste(
      "Database not responding. Check internet connection and try again."
    ))
  }