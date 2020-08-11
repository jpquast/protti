#' Query from url
#'
#' Downloads data from url using the \code{fread} function from the \code{data.table} package. If an error occures during the query (for example due to no connection) the function waits 3 seconds and tries again.  
#'
#' @param url The url of the website that contains the table that should be downloaded.
#' @param max_tries Number of times the function tries to download the data in case an error occures.
#' @param silent Logical, if TRUE no individual messages are printed after each try that failed.
#' @param sep The separator of the table at the target url. Default is tab.
#' @param quote Logical, determines if quoting in character vectors should be considered or not. If TRUE, " and ' are escaped with \. 
#'
#' @return A data frame that contains the table from the url.
#'
#' @examples
#' \dontrun{
#' try_query("http://www.uniprot.org/uniprot/?query=id:P36578&format=tab")
#' }
try_query <-
  function (url, max_tries = 5, silent = TRUE, sep = "\t", quote = FALSE)
  {
    for (i in 1:max_tries) {
      result <- tryCatch({
        if(quote == FALSE) {utils::read.delim(url, quote = "", stringsAsFactors = FALSE, fill = FALSE, sep = sep)
        } else {utils::read.delim(url, stringsAsFactors = FALSE, fill = FALSE, sep = sep)}
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