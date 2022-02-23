#' Query from URL
#'
#' Downloads data table from URL. If an error occurs during the query (for example due to no
#' connection) the function waits 3 seconds and tries again. If no result could be obtained
#' after the given number of tries a message indicating the problem is returned.
#'
#' @param url a character value of an URL to the website that contains the table that should be
#' downloaded.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs.
#' @param silent a logical value that specifies if individual messages are printed after each try
#' that failed.
#' @param type a character value that specifies the type of data at the target URL. Options are
#' all options that can be supplied to httr::content, these include e.g.
#' "text/tab-separated-values", "application/json" and "txt/csv". Default is "text/tab-separated-values".
#' Default is "tab-separated-values".
#' @param ... other parameters supplied to the parsing function used by httr::content.
#'
#' @importFrom curl has_internet
#'
#' @return A data frame that contains the table from the url.
try_query <-
  function(url, max_tries = 5, silent = TRUE, type = "text/tab-separated-values", ...) {
    if (!requireNamespace("httr", quietly = TRUE)) {
      message("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }

    if (!curl::has_internet()) {
      if (!silent) message("No internet connection.")
      return(invisible("No internet connection"))
    }

    query_result <- NULL
    try_n <- 0
    while (class(query_result) != "response" & try_n < max_tries) {
      query_result <- tryCatch(httr::GET(url, httr::timeout(30)),
        error = function(e) conditionMessage(e),
        warning = function(w) conditionMessage(w)
      )
      try_n <- try_n + 1
      if (!silent & class(query_result) != "response") {
        message(paste0(try_n, ". try failed!"))
      }
      if (class(query_result) != "response") {
        Sys.sleep(3)
      }
    }

    if (class(query_result) != "response") {
      if (!silent) message(query_result)
      return(invisible(query_result))
    }

    if (httr::http_error(query_result)) {
      if (!silent) httr::message_for_status(query_result)
      return(invisible(httr::http_status(query_result)$message))
    }

    result <- suppressMessages(httr::content(query_result, type = type, encoding = "UTF-8", ...))

    return(result)
  }
