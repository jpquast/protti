#' Query from URL
#'
#' Downloads data table from URL. If an error occurs during the query (for example due to no connection) the function waits 3 seconds
#' and tries again. If no result could be obtained after the given number of tries a message indicating the problem is returned.
#'
#' @param url a character vector of an URL to the website that contains the table that should be downloaded.
#' @param max_tries a numeric vector specifying the number of times the function tries to download the data in case an error occurs.
#' @param silent a logical, if TRUE no individual messages are printed after each try that failed.
#' @param header a logical, indicates if the first row of the data frame contains variable names.
#' @param sep a character vector, specifying the separator of the table at the target URL. Options are "tab-separated-values" or "csv".
#' Default is "tab-separated-values".
#'
#' @importFrom curl has_internet
#'
#' @return A data frame that contains the table from the url.
try_query <-
  function(url, max_tries = 5, silent = TRUE, header = TRUE, sep = "tab-separated-values") {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }

    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
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
      message(query_result)
      return(invisible(NULL))
    }

    if (httr::http_error(query_result)) {
      httr::message_for_status(query_result)
      return(invisible(NULL))
    }

    result <- suppressMessages(httr::content(query_result, type = paste0("text/", sep), encoding = "UTF-8", col_names = header))

    return(result)
  }
