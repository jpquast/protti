#' Query from URL
#'
#' Downloads data table from URL. If an error occurs during the query (for example due to no
#' connection) the function waits 3 seconds and tries again. If no result could be obtained
#' after the given number of tries a message indicating the problem is returned.
#'
#' @param url a character value of an URL to the website that contains the table that should be
#' downloaded.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. By default does not try again if database times out see `try_if_timout`.
#' @param try_if_timeout a logical value specifying if the function should try again in case of a timeout.
#' The default is `FALSE`, which prevents queries from taking very long even though the database is not
#' responding at the time. In some cases it makes sense to set this to `TRUE`, increase the number
#' of tries and reduce the timeout time.
#' @param silent a logical value that specifies if individual messages are printed after each try
#' that failed.
#' @param type a character value that specifies the type of data at the target URL. Options are
#' all options that can be supplied to httr::content, these include e.g.
#' "text/tab-separated-values", "application/json" and "txt/csv". Default is "text/tab-separated-values".
#' Default is "tab-separated-values".
#' @param timeout a numeric value that specifies the maximum request time. Default is 60 seconds.
#' @param accept a character value that specifies the type of data that should be sent by the API if
#' it uses content negotiation. The default is NULL and it should only be set for APIs that use
#' content negotiation.
#' @param ... other parameters supplied to the parsing function used by httr::content.
#'
#' @importFrom curl has_internet
#' @importFrom httr GET timeout http_error message_for_status http_status content accept
#'
#' @return A data frame that contains the table from the url.
try_query <-
  function(url, max_tries = 5, try_if_timeout = FALSE, silent = TRUE, type = "text/tab-separated-values", timeout = 60, accept = NULL, ...) {
    # Check if there is an internet connection first
    if (!curl::has_internet()) {
      if (!silent) message("No internet connection.")
      return(invisible("No internet connection"))
    }

    query_result <- "empty"
    try_n <- 0
    # Note: The handling of retries in case of a timeout could be adjusted in the future.
    # For now I introduced try_if_timeout, but potentially it makes sense to always retry even
    # in case of a timeout. Then however, the number of retries needs to be adjusted for some
    # functions that retrieve a lot of data.
    while (!is(query_result, "response") &
      try_n < max_tries &
      # this ifelse stops requery if they timeout except for if try_if_timeout is TRUE
      !ifelse(is(query_result, "character"),
        stringr::str_detect(query_result, pattern = "Timeout was reached") & !try_if_timeout,
        FALSE
      )
      ) {
      if (!missing(accept)) {
        # with accept set
        query_result <- tryCatch(httr::GET(url, httr::accept(accept), httr::timeout(timeout)),
          error = function(e) conditionMessage(e),
          warning = function(w) conditionMessage(w)
        )
      } else {
        # without accept set (the usual case)
        query_result <- tryCatch(httr::GET(url, httr::timeout(timeout)),
          error = function(e) conditionMessage(e),
          warning = function(w) conditionMessage(w)
        )
      }
      try_n <- try_n + 1
      if (!silent & !is(query_result, "response")) {
        message(paste0(try_n, ". try failed!"))
      }
      if (!is(query_result, "response")) {
        Sys.sleep(3)
      }
    }

    # Check again if there is an internet connection, if not then the correct error is returned
    if (!curl::has_internet()) {
      if (!silent) message("No internet connection.")
      return(invisible("No internet connection"))
    }

    if (!is(query_result, "response")) {
      if (!silent) message(query_result)
      return(invisible(query_result))
    }

    if (httr::http_error(query_result)) {
      if (!silent) httr::message_for_status(query_result)
      return(invisible(httr::http_status(query_result)$message))
    }

    # Record readr progress variable to set back later
    readr_show_progress <- getOption("readr.show_progress")
    on.exit(options(readr.show_progress = readr_show_progress))
    # Change variable to not show progress if readr is used
    options(readr.show_progress = FALSE)

    result <- suppressMessages(httr::content(query_result, type = type, encoding = "UTF-8", ...))

    return(result)
  }
