#' Query from URL
#'
#' Downloads data table from URL. If an error occurs during the query (for example due to no
#' connection) the function waits 3 seconds and tries again. If no result could be obtained
#' after the given number of tries a message indicating the problem is returned.
#'
#' @param url a character value of an URL to the website that contains the table that should be
#' downloaded.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. Default is 5.
#' @param silent a logical value that specifies if individual messages are printed after each try
#' that failed.
#' @param type a character value that specifies the type of data at the target URL. Options are
#' all options that can be supplied to httr::content, these include e.g.
#' "text/tab-separated-values", "application/json" and "txt/csv". Default is "text/tab-separated-values".
#' @param timeout a numeric value that specifies the maximum request time. Default is 60 seconds.
#' @param accept a character value that specifies the type of data that should be sent by the API if
#' it uses content negotiation. The default is NULL and it should only be set for APIs that use
#' content negotiation.
#' @param ... other parameters supplied to the parsing function used by httr::content.
#'
#' @importFrom curl has_internet
#' @importFrom httr GET timeout http_error message_for_status http_status content accept
#' @importFrom readr read_tsv read_csv
#'
#' @return A data frame that contains the table from the url.
try_query <-
  function(url, max_tries = 5, silent = TRUE, type = "text/tab-separated-values", timeout = 60, accept = NULL, ...) {
    # Check timeout time
    if (timeout < 1) {
      stop("The timeout cannot be less than 1 second.")
    }

    # Check if there is an internet connection first
    if (!curl::has_internet()) {
      if (!silent) message("\nNo internet connection.")
      return(invisible("No internet connection"))
    }

    query_result <- "empty"
    try_n <- 0
    retry_delay <- 3 # seconds
    while (!is(query_result, "response") &
      try_n < max_tries
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

      if (inherits(query_result, "response")) {
        break
      } else if (stringr::str_detect(query_result, "Timeout was reached")) {
        if (!silent) message(paste0("\nAttempt ", try_n + 1, "/", max_tries, " failed: Timeout was reached. Retrying..."))
      } else {
        if (!silent) message(paste0("\nAttempt ", try_n + 1, "/", max_tries, " failed: ", query_result, ". Retrying..."))
      }

      try_n <- try_n + 1

      Sys.sleep(retry_delay)
      retry_delay <- retry_delay * 2 # Exponential backoff
    }

    # Check again if there is an internet connection, if not then the correct error is returned
    if (!curl::has_internet()) {
      if (!silent) message("\nNo internet connection.")
      return(invisible("No internet connection"))
    }

    # If response was an error return that error message
    if (inherits(query_result, "response") && httr::http_error(query_result)) {
      if (!silent) httr::message_for_status(query_result)
      return(invisible(httr::http_status(query_result)$message))
    }

    # Handle other types of errors separately from query errors
    if(inherits(query_result, "character")) {
      if (!silent) message(query_result)
      return(invisible(query_result))
    }

    # Record readr progress variable to set back later
    readr_show_progress <- getOption("readr.show_progress")
    on.exit(options(readr.show_progress = readr_show_progress))
    # Change variable to not show progress if readr is used
    options(readr.show_progress = FALSE)

    # Retrieve the content as raw bytes using httr::content
    raw_content <- httr::content(query_result, type = "raw")
    # Check for gzip magic number (1f 8b) before decompression
    compressed <- length(raw_content) >= 2 && raw_content[1] == as.raw(0x1f) && raw_content[2] == as.raw(0x8b)

    # Check if the content is gzip compressed
    if (!is.null(query_result$headers[["content-encoding"]]) && query_result$headers[["content-encoding"]] == "gzip" && compressed) {
      # Decompress the raw content using base R's `memDecompress`
      decompressed_content <- memDecompress(raw_content, type = "gzip")

      # Convert the raw bytes to a character string
      text_content <- rawToChar(decompressed_content)

      # Read the decompressed content based on the specified type
      if (type == "text/tab-separated-values") {
        result <- readr::read_tsv(text_content, ...)
      } else if (type == "text/html") {
        result <- xml2::read_html(text_content, ...)
      } else if (type == "text/xml") {
        result <- xml2::read_xml(text_content, ...)
      } else if (type == "text/csv" || type == "txt/csv") {
        result <- readr::read_csv(text_content, ...)
      } else if (type == "application/json") {
        result <- jsonlite::fromJSON(text_content, ...) # Using jsonlite for JSON parsing
      } else if (type == "text") {
        result <- text_content # Return raw text as-is
      } else {
        stop("Unsupported content type: ", type)
      }
    } else {
      result <- suppressMessages(httr::content(query_result, type = type, encoding = "UTF-8", ...))
    }

    return(result)
  }
