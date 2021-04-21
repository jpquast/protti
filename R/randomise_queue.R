#' Randomise samples in MS queue
#'
#' This function randomises the order of samples in an MS queue. QC and Blank samples are left in place. It is also possible 
#' to randomise only parts of the queue. Before running this make sure to set a specific seed with the \code{set.seed()} function. 
#' This ensures that the randomisation of the result is consistent if the function is run again.
#'
#' @param data optional, a data frame containing a queue. If not provided a queue file can be chosen interactively.
#' @param rows optional, a range of rows in for which samples should be randomized.
#' @param export logical, if TRUE a \code{"randomised_queue.csv"} file will be saved in the working directory. If FALSE a data frame will be returned.
#'
#' @return If \code{export = TRUE} a \code{"randomised_queue.csv"} file will be saved in the working directory. If \code{export = FALSE} a data frame that contains the randomised queue is returned.
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' randomise_queue(data = data, rows = 195:235, export = TRUE)
#' }
randomise_queue <-
  function(data = NULL,
           rows = NULL,
           export = FALSE) {

    # load data interactively if no data is provided in the data argument
    if (is.null(data)) {
      path <- file.choose(".")
      data <- data.table::fread(path, skip = 1)
    }

    unused_data <- NULL

    # if rows argument is provided subset dataset
    if (!is.null(rows)) {
      unused_data <- data %>%
        dplyr::slice(-({{ rows }}))

      data <- data %>%
        dplyr::slice({{ rows }}) %>%
        dplyr::mutate(n = row_number())
    } else {
      data <- data %>%
        dplyr::mutate(n = row_number())
    }

    # separate qc and blank samples
    qc_samples <- data %>%
      dplyr::filter(.data$`Sample Type` == "QC" | .data$`Sample Type` == "Blank")

    # randomise sample order
    samples <- data %>%
      dplyr::filter(.data$`Sample Type` == "Unknown") %>%
      dplyr::mutate(n = sample(n))

    # add back qc and blank samples and also parts of the queue that were not randomised
    result <- qc_samples %>%
      dplyr::bind_rows(samples) %>%
      dplyr::arrange(n) %>%
      dplyr::select(-n)

    result <- unused_data %>%
      dplyr::bind_rows(result)

    # export queue and add bracket type = 4 line in front of data for proper queue import in xcalibur
    if (export == TRUE) {
      write("Bracket Type=4", file = "randomised_queue.csv")
      data.table::fwrite(result, file = "randomised_queue.csv", append = TRUE, col.names = TRUE)
    } else {
      result
    }
  }
