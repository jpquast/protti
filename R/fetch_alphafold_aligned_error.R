#' Fetch AlphaFold aligned error
#'
#' Fetches the aligned error for AlphaFold predictions for provided proteins.
#' The aligned error is useful for assessing inter-domain accuracy. In detail it
#' represents the expected position error at residue x (scored residue), when
#' the predicted and true structures are aligned on residue y (aligned residue).
#'
#' @param uniprot_ids a character vector of UniProt identifiers for which predictions
#' should be fetched.
#' @param error_cutoff a numeric value specifying the maximum position error (in Angstroms) that should be retained.
#' setting this value to a low number reduces the size of the retrieved data. Default is 20.
#' @param timeout a numeric value specifying the time in seconds until the download times out.
#' The default is 30 seconds.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. The default is 1.
#' @param return_data_frame a logical value; if `TRUE` a data frame instead of a list
#' is returned. It is recommended to only use this if information for few proteins is retrieved.
#' Default is `FALSE`.
#' @param show_progress a logical value; if `TRUE` a progress bar will be shown.
#' Default is `TRUE`.
#'
#' @return A list that contains aligned errors for AlphaFold predictions. If return_data_frame is
#' TRUE, a data frame with this information is returned instead. The data frame contains the
#' following columns:
#'
#' * scored_residue: The error for this position is calculated based on the alignment to the
#' aligned residue.
#' * aligned_residue: The residue that is aligned for the calculation of the error of the scored
#' residue
#' * error: The predicted aligned error computed by alpha fold.
#' * accession: The UniProt protein identifier.
#'
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @importFrom tibble tibble
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' aligned_error <- fetch_alphafold_aligned_error(
#'   uniprot_ids = c("F4HVG8", "O15552"),
#'   error_cutoff = 5,
#'   return_data_frame = TRUE
#' )
#'
#' head(aligned_error, n = 10)
#' }
fetch_alphafold_aligned_error <- function(uniprot_ids = NULL,
                                          error_cutoff = 20,
                                          timeout = 30,
                                          max_tries = 1,
                                          return_data_frame = FALSE,
                                          show_progress = TRUE) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  # remove NAs from UniProt IDs
  uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]

  batches <- purrr::map(
    .x = uniprot_ids,
    .f = ~ paste0("https://alphafold.ebi.ac.uk/files/AF-", .x, "-F1-predicted_aligned_error_v4.json")
  )

  names(batches) <- uniprot_ids

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(batches),
      format = "Fetching AlphaFold aligned error [:bar] :current/:total (:percent) :eta"
    )
  }

  query_result <- purrr::map(
    .x = batches,
    .f = ~ {
      # query information from database
      query <- try_query(.x,
        type = "application/json",
        timeout = timeout,
        max_tries = max_tries
      )

      if (show_progress == TRUE) {
        pb$tick()
      }

      if (methods::is(query, "character")) {
        return(invisible(query))
      }

      initial_result <- query[[1]][["predicted_aligned_error"]]

      final_result <- purrr::imap_dfr(
        .x = initial_result,
        .f = ~ {
          tibble::tibble(
            aligned_residue = .y,
            scored_residue = 1:length(initial_result),
            error = as.numeric(.x)
          ) %>%
            dplyr::filter(error <= error_cutoff)
        }
      )
      final_result
    }
  )

  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>%
    purrr::keep(.p = ~ is.character(.x))

  if (length(error_list) != 0) {
    error_table <- tibble::tibble(
      id = names(error_list),
      error = unlist(error_list)
    ) %>%
      dplyr::distinct()

    if (any(stringr::str_detect(error_table$error, pattern = "Timeout"))){
      message('Consider increasing the "timeout" or "max_tries" argument. \n')
    }

    message("The following IDs have not been retrieved correctly:")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output

  query_result <- query_result %>%
    purrr::keep(.p = ~ !is.character(.x))

  if (length(query_result) == 0) {
    message("No valid information could be retrieved!")
    return(invisible(NULL))
  }

  if (return_data_frame == FALSE) {
    return(query_result)
  } else {
    query_result_df <- purrr::map2_dfr(
      .x = query_result,
      .y = names(query_result),
      .f = ~ {
        .x %>%
          dplyr::mutate(accession = .y)
      }
    )
    return(query_result_df)
  }
}
