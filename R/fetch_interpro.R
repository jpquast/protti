#' Fetch domain and residue information from InterPro
#'
#' Fetches either domain level information with e.g. gene ontology annotations or
#' residue level information from the InterPro database.
#'
#' @param uniprot_ids a character vector of UniProt accession numbers.
#' @param return_residue_info a logical value that specifies if either domain or residue information
#' should be returned by the function. The default is `FALSE`.
#' @param manual_query optional, a character value that is a custom query to the InterPro database.
#' This query is pastes after "https://www.ebi.ac.uk/interpro/api/" and before "&page_size=200".
#' The raw data of the query is returned as a list.
#' @param page_size a numeric value that specifies the number of entries that should be retrieved
#' per page of a request. The function anyway iterates through all pages, but this parameters allows you
#' to finetune the number of iterations and thus number of requests to the database. Default is 200.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. The default is 3.
#' @param timeout a numeric value that specifies the maximum request time per try. Default is 20 seconds.
#' @param show_progress a logical value that determines if a progress bar will be shown. Default
#' is TRUE.
#'
#' @return A data frame that contains either domain or residue level information for the provided
#' UniProt IDs.
#' @import dplyr
#' @import progress
#' @import purrr
#' @importFrom tidyr crossing
#' @importFrom stringr str_to_upper
#' @importFrom curl has_internet
#' @export
#'
#' @examples
#' \donttest{
#' uniprot_ids <- c("P36578", "O43324", "Q00796", "O32583")
#'
#' domain_info <- fetch_interpro(uniprot_ids = uniprot_ids)
#'
#' head(domain_info)
#'
#' residue_info <- fetch_interpro(
#'   uniprot_ids = uniprot_ids,
#'   return_residue_info = TRUE
#' )
#'
#' head(residue_info)
#' }
fetch_interpro <- function(uniprot_ids = NULL,
                           return_residue_info = FALSE,
                           manual_query = NULL,
                           page_size = 200,
                           max_tries = 3,
                           timeout = 20,
                           show_progress = TRUE) {
  # Define the progress bar settings
  progress_option_residue <- if (show_progress) {
    list(
      type = "iterator",
      format = "Fetching Residue Positions {cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total}) ETA: {cli::pb_eta}",
      clear = TRUE
    )
  } else {
    NULL
  }
  progress_option_domain <- if (show_progress) {
    list(
      type = "iterator",
      format = "Fetching InterPro Domains {cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total}) ETA: {cli::pb_eta}",
      clear = TRUE
    )
  } else {
    NULL
  }

  # check for internet connection
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  # check if any information for retrieval is provided
  if (missing(manual_query) & missing(uniprot_ids)) {
    stop('Please provide information about what to retrieve to either "uniprot_ids" or "manual_query".')
  }

  # check if uniprot_ids were provided
  if (!missing(uniprot_ids)) {
    uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]
    if (length(uniprot_ids) == 0) {
      stop("Please provide valid UniProt IDs.")
    }
    query_url_part <- paste0("entry/interpro/protein/uniprot/", uniprot_ids, "?")
  }

  # check if manual_query was provided
  # if it was provided it will always be used over protein IDs
  if (!missing(manual_query)) {
    # remove NA values
    query_url_part <- manual_query[!is.na(manual_query)]
    if (length(manual_query) == 0) {
      stop("Please provide a valid manual query based on the interpro documentation.\n
           The query is pasted after the general API URL: https://www.ebi.ac.uk/interpro/api/ \n
           A valid query is eg. entry/interpro/protein/uniprot/O32583?")
    }
  }

  # query data
  if (!return_residue_info | !missing(manual_query)) {
    query_result <- purrr::map(
      .x = query_url_part,
      .f = ~ {
        # Initialize the next page URL with the base query URL
        next_url <- paste0("https://www.ebi.ac.uk/interpro/api/", .x, "&page_size=", page_size)
        all_results <- list()

        # progress bar variable
        total_pages <- 1
        current_page <- 1

        while (!is.null(next_url)) {
          # Call the protti function to query the API
          query <- protti:::try_query(next_url,
            type = "application/json",
            max_tries = max_tries,
            timeout = timeout
          )

          # Store the results from the current page
          all_results[[current_page]] <- query

          # initialize progress bar
          if (show_progress == TRUE & current_page == 1 & length(query_url_part) == 1 & (!is(query, "character") && !is.null(query$count))) {
            pb <- progress::progress_bar$new(
              total = ceiling(query$count / page_size),
              format = "  Fetching Pages [:bar] :current/:total (:percent) :eta"
            )
          }

          # tick progress bar
          if (show_progress == TRUE & length(query_url_part) == 1 & (!is(query, "character") && !is.null(query$count)) && !pb$finished) {
            pb$tick()
          }

          if (!is(query, "character")) {
            # Get the next page URL from the "next" field
            next_url <- query$`next`
          } else {
            next_url <- NULL
          }

          # add to page counter
          current_page <- current_page + 1
        }
        all_results
      },
      .progress = progress_option_domain
    )
  }

  if (!missing(manual_query)) {
    # Check if any element is of type character to report issues with the fetching to the user.
    if (any(purrr::map_lgl(query_result[[1]], is.character))) {
      issue <- purrr::keep(query_result[[1]], is.character)[[1]]

      message("\nThere was an issue fetching from InterPro:\n", issue)
      return(NULL)
    } else {
      # return raw data for manual queries
      return(query_result)
    }
  }

  # check if any querys were not completed and give a message about them
  if(exists("query_result") && any(purrr::map_lgl(query_result, .f = ~ {length(.x) == 0}))){
    problematic_ids <- uniprot_ids[purrr::map_lgl(query_result, .f = ~ {length(.x) == 0})]

    message("\nThe following IDs do not exist in the InterPro database:\n", paste0(problematic_ids, collapse = "\n"))
  }

  # If any element is character stop the function and report the issue
  if (exists("query_result") && any(purrr::map_lgl(query_result, .f = ~ {
    if (length(.x) == 0){
      FALSE
    } else {
      purrr::map_lgl(.x, is.character)
    }
  }))) {
    problematic_ids <- uniprot_ids[purrr::map_lgl(query_result, .f = ~ {
      if (length(.x) == 0){
        FALSE
      } else {
        any(purrr::map_lgl(.x, is.character))
      }
    })]
    issue <- purrr::keep(query_result, .p = ~ {
      if (length(.x) == 0){
        FALSE
      } else {
        purrr::map_lgl(.x, is.character)
      }
    }) %>%
      unlist()

    message("\nThere was an issue fetching from InterPro:\n", paste0(problematic_ids, ": ", issue, "\n"))
    return(NULL)
  }

  # if residue level information should be returned, fetch it here
  if (return_residue_info) {
    query_residue <- purrr::map_dfr(
      .x = uniprot_ids,
      .f = ~ {
        # Initialize the next page URL with the base query URL
        next_url <- paste0("https://www.ebi.ac.uk/interpro/api/protein/uniprot/", .x, "?residues")

        query <- protti:::try_query(next_url,
          type = "application/json",
          max_tries = max_tries,
          timeout = timeout
        )

        if (is(query, "character")) {
          result <- data.frame(
            accession = .x,
            issue = query
          )
        } else {
          if (length(query) == 0) {
            result <- data.frame(accession = .x)
          } else {
            # map for each database accession
            result <- purrr::map_dfr(
              .x = query,
              .f = ~ {
                # map for each location
                location <- purrr::map_dfr(
                  .x = .x[["locations"]],
                  .f = ~ {
                    # map for each fragment
                    fragment <- purrr::map_dfr(
                      .x = .x[["fragments"]],
                      .f = ~ {
                        # lowest level
                        data.frame(
                          start = .x[["start"]],
                          end = .x[["end"]],
                          residues = .x[["residues"]]
                        )
                      }
                    )

                    fragment %>%
                      dplyr::mutate(fragment_description = .x[["description"]])
                  }
                )
                location %>%
                  dplyr::mutate(
                    source_database = .x[["source_database"]],
                    source_accession = .x[["accession"]],
                    source_name = .x[["name"]]
                  )
              }
            ) %>%
              mutate(accession = .x)
          }
        }

        result
      },
      .progress = progress_option_residue
    )

    if ("issue" %in% colnames(query_residue)) {
      problematic_ids <- query_residue %>%
        dplyr::filter(!is.na(issue))

      message("\nThere was an issue fetching from InterPro:\n", paste0(problematic_ids$accession, ": ", problematic_ids$issue, "\n"))
      return(NULL)
    }

    return(query_residue)
  }

  # start with unnesting and flattening the data
  flattened_query_result <- purrr::map(
    .x = query_result,
    .f = ~ {
      # page level
      purrr::map(
        .x = .x,
        .f = ~ {
          .x[["results"]]
        }
      )
    }
  ) %>%
    purrr::flatten() %>%
    purrr::flatten()

  # unnest domain level data
  domain_result <- purrr::map_dfr(
    .x = flattened_query_result,
    .f = ~ {
      go_terms <- purrr::map_dfr(
        .x = .x[["metadata"]][["go_terms"]],
        .f = ~ {
          data.frame(
            go_id = .x[["identifier"]],
            go_name = .x[["name"]],
            go_code = .x[["category"]][["code"]],
            go_type = .x[["category"]][["name"]]
          )
        }
      )

      metadata <- data.frame(
        identifier = .x[["metadata"]][["accession"]],
        identifier_name = .x[["metadata"]][["name"]],
        identifier_source_database = .x[["metadata"]][["source_database"]],
        identifier_type = .x[["metadata"]][["type"]]
      )

      if (nrow(go_terms) != 0) {
        metadata <- metadata %>%
          bind_cols(go_terms)
      }

      # Protein level
      protein_level <- purrr::map_dfr(
        .x = .x[["proteins"]],
        .f = ~ {
          # Entry protein location level
          entry_protein_locations_level <- purrr::map_dfr(
            .x = .x[["entry_protein_locations"]],
            .f = ~ {
              # Fragment Level
              protein_info <- purrr::map_dfr(
                .x = .x[["fragments"]],
                .f = ~ {
                  data.frame(
                    start = .x[["start"]],
                    end = .x[["end"]],
                    dc_status = .x[["dc-status"]],
                    representative = .x[["representative"]]
                  )
                }
              ) %>%
                dplyr::mutate(
                  model = .x[["model"]],
                  score = .x[["score"]]
                )
            }
          ) %>%
            dplyr::mutate(accession = stringr::str_to_upper(.x[["accession"]]))
        }
      )

      tidyr::crossing(metadata, protein_level)
    }
  )

  return(domain_result)
}
