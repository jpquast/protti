#' Analyse protein interaction network for significant hits
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `analyse_functional_network()`.
#'
#' @return A network plot displaying interactions of the provided proteins. If
#' \code{binds_treatment} was provided halos around the proteins show which proteins interact with
#' the treatment. If \code{plot = FALSE} a data frame with interaction information is returned.
#' @keywords internal
#' @export
network_analysis <-
  function(...) {
    # This function has been renamed and is therefore deprecated.
    lifecycle::deprecate_warn("0.2.0",
      "network_analysis()",
      "analyse_functional_network()",
      details = "This function has been renamed."
    )

    analyse_functional_network(...)
  }
#' Analyse protein interaction network for significant hits
#'
#' The STRING database provides a resource for known and predicted protein-protein interactions.
#' The type of interactions include direct (physical) and indirect (functional) interactions.
#' Through the R package \code{STRINGdb} this resource if provided to R users. This function
#' provides a convenient wrapper for \code{STRINGdb} functions that allow an easy use within the
#' protti pipeline.
#'
#' @param data a data frame that contains significantly changing proteins (STRINGdb is only able
#' to plot 400 proteins at a time so do not provide more for network plots). Information about
#' treatment binding can be provided and will be displayed as colorful halos around the proteins
#' in the network.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' accession numbers.
#' @param string_id a character column in the \code{data} data frame that contains STRING database
#' identifiers. These can be obtained from UniProt.
#' @param organism_id a numeric value specifying an organism ID (NCBI taxon-ID). This can be
#' obtained from
#' \href{https://string-db.org/cgi/input?sessionId=bpvps5GS2As6&input_page_show_search=on}{here}.
#' H. sapiens: 9606, S. cerevisiae: 4932, E. coli: 511145.
#' @param version a character value that specifies the version of STRINGdb to be used.
#' Default is 12.0.
#' @param score_threshold a numeric value specifying the interaction score that based on
#' \href{https://string-db.org/cgi/info?sessionId=bBP5N4cIf0PA&footer_active_subpage=scores}{STRING}
#' has to be between 0 and 1000. A score closer to 1000 is related to a higher confidence for the
#' interaction. The default value is 900.
#' @param binds_treatment a logical column in the \code{data} data frame that indicates if the
#' corresponding protein binds to the treatment. This information can be obtained from different
#' databases, e.g UniProt.
#' @param halo_color optional, character value with a color hex-code. This is the color of the
#' halo of proteins that bind the treatment.
#' @param plot a logical that indicates whether the result should be plotted or returned as a table.
#'
#' @return A network plot displaying interactions of the provided proteins. If
#' \code{binds_treatment} was provided halos around the proteins show which proteins interact with
#' the treatment. If \code{plot = FALSE} a data frame with interaction information is returned.
#'
#' @importFrom dplyr distinct pull mutate filter rename
#' @importFrom rlang .data ensym !! as_name enquo
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' data <- data.frame(
#'   uniprot_id = c(
#'     "P0A7R1",
#'     "P02359",
#'     "P60624",
#'     "P0A7M2",
#'     "P0A7X3",
#'     "P0AGD3"
#'   ),
#'   xref_string = c(
#'     "511145.b4203;",
#'     "511145.b3341;",
#'     "511145.b3309;",
#'     "511145.b3637;",
#'     "511145.b3230;",
#'     "511145.b1656;"
#'   ),
#'   is_known = c(
#'     TRUE,
#'     TRUE,
#'     TRUE,
#'     TRUE,
#'     TRUE,
#'     FALSE
#'   )
#' )
#'
#' # Perform network analysis
#' network <- analyse_functional_network(
#'   data,
#'   protein_id = uniprot_id,
#'   string_id = xref_string,
#'   organism_id = 511145,
#'   binds_treatment = is_known,
#'   plot = TRUE
#' )
#'
#' network
#' }
analyse_functional_network <- function(data,
                                       protein_id,
                                       string_id,
                                       organism_id,
                                       version = "12.0",
                                       score_threshold = 900,
                                       binds_treatment = NULL,
                                       halo_color = NULL,
                                       plot = TRUE) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    message(strwrap("Package \"STRINGdb\" is needed for this function to work. Please install it.",
      prefix = "\n", initial = ""
    ), call. = FALSE)
    return(invisible(NULL))
  }

  # Ensure data frame is not empty and columns exist
  if (nrow(data) == 0) {
    stop("The input data frame is empty.")
  }

  required_columns <- c(ensym(protein_id), ensym(string_id))
  missing_columns <- required_columns[!required_columns %in% colnames(data)]
  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from the input data frame: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (plot && nrow(data) > 400) {
    stop("Please only provide the top 400 significant proteins for plots! STRING cannot plot more at once.")
  }


  data <- data %>%
    dplyr::distinct({{ protein_id }}, {{ string_id }}, {{ binds_treatment }})
  if (length(unique(dplyr::pull(data, !!ensym(protein_id)))) != nrow(data)) {
    stop(strwrap("Please provide unique annotations for each protein! The number of proteins
    does not match the number of rows in your data.", prefix = "\n", initial = ""))
  }

  STRINGdb <- get("STRINGdb", envir = loadNamespace("STRINGdb"))

  string_db <- STRINGdb$new(
    version = version,
    species = organism_id, # Check on String database to get the right code (E.coli K12: 511145)
    score_threshold = score_threshold, # Cutoff score to consider something an interaction
    input_directory = ""
  )

  input <- data %>%
    dplyr::mutate({{ string_id }} := stringr::str_extract({{ string_id }}, pattern = ".+[^;]")) %>%
    tidyr::drop_na({{ string_id }})

  string_ids <- dplyr::pull(input, !!ensym(string_id))

  payload_id <- NULL

  if (!missing(binds_treatment)) {
    if (is.null(halo_color)) {
      halo_color <- "#5680C1"
    }

    coloring <- input %>%
      dplyr::filter({{ binds_treatment }}) %>%
      dplyr::mutate(color = halo_color)

    payload_id <- string_db$post_payload(dplyr::pull(coloring, {{ string_id }}),
      colors = coloring$color
    )
  }


  tryCatch(
    {
      if (plot) {
        string_db$plot_network(string_ids, payload_id = payload_id)
      } else {
        mapping <- input %>%
          dplyr::distinct({{ protein_id }}, {{ string_id }})

        interactions <- string_db$get_interactions(string_ids) %>%
          dplyr::left_join(mapping, by = c("from" = rlang::as_name(rlang::enquo(string_id)))) %>%
          dplyr::rename(from_protein = {{ protein_id }}) %>%
          dplyr::left_join(mapping, by = c("to" = rlang::as_name(rlang::enquo(string_id)))) %>%
          dplyr::rename(to_protein = {{ protein_id }}) %>%
          dplyr::distinct()
        return(interactions)
      }
    },
    error = function(e) {
      message("An error occurred during the interaction network analysis: ", e$message)
      return(invisible(NULL))
    }
  )
}
