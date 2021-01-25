#' Analyse protein interaction network for significant hits
#'
#' The STRING database provides a resource for known and predicted protein-protein interactions. The type of
#' interactions include direct (physical) and indirect (functional) interactions. Through the R package
#' \code{STRINGdb} this resource if provided to R users. This function provides a convenient wrapper for
#' \code{STRINGdb} functions that allow an easy use within the protti pipeline.
#'
#' @param data A data frame that contains significantly changing proteins (STRINGdb is only able to plot
#' 400 proteins at a time so do not provide more for network plots). Information about treatment binding can be provided and will be displayed as colorful halos around
#' the proteins in the network.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param string_id The name of the column containing STRING database identifiers. These can be obtained from UniProt.
#' @param organism_id Numeric organism ID (NCBI taxon-ID). This can be obtained from \href{https://string-db.org/cgi/input?sessionId=bpvps5GS2As6&input_page_show_search=on}{here}.
#' H. sapiens: 9606, S. cerevisiae: 4932, E. coli: 511145.
#' @param score_threshold The interaction score based on \href{https://string-db.org/cgi/info?sessionId=bBP5N4cIf0PA&footer_active_subpage=scores}{STRING}
#' has to be between 0 and 1000. A score closer to 1000 is related to a higher confidence for the interaction. The default value is 900.
#' @param binds_treatment The name of the column containing a logical indicating if the corresponding
#' protein binds to the treatment. This information can be obtained from different databases, e.g UniProt.
#' @param halo_color Optional, A character vector with a color hex-code. This is the color of the halo of proteins that bind the treatment.
#' @param plot A logical indicating whether the result should be plotted or returned as a table.
#'
#' @return A network plot displaying interactions of the provided proteins. If \code{binds_treatment} was provided halos around
#' the proteins show which proteins interact with the treatment. If \code{plot = FALSE} a table with interaction information is
#' returned.
#'
#' @import STRINGdb
#' @importFrom dplyr distinct pull mutate filter rename
#' @importFrom rlang .data ensym !! as_name enquo
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' \dontrun{
#' network_analysis(
#'   data,
#'   protein_id = pg_protein_accessions,
#'   string_id = database_string,
#'   organism_id = 511145,
#'   binds_treatment = is_known,
#'   plot = TRUE
#' )
#' }
network_analysis <- function(data, protein_id, string_id, organism_id, score_threshold = 900, binds_treatment = NULL, halo_color = NULL, plot = TRUE) {
  data <- data %>%
    dplyr::distinct({{ protein_id }}, {{ string_id }}, {{ binds_treatment }})

  if (length(unique(dplyr::pull(data, !!ensym(protein_id)))) != nrow(data)) stop("Please provide unique annotations for each protein! The number of proteins does not match the number of rows in your data.")

  string_db <- STRINGdb$new(
    version = "11",
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
    if (missing(halo_color)) {
      coloring <- input %>%
        dplyr::filter({{ binds_treatment }}) %>%
        dplyr::mutate(color = "#5680C1")
    } else {
      coloring <- input %>%
        dplyr::filter({{ binds_treatment }}) %>%
        dplyr::mutate(color = halo_color)
    }
    payload_id <- string_db$post_payload(coloring$database_string,
      colors = coloring$color
    )
  }
  if (plot == TRUE) {
    if (length(unique(dplyr::pull(data, !!ensym(protein_id)))) > 400) stop("Please only provide the top 400 significant proteins for plots! String cannot plot more at once.")
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
}
