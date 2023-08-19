#' Predict protein domains of AlphaFold predictions
#'
#' Uses the predicted aligned error (PAE) of AlphaFold predictions to find possible protein domains.
#' A graph-based community clustering algorithm (Leiden clustering) is used on the predicted error
#' (distance) between residues of a protein in order to infer pseudo-rigid groups in the protein. This is
#' for example useful in order to know which parts of protein predictions are likely in a fixed relative
#' position towards each other and which might have varying distances.
#' This function is based on python code written by Tristan Croll. The original code can be found on his
#' \href{https://github.com/tristanic/pae_to_domains}{GitHub page}.
#'
#' @param pae_list a list of proteins that contains aligned errors for their AlphaFold predictions.
#' This list can be retrieved with the `fetch_alphafold_aligned_error()` function. It should contain a
#' column containing the scored residue (`scored_residue`), the aligned residue (`aligned_residue`) and
#' the predicted aligned error (`error`).
#' @param pae_power a numeric value, each edge in the graph will be weighted proportional to (`1 / pae^pae_power`).
#' Default is `1`.
#' @param pae_cutoff a numeric value, graph edges will only be created for residue pairs with `pae < pae_cutoff`.
#' Default is `5`.
#' @param graph_resolution a numeric value that regulates how aggressive the clustering algorithm is. Smaller values
#' lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.
#' Higher values lead to stricter (i.e. smaller) clusters. The value is provided to the Leiden clustering algorithm
#' of the `igraph` package as `graph_resolution / 100`. Default is `1`.
#' @param return_data_frame a logical value; if `TRUE` a data frame instead of a list
#' is returned. It is recommended to only use this if information for few proteins is retrieved.
#' Default is `FALSE`.
#' @param show_progress a logical value that specifies if a progress bar will be shown. Default
#' is `TRUE`.
#'
#' @return A list of the provided proteins that contains domain assignments for each residue. If `return_data_frame` is
#' `TRUE`, a data frame with this information is returned instead. The data frame contains the
#' following columns:
#' \itemize{
#' \item{residue: }{The protein residue number.}
#' \item{domain: }{A numeric value representing a distinct predicted domain in the protein.}
#' \item{accession: }{The UniProt protein identifier.}
#' }
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # Fetch aligned errors
#' aligned_error <- fetch_alphafold_aligned_error(
#'   uniprot_ids = c("F4HVG8", "O15552"),
#'   error_cutoff = 4
#' )
#'
#' # Predict protein domains
#' af_domains <- predict_alphafold_domain(
#'   pae_list = aligned_error,
#'   return_data_frame = TRUE
#' )
#'
#' head(af_domains, n = 10)
#' }
predict_alphafold_domain <- function(pae_list,
                                     pae_power = 1,
                                     pae_cutoff = 5,
                                     graph_resolution = 1,
                                     return_data_frame = FALSE,
                                     show_progress = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
    return(invisible(NULL))
  }

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(pae_list),
      format = "Predicting protein domains [:bar] :current/:total (:percent) :eta"
    )
  }

  # Calculate domains
  result <- pae_list %>%
    purrr::map(.f = ~ {
      if (show_progress == TRUE) {
        pb$tick()
      }
      # Create aligned error matrix
      aligned_error_matrix <- .x %>%
        dplyr::select(c("scored_residue", "aligned_residue", "error")) %>%
        # prevent the creation of Inf weights. Convert all 0 values to a low value instead.
        dplyr::mutate(error = ifelse(error == 0, 0.001, .data$error)) %>%
        tidyr::pivot_wider(
          names_from = "scored_residue",
          values_from = "error"
        ) %>%
        dplyr::select(-c("aligned_residue")) %>%
        as.matrix()

      # Calculate all weights
      weights <- 1 / aligned_error_matrix^pae_power

      # Create index list for edges
      indexes <- which(aligned_error_matrix < pae_cutoff, arr.ind = TRUE)
      indexes_sorted <- as.numeric(unlist(stringr::str_split(paste0(indexes[, 1], ",", indexes[, 2], collapse = ","), pattern = ",")))

      # Adjust weights to include only those falling within the cutoff
      sel_weight <- weights[indexes]

      # Create graph
      g <- igraph::make_empty_graph(directed = FALSE)

      # Create vertices based on numer of residues
      g <- igraph::add_vertices(g, nv = dim(aligned_error_matrix)[1])

      # Add edges to graph
      g <- igraph::add_edges(g, indexes_sorted)

      # Assign weights as "weight" attribute to edges
      g <- igraph::set_edge_attr(g, "weight", value = sel_weight)

      # Cluster with "Leiden"
      vc <- igraph::cluster_leiden(g, weights = NULL, resolution_parameter = graph_resolution / 100)

      tibble::tibble(
        residue = 1:dim(aligned_error_matrix)[1],
        domain = as.numeric(igraph::membership(vc))
      )
    })

  if (return_data_frame == FALSE) {
    return(result)
  } else {
    result_df <- purrr::map2_dfr(
      .x = result,
      .y = names(result),
      .f = ~ {
        .x %>%
          dplyr::mutate(accession = .y)
      }
    )
    return(result_df)
  }
}
