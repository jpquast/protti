#' Extract metal-bind protein information from UniProt 
#'
#' Information of metal binding proteins is extracted from UniProt data retreived with \code{fetch_uniprot}. ChEBI ID's, potential
#' sub-ID's for metal cations, binding site locations in the protein and sub-ID evidence level (based on metal presence as cofactor) are extracted.
#'
#' @param data A data frame containing at least the input columns.
#' @param protein_id The name of the column containing protein identifiers.
#' @param feature_metal_binding The name of the column containing feature metal binding information from UniProt.
#' @param chebi_cofactor The name of the column containing ChEBI cofactor information from UniProt.
#' @param chebi_catalytic_activity The name of the column containing ChEBI catalytic activity information from UniProt.
#' @param chebi_data Optinal, a data frame that can be manually obtained with \code{fetch_chebi()}. If not provided it will be fetched within the function. 
#' If the funciton is run many times it is recomended to provide the data frame to save time.
#' @param chebi_relation_data Optinal, a data frame that can be manually obtained with \code{fetch_chebi(relation = TRUE)}. If not provided it will be fetched within the function. 
#' If the funciton is run many times it is recomended to provide the data frame to save time.
#' 
#' @return A data frame containing information on protein metal binding state. It contains the following types of columns (the naming might vary based on the input): 
#' \itemize{
#' \item{\code{protein_id}: }{UniProt protein identifier.}
#' \item{\code{source}: }{The source of the information, can be either \code{feature_metal_binding}, \code{chebi_cofactor} or \code{chebi_catalytic_activity}.}
#' \item{\code{ids}: }{ChEBI ID assigned to protein and binding site based on \code{metal_type} column name. These are general ID's that have sub-ID's. Thus, they generally describe the type of metal ion bound to the protein.}
#' \item{\code{metal_position}: }{Amino acid position within the protein that is involved in metal binding.}
#' \item{\code{metal_type}: }{Metal name extracted from \code{feature_metal_binding} information. This is the name that is used as a search pattern in order to assign a ChEBI ID with the \code{split_metal_name} helper function within this function.}
#' \item{\code{sub_ids}: }{ChEBI ID that is a sub-ID (incoming) of the ID in the \code{ids} column. Thus, they more specifically describe the potential nature of the metal ion.}
#' \item{\code{main_id_name}: }{Official ChEBI name associated with the ID in the \code{ids} column.}
#' \item{\code{multi_evidence}: }{If there is overlapping information in \code{feature_metal_binding} and \code{chebi_cofactor}, only \code{feature_metal_binding} is retained and multi_evidence is TRUE.}
#' \item{\code{sub_id_name}: }{Official ChEBI name associated with the ID in the \code{sub_ids} column.}
#' } 
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' extract_metal_binders(
#' data
#' )
#' }
extract_metal_binders <-
  function(data,
           protein_id = id,
           feature_metal_binding = feature_metal_binding,
           chebi_cofactor = chebi_cofactor,
           chebi_catalytic_activity = chebi_catalytic_activity,
           chebi_data = NULL,
           chebi_relation_data = NULL) {
    # Download chebi database if not provided
    if (is.null(chebi_data)) {
      chebi <- fetch_chebi()
    } else {
      chebi <- chebi_data
    }
    # Download chebi relation dataset if not provided
    if (is.null(chebi_relation_data)) {
      chebi_relation <- fetch_chebi(relation = TRUE)
    } else {
      chebi_relation <- chebi_relation_data
    }
    # Subset chebi data to sub id's of "metal cation" and "iron sulfur cluster"
    metal_cation <-
      find_all_subs(chebi_relation, id = "25213", type = "is_a")[[1]]
    iron_sulfur_cluster <-
      find_all_subs(chebi_relation, id = "30408", type = "is_a")[[1]]
    
    chebi_metal <- chebi %>%
      dplyr::filter(id %in% metal_cation)
    
    chebi_iron_sulfur <- chebi %>%
      dplyr::filter(id %in% iron_sulfur_cluster) %>%
      dplyr::mutate(name = stringr::str_replace_all(
        .data$name,
        pattern = stringr::regex("iron", ignore_case = TRUE),
        replacement = ""
      ))
    
    chebi_metal <- chebi_metal %>%
      dplyr::bind_rows(chebi_iron_sulfur)
    
    # Extract chebi standard id name pairs
    chebi_names <- chebi %>%
      dplyr::filter(.data$type_name == "STANDARD") %>%
      dplyr::distinct(.data$id, .data$name) %>%
      dplyr::mutate(id = as.character(.data$id))
    
    # Extract metal binding information
    result <- data %>%
      dplyr::distinct({{protein_id}}, {{feature_metal_binding}}, {{chebi_cofactor}}, {{chebi_catalytic_activity}}) %>%
      tidyr::pivot_longer(-{{protein_id}}, names_to = "source", values_to = "ids") %>%
      dplyr::filter(.data$ids != "") %>%
      dplyr::mutate(feature_metal_binding_sep = stringr::str_extract_all(.data$ids, pattern = "METAL.+?(?=METAL)|METAL.+$")) %>%
      dplyr::mutate(
        ids = ifelse(
          as.character(.data$feature_metal_binding_sep) == "character(0)",
          .data$ids,
          .data$feature_metal_binding_sep
        )
      ) %>%
      dplyr::select(-.data$feature_metal_binding_sep) %>%
      tidyr::unnest(.data$ids) %>%
      dplyr::mutate(cofactor = stringr::str_extract_all(.data$ids, pattern = "(?<=\\[CHEBI:).+?(?=\\])")) %>%
      dplyr::mutate(
        ids = ifelse(
          as.character(.data$cofactor) == "character(0)",
          .data$ids,
          .data$cofactor
        )
      ) %>%
      dplyr::select(-.data$cofactor) %>%
      tidyr::unnest(.data$ids) %>%
      dplyr::mutate(metal_position = ifelse(source == "feature_metal_binding",  as.numeric(stringr::str_extract(.data$ids, pattern = "\\d+")), NA)) %>%
      dplyr::mutate(metal_type = stringr::str_extract(.data$ids, pattern = "(?<=note=\\\").+?(?=\\\"|\\;)" )) %>%
      dplyr::mutate(metal_type = stringr::str_remove_all(.data$metal_type, pattern = "\\s\\d$")) %>%
      dplyr::mutate(pattern = split_metal_name(.data$metal_type)) %>%
      dplyr::mutate(chebi_ids = find_chebis(chebi_metal, .data$pattern)) %>%
      dplyr::mutate(
        ids = ifelse(
          as.character(.data$chebi_ids) == "character(0)",
          .data$ids,
          .data$chebi_ids
        )
      ) %>%
      dplyr::select(-c(.data$pattern, .data$chebi_ids)) %>%
      tidyr::unnest(.data$ids) %>%
      dplyr::mutate(is_metal = .data$ids %in% as.character(metal_cation) | .data$ids %in% as.character(iron_sulfur_cluster)) %>%
      dplyr::filter(.data$is_metal == TRUE) %>%
      dplyr::mutate(sub_ids = find_all_subs(chebi_relation, id = .data$ids, type = "is_a")) %>%
      dplyr::mutate(null = purrr::map(.data$sub_ids, ~is.null(.))) %>%
      dplyr::mutate(sub_ids = ifelse((source == rlang::as_name(rlang::enquo(chebi_cofactor)) | source == rlang::as_name(rlang::enquo(chebi_catalytic_activity))) & .data$null == TRUE, as.integer(.data$ids), .data$sub_ids)) %>%
      dplyr::select(-.data$null) %>%
      # unnest removes all NULL values in sub_ids. This might cause some ID's without sub ID's to get lost.
      tidyr::unnest(.data$sub_ids) %>%
      dplyr::left_join(chebi_names, by = c("ids" = "id")) %>%
      dplyr::rename(main_id_name = .data$name) %>%
      dplyr::distinct() %>%
      dplyr::group_by(.data$protein_name, .data$metal_type) %>%
      dplyr::mutate(n_position = dplyr::n_distinct(.data$metal_position, na.rm = TRUE)) %>%
      dplyr::group_by(.data$protein_name, .data$sub_ids) %>%
      dplyr::mutate(n_ids = dplyr::n()) %>%
      dplyr::mutate(multi_evidence = ifelse(.data$n_ids/.data$n_position > 1 & .data$n_position != 0, TRUE, FALSE)) %>%
      dplyr::filter(!(.data$n_position == 0 & .data$n_ids > 1)) %>%
      dplyr::select(-c(.data$is_metal, .data$n_position, .data$n_ids)) %>%
      dplyr::mutate(sub_ids = as.character(.data$sub_ids)) %>%
      dplyr::left_join(chebi_names, by = c("sub_ids" = "id")) %>%
      dplyr::rename(sub_id_name = .data$name)
    
    result
  } 