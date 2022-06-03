#' Extract metal-bind protein information from UniProt
#'
#' Information of metal binding proteins is extracted from UniProt data retrieved with
#' \code{fetch_uniprot}. ChEBI IDs, potential sub-IDs for metal cations, binding site locations
#' in the protein and sub-ID evidence level (based on metal presence as cofactor) are extracted.
#'
#' @param data a data frame containing at least the input columns.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' identifiers.
#' @param feature_metal_binding a character column in the \code{data} data frame that contains the
#' feature metal binding information from UniProt.
#' @param chebi_cofactor a character column in the \code{data} data frame that contains the ChEBI
#' cofactor information from UniProt.
#' @param chebi_catalytic_activity a character column in the \code{data} data frame that contains
#' the ChEBI catalytic activity information from UniProt.
#' @param comment_cofactor a character cholumn in the \code{data} data frame that contains
#' the comment cofactor information from UniProt.
#' @param go_molecular_function a character column in the \code{data} data frame that contains 
#' the gene ontology "molecular function" column information from UniProt.
#' @param chebi_data optional, a data frame that can be manually obtained with \code{fetch_chebi()}.
#' If not provided it will be fetched within the function. If the function is run many times it
#' is recommended to provide the data frame to save time.
#' @param chebi_relation_data optional, a data frame that can be manually obtained with
#' \code{fetch_chebi(relation = TRUE)}. If not provided it will be fetched within the function.
#' If the function is run many times it is recommended to provide the data frame to save time.
#'
#' @return A data frame containing information on protein metal binding state. It contains the
#' following types of columns (the naming might vary based on the input):
#' \itemize{
#' \item{\code{protein_id}: }{UniProt protein identifier.}
#' \item{\code{metal_type}: }{Metal name extracted from \code{feature_metal_binding} information.
#' This is the name that is used as a search pattern in order to assign a ChEBI ID with the
#' \code{split_metal_name} helper function within this function.}
#' \item{\code{metal_position}: }{Amino acid position within the protein that is involved in metal
#' binding.}
#' \item{\code{metal_id}: }{This is a unique identifier for the metal of a protein. If the protein
#' binds only one metal this will be 1. If the protein e.g. binds 2 magnesium ions, one will have 
#' the ID 1, while the other has the ID 2.}
#' \item{\code{ids}: }{ChEBI ID assigned to protein and binding site based on \code{metal_type}
#' column name. These are general IDs that have sub-IDs. Thus, they generally describe the type of
#' metal ion bound to the protein.}
#' \item{\code{main_id_name}: }{Official ChEBI name associated with the ID in the \code{ids} column.}
#' \item{\code{sub_ids}: }{ChEBI ID that is a sub-ID (incoming) of the ID in the \code{ids} column.
#' Thus, they more specifically describe the potential nature of the metal ion.}
#' \item{\code{sub_id_name}: }{Official ChEBI name associated with the ID in the \code{sub_ids}
#' column.}
#' \item{\code{source}: }{The source of the information, can be either \code{feature_metal_binding},
#' \code{chebi_cofactor}, \code{comment_cofactor} or \code{chebi_catalytic_activity}.}
#' \item{\code{evidence}: }{Evidence level of this annotation. This is usually an evidence & conclusion 
#' ontology ID (ECO ID). Sometimes a publication ID is included.}
#' \item{\code{note}: }{A note specifying the annotation further. Can for example specify how many
#' metals are bound per protein.}
#' }
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import stringr
#' @importFrom stats na.omit
#' @importFrom stringi stri_trans_totitle stri_opts_brkiter
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' data <- fetch_uniprot(
#'   uniprot_ids = c("Q03640", "Q03778", "P22276"),
#'   columns = c(
#'     "feature(METAL BINDING)",
#'     "chebi(Cofactor)",
#'     "chebi(Catalytic activity)",
#'     "comment(COFACTOR)",
#'     "go(molecular function)"
#'   )
#' )
#'
#' # Extract metal binding information
#' metal_info <- extract_metal_binders(
#'   data = data,
#'   protein_id = id,
#'   feature_metal_binding = feature_metal_binding,
#'   chebi_cofactor = chebi_cofactor,
#'   chebi_catalytic_activity = chebi_catalytic_activity,
#'   comment_cofactor = comment_cofactor,
#'   go_molecular_function = go_molecular_function
#' )
#'
#' metal_info
#' }
extract_metal_binders <-
  function(data,
           protein_id = id,
           feature_metal_binding = feature_metal_binding,
           chebi_cofactor = chebi_cofactor,
           chebi_catalytic_activity = chebi_catalytic_activity,
           comment_cofactor = comment_cofactor,
           go_molecular_function = go_molecular_function,
           chebi_data = NULL,
           chebi_relation_data = NULL) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }
    if (!requireNamespace("stringi", quietly = TRUE)) {
      message("Package \"stringi\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }
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
    # Subset chebi data to sub IDs of "metal cation" and "iron sulfur cluster"
    metal_cation <-
      find_all_subs(chebi_relation, ids = "25213", types = "is_a")[[1]]
    inorganic_monovalent_cation <- find_all_subs(chebi_relation, ids = "60242", types = "is_a")[[1]]
    inorganic_metal_cation <- inorganic_monovalent_cation[!inorganic_monovalent_cation %in%
      c("29234", "29233", "29120", "28938", "24636", "15378")]
    # first is metal ions, second is monovalent inorganic cations. non-metals are excluded
    metal_cation <- c(metal_cation, inorganic_metal_cation)
    iron_sulfur_cluster <-
      find_all_subs(chebi_relation, ids = "30408", types = "is_a")[[1]]

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
    
    # Extract go_molecular_function information separately
    # metal_pattern_go might not be complete, if something is added also add to the renaming below.
    ion_pattern_go <- " ion | cation "
    special_name_pattern <- "ferric iron|ferrous iron|ferric|ferrous|cupric|cuprous" # also used below
    metal_pattern_go <- paste0(special_name_pattern, "|copper|iron|zinc|magnesium|calcium|cobalamin|heme|molybdopterin|nickel")
    iron_sulfur_pattern_go <- "4 iron, 4 sulfur cluster|2 iron, 2 sulfur cluster|3 iron, 4 sulfur cluster|iron-sulfur cluster binding"
    source_go <-  rlang::as_name(rlang::enquo(go_molecular_function)) # can't be used directly in mutate
    
    go_data <- data %>% 
      dplyr::distinct({{ protein_id }},
                      {{ go_molecular_function }}) %>% 
      tidyr::drop_na({{ go_molecular_function }}) %>% 
      dplyr::mutate({{ go_molecular_function }} := stringr::str_split({{ go_molecular_function }}, pattern = "; ")) %>% 
      tidyr::unnest({{ go_molecular_function }}) %>% 
      dplyr::filter(stringr::str_detect({{ go_molecular_function }}, pattern = paste0(ion_pattern_go, "|", metal_pattern_go, "|", iron_sulfur_pattern_go)) & 
                      !stringr::str_detect({{ go_molecular_function }}, pattern = "phosphate ion|chloride ion|organic cation")) %>% 
      dplyr::mutate(metal_type = stringr::str_extract({{ go_molecular_function }}, pattern = paste0("[:alpha:]+(?= ion)|", metal_pattern_go))) %>% 
      dplyr::mutate(iron_sulfur_type = stringr::str_extract({{ go_molecular_function }}, pattern = iron_sulfur_pattern_go)) %>% 
      dplyr::mutate(metal_type = stringr::str_split(ifelse(!is.na(.data$iron_sulfur_type), 
                                                  paste0(.data$metal_type, ";", .data$iron_sulfur_type), 
                                                  .data$metal_type), pattern = ";")) %>% 
      tidyr::unnest(.data$metal_type) %>%
      dplyr::select(-.data$iron_sulfur_type) %>% 
      dplyr::mutate(source = source_go,
                    metal_type = stringi::stri_trans_totitle(.data$metal_type,
                                                             opts_brkiter = stringi::stri_opts_brkiter(type = "sentence"))) %>%
      dplyr::mutate(metal_type = dplyr::case_when(.data$metal_type == "Metal" ~ "A metal cation",
                                                  # could also be "Iron-sulfur (Fe-S)" but feature_metal_binding also only contains Iron-sulfur
                                                  .data$metal_type == "Iron-sulfur cluster binding" ~ "Iron-sulfur", 
                                                  .data$metal_type == "4 iron, 4 sulfur cluster" ~ "Iron-sulfur (4Fe-4S)",
                                                  .data$metal_type == "2 iron, 2 sulfur cluster" ~ "Iron-sulfur (2Fe-2S)",
                                                  .data$metal_type == "3 iron, 4 sulfur cluster" ~ "Iron-sulfur (3Fe-4S)",
                                                  .data$metal_type == "Cobalamin" ~ "Cobalt",
                                                  .data$metal_type == "Ferric iron" | .data$metal_type == "Ferric" ~ "Ferric ion",
                                                  .data$metal_type == "Ferrous iron" | .data$metal_type == "Ferrous" ~ "Ferrous ion",
                                                  .data$metal_type == "Heme" ~ "Iron",
                                                  .data$metal_type == "Molybdopterin" | .data$metal_type == "Molybdate" ~ "Molybdenum",
                                                  .data$metal_type == "Cupric" ~ "Cupric ion",
                                                  .data$metal_type == "Cuprous" ~ "Cuprous ion",
                                                  .data$metal_type == "Gated" |
                                                  .data$metal_type == "Coupled" |
                                                  .data$metal_type == "Mechanosensitive" |
                                                  .data$metal_type == "Type" |
                                                  .data$metal_type == "Sensing" ~ NA_character_,
                                                  TRUE ~ .data$metal_type)) %>% 
      dplyr::group_by(.data$id) %>% 
      # remove A metal cation in case that there are more specific GO terms
      dplyr::filter(!(dplyr::n() > 1 & .data$metal_type == "A metal cation" & length(unique(stats::na.omit(.data$metal_type))) > 1))

    # Extract metal binding information
    result <- data %>%
      dplyr::distinct(
        {{ protein_id }},
        {{ feature_metal_binding }},
        {{ chebi_cofactor }},
        {{ chebi_catalytic_activity }},
        {{ comment_cofactor }}
      ) %>%
      tidyr::pivot_longer(-{{ protein_id }},
        names_to = "source",
        values_to = "ids"
      ) %>%
      dplyr::filter(.data$ids != "") %>%
      dplyr::mutate(feature_metal_binding_sep = stringr::str_extract_all(
        .data$ids,
        pattern = "METAL.+?(?=METAL)|METAL.+$"
      )) %>%
      dplyr::mutate(
        ids = ifelse(
          as.character(.data$feature_metal_binding_sep) == "character(0)",
          .data$ids,
          .data$feature_metal_binding_sep
        )
      ) %>%
      dplyr::select(-.data$feature_metal_binding_sep) %>%
      tidyr::unnest(.data$ids) %>%
      # split information if there are mutliple cofactor entries
      # this ensures that notes are still associated with the right entries
      dplyr::mutate(ids = ifelse(source == rlang::as_name(rlang::enquo(comment_cofactor)),
                                       stringr::str_extract_all(
                                         .data$ids,
                                         pattern = "(?<=COFACTOR:).+?(?=COFACTOR|$)"),
                                       .data$ids)) %>% 
      tidyr::unnest(.data$ids) %>% 
      # extract notes from comment_cofactor
      dplyr::mutate(note = ifelse(source == rlang::as_name(rlang::enquo(comment_cofactor)), 
                                  stringr::str_extract(
                                    .data$ids,
                                    pattern = "(?<=Note\\=).+?(?=;)"
                                  ),
                                  NA_character_)) %>% 
      tidyr::unnest(.data$note) %>% 
      dplyr::mutate(ids = ifelse(source == rlang::as_name(rlang::enquo(comment_cofactor)), 
                                                  stringr::str_extract_all(
        .data$ids,
        pattern = "(?<=Name\\=).+?(?=Name|Note|COFACTOR|$)"
      ),
      .data$ids)) %>% 
      tidyr::unnest(.data$ids) %>%
      # extract evidence from feature_metal_binding and comment_cofactor
      dplyr::mutate(evidence = stringr::str_extract(.data$ids, pattern = regex('(?<=evidence=).+(?=;|$)', ignore_case = TRUE))) %>% 
      dplyr::mutate(evidence = str_remove_all(.data$evidence, pattern = '"|;|\\{|\\}')) %>% 
      dplyr::mutate(cofactor = stringr::str_extract_all(
        .data$ids,
        pattern = "(?<=CHEBI:)\\d+"
      )) %>%
      dplyr::mutate(
        ids = ifelse(
          as.character(.data$cofactor) == "character(0)",
          .data$ids,
          .data$cofactor
        )
      ) %>%
      dplyr::select(-.data$cofactor) %>%
      tidyr::unnest(.data$ids) %>%
        # extract metal positions from feature_metal_binding
      dplyr::mutate(metal_position = ifelse(source == rlang::as_name(rlang::enquo(feature_metal_binding)),
        as.numeric(stringr::str_extract(.data$ids, pattern = "\\d+")),
        NA
      )) %>%
        # extract the metal type from the "note" section of feature_metal_binding
      dplyr::mutate(metal_type = stringr::str_extract(
        .data$ids,
        pattern = "(?<=note=\\\").+?(?=\\\"|\\;)"
      )) %>%
        # extract the metal ID from feature_metal_binding e.g. if there are mutliple metals per protein
      dplyr::mutate(metal_id = stringr::str_extract(.data$metal_type, pattern = "\\d$")) %>% 
      dplyr::mutate(metal_id = ifelse(is.na(.data$metal_id) & source == rlang::as_name(rlang::enquo(feature_metal_binding)), 1, .data$metal_id)) %>% 
       # clean up the metal type
      dplyr::mutate(metal_type = stringr::str_remove_all(.data$metal_type, pattern = "\\s\\d$")) %>%
      dplyr::bind_rows(go_data) %>% 
      # create a pattern of the metal name to look for matching ChEBI IDs
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
      dplyr::mutate(is_metal = .data$ids %in% as.character(metal_cation) |
        .data$ids %in% as.character(iron_sulfur_cluster) | is.na(.data$ids)) %>%
      dplyr::filter(.data$is_metal == TRUE) %>%
      # look for ChEBI sub IDs. Those are IDs that have the ID in the ids column as a parent
      dplyr::mutate(sub_ids = find_all_subs(chebi_relation, ids = .data$ids, types = "is_a")) %>%
      # find out which rows are NULL and dont have any sub IDs
      dplyr::mutate(null = purrr::map(.data$sub_ids, ~ is.null(.))) %>%
      # Do not include go terms and feature_metal_binding in transfer of IDs to subIDs to clean up. 
      # But keep specific metal_types for GO and feature to prevent them from being lost.
      dplyr::mutate(sub_ids = ifelse(((source == rlang::as_name(rlang::enquo(chebi_cofactor)) |
        source == rlang::as_name(rlang::enquo(chebi_catalytic_activity)) | 
          source == rlang::as_name(rlang::enquo(comment_cofactor))) &
        .data$null == TRUE) | 
          is.na(.data$ids) |
          (source == source_go &
          str_detect(.data$metal_type, pattern = " ion$")) |
            (source == rlang::as_name(rlang::enquo(feature_metal_binding)) &
               str_detect(.data$metal_type, pattern = "\\+\\)")), 
        as.integer(.data$ids), 
        .data$sub_ids)) %>%
      dplyr::select(-.data$null) %>%
      # unnest removes all NULL values in sub_ids. This might cause some IDs without sub IDs to get lost. 
      # some are saved with: c("Ferrous ion", "Ferric ion", "Cupric ion", "Cuprous ion") but there could be more
      tidyr::unnest(.data$sub_ids) %>%
      dplyr::left_join(chebi_names, by = c("ids" = "id")) %>%
      dplyr::rename(main_id_name = .data$name) %>%
      dplyr::distinct() %>%
      # combined chebi_cofactor and comment_cofactor information as much as possible
      dplyr::mutate(is_chebi_cofactor_comment_cofactor = ifelse(.data$source == rlang::as_name(rlang::enquo(chebi_cofactor)) | 
                                                                .data$source == rlang::as_name(rlang::enquo(comment_cofactor)),
                                                                TRUE,
                                                                FALSE)) %>% 
      dplyr::group_by({{ protein_id }}, .data$sub_ids, .data$is_chebi_cofactor_comment_cofactor) %>% 
      dplyr::mutate(n_id = dplyr::n()) %>% 
      dplyr::mutate(source = ifelse(.data$n_id == 2 & 
                                    .data$is_chebi_cofactor_comment_cofactor,
                                    paste0(rlang::as_name(rlang::enquo(chebi_cofactor)), " & ", rlang::as_name(rlang::enquo(comment_cofactor))),
                                    .data$source)) %>% 
      dplyr::group_by({{ protein_id }}, .data$sub_ids, .data$source) %>% 
      dplyr::mutate(evidence_all_na = all(is.na(.data$evidence))) %>% 
      dplyr::mutate(note_all_na = all(is.na(.data$note))) %>% 
      # make sure go_molecular_function is annotated for all rows with same sub IDS
      dplyr::group_by({{ protein_id }}, .data$sub_ids) %>% 
      dplyr::mutate({{ go_molecular_function }} := paste0(unique({{ go_molecular_function }}[!is.na({{ go_molecular_function }})]), collapse = ";")) %>% 
      dplyr::distinct() %>% 
      dplyr::filter(!(.data$source == paste0(rlang::as_name(rlang::enquo(chebi_cofactor)), " & ", rlang::as_name(rlang::enquo(comment_cofactor))) & 
                        dplyr::case_when(.data$evidence_all_na == TRUE & .data$note_all_na == TRUE ~ FALSE,
                                         .data$evidence_all_na == FALSE & .data$note_all_na == TRUE ~ is.na(.data$evidence),
                                         .data$evidence_all_na == TRUE & .data$note_all_na == FALSE ~ is.na(.data$note),
                                         .data$evidence_all_na == FALSE & .data$note_all_na == FALSE ~ is.na(.data$evidence) & is.na(.data$note)
                               ))) %>% 
      dplyr::group_by({{ protein_id }}, .data$sub_ids, {{ go_molecular_function }}) %>%
      # filter out GO annotations that are also part of other sources
      dplyr::filter(!(dplyr::n() > 1 & .data$source == source_go & !all(.data$source == source_go))) %>% 
      dplyr::group_by({{ protein_id }}, .data$sub_ids) %>% 
      # Combine information as much as possible
      dplyr::mutate(ids = paste0(unique(.data$ids[!is.na(.data$ids)]), collapse = ";"),
                    main_id_name = paste0(unique(.data$main_id_name[!is.na(.data$main_id_name)]), collapse = ";"),
                    source = paste0(unique(.data$source[!is.na(.data$source)]), collapse = " & "),
                    evidence = paste0(unique(.data$evidence[!is.na(.data$evidence)]), collapse = ";"),
                    note = paste0(unique(.data$note[!is.na(.data$note)]), collapse = ";")) %>% 
      dplyr::filter(!(dplyr::n() > 1 & is.na(.data$metal_type))) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # If a GO annotation is part of another source add it to source
      dplyr::mutate(source = ifelse({{ go_molecular_function }} != "" & .data$source != source_go, 
                                paste0(source, " & ", source_go), 
                                .data$source)) %>% 
      dplyr::select(-c(.data$n_id, .data$evidence_all_na, .data$note_all_na, .data$is_chebi_cofactor_comment_cofactor, .data$is_metal)) %>% 
      dplyr::mutate(sub_ids = as.character(.data$sub_ids)) %>%
      dplyr::left_join(chebi_names, by = c("sub_ids" = "id")) %>%
      dplyr::rename(sub_id_name = .data$name) %>% 
      dplyr::mutate(n_sources = stringr::str_count(.data$source, pattern = "&") + 1) %>% 
      # remove any rows that are part of a group with chebi IDs but that do not contain any chebi IDs themselves
      dplyr::mutate(is_iron_cation = stringr::str_detect(.data$main_id_name, pattern = "iron cation")) %>% # to not remove Iron IDs from FeS clusters
      dplyr::group_by(.data$id, .data$metal_type, .data$is_iron_cation) %>% 
      dplyr::mutate(any_chebi = any(stringr::str_detect(.data$source, pattern = "chebi"))) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(ifelse(.data$any_chebi,
                     (stringr::str_detect(.data$source, pattern = "chebi") |
                         # If GO contains specific iron and copper related patterns do not filter out.
                         # This is to prevent accidental removal of specific information about the bound metal
                         stringr::str_detect({{go_molecular_function}}, pattern = special_name_pattern)),
                     TRUE)) %>% 
      dplyr::select({{ protein_id }}, 
                    .data$metal_type, 
                    .data$metal_position, 
                    .data$metal_id, 
                    .data$ids, 
                    .data$main_id_name, 
                    .data$sub_ids, 
                    .data$sub_id_name,
                    {{ go_molecular_function }},
                    .data$source, 
                    .data$n_sources,
                    .data$evidence, 
                    .data$note) %>% 
      dplyr::arrange({{ protein_id }})

    result
  }