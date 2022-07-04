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
#' #metal_info <- extract_metal_binders(
#'  # data = data,
#'  # protein_id = id,
#'  # feature_metal_binding = feature_metal_binding,
#'  # chebi_cofactor = chebi_cofactor,
#'  # chebi_catalytic_activity = chebi_catalytic_activity,
#'  # comment_cofactor = comment_cofactor,
#'  # go_molecular_function = go_molecular_function
#' #)
#'
#' #metal_info
#' }
extract_metal_binders <-
  function(data_uniprot,
           data_quickgo,
           data_chebi = NULL,
           data_chebi_relation = NULL,
           data_eco = NULL,
           data_eco_relation = NULL,
           show_progress = TRUE) {
    # Check if required R packages are installed
    if (!requireNamespace("igraph", quietly = TRUE)) {
      message("Package \"igraph\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }
    if (!requireNamespace("stringi", quietly = TRUE)) {
      message("Package \"stringi\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }
    # Check if provided data has the right format
    # data_uniprot
    if(!("feature_metal_binding" %in% colnames(data_uniprot) &
         "comment_cofactor" %in% colnames(data_uniprot) & 
         "comment_catalytic_activity" %in% colnames(data_uniprot))){
      stop('Please include at least the columns "feature_metal_binding", "comment_cofactor" and "comment_catalytic_activity" in "data_uniprot"!')
    }
    # data_quickgo
    if(!("gene_product_db" %in% colnames(data_quickgo) &
         "gene_product_id" %in% colnames(data_quickgo) &
         "go_name" %in% colnames(data_quickgo) & 
         "go_term" %in% colnames(data_quickgo) &
         "go_aspect" %in% colnames(data_quickgo) &
         "eco_id" %in% colnames(data_quickgo) &
         "go_evidence_code" %in% colnames(data_quickgo) &
         "reference" %in% colnames(data_quickgo) &
         "with_from" %in% colnames(data_quickgo))){
      stop('Please include at least the columns "gene_product_db", "gene_product_id", "go_name", "go_term", "go_aspect", "go_evidence_code", "reference" and "with_from" in "data_quickgo"!')
    }
    # Download ChEBI database if not provided
    if (missing(data_chebi)) {
      if (show_progress == TRUE) {
        message("Downloading ChEBI information ... ", appendLF = FALSE)
        start_time <- Sys.time()
      }
      
      data_chebi <- fetch_chebi()
      # If database could not be retrieved let the function return NULL.
      if (is.null(data_chebi)){
        return(invisible(NULL))
      }
      
      if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
      }
    } 
    # Download chebi relation dataset if not provided
    if (missing(data_chebi_relation)) {
      if (show_progress == TRUE) {
        message("Downloading ChEBI relational information ... ", appendLF = FALSE)
        start_time <- Sys.time()
      }
      
      data_chebi_relation <- fetch_chebi(relation = TRUE)
      
      # If database could not be retrieved let the function return NULL.
      if (is.null(data_chebi_relation)){
        return(invisible(NULL))
      }
      
      if (show_progress == TRUE) {
        message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
      }
    } 
    
    # Prepare ChEBI dataframe for annotation
    data_chebi_filtered <- data_chebi %>% 
      dplyr::filter(.data$type_name == "STANDARD") %>% 
      dplyr::mutate(chebi_id = as.character(.data$id)) %>% 
      dplyr::distinct(chebi_id, definition, name, formula, mass, charge) 
    
    # Retrieve ECO annotation information if not provided
    if (missing(data_eco)) {
      if (show_progress == TRUE) {
        message("Downloading ECO annotation information from QuickGO ... ", appendLF = FALSE)
        start_time <- Sys.time()
      }
      
      eco_data <- fetch_eco(show_progress = FALSE)
      
      # If database could not be retrieved let the function return NULL.
      if (is.null(eco_chebi)){
        return(invisible(NULL))
      }
      
      if (show_progress == TRUE) {
        message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
      }
    }
    
    # Retrieve ECO relation information if not provided
    if (missing(data_eco_relation)) {
      if (show_progress == TRUE) {
        message("Downloading ECO relation information from QuickGO ... ", appendLF = FALSE)
        start_time <- Sys.time()
      }
      
      eco_relation <- fetch_eco(return_relation = TRUE, show_progress = FALSE)
      
      if (show_progress == TRUE) {
        message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
      }
    }
    
    # Create two vectors that contain all IDs related to either manual or automatic assertion
    manual_eco <- eco_relation %>% 
      find_all_subs(ids = c("ECO:0000352"), 
                    main_id = main_id, 
                    sub_id = sub_id, 
                    type = relation, 
                    accepted_types = "all") %>% 
      unlist()
    
    automatic_eco <- eco_relation %>% 
      find_all_subs(ids = c("ECO:0000501"), 
                    main_id = main_id, 
                    sub_id = sub_id, 
                    type = relation, 
                    accepted_types = "all") %>% 
      unlist()
    
    # Extract feature_metal_binding information from UniProt
    if (show_progress == TRUE) {
      message("Extract feature_metal_binding information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    fmb_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$id, .data$feature_metal_binding) %>% 
      tidyr::drop_na(.data$feature_metal_binding) %>% 
      # Extract names to be compatible with manual annotation
      dplyr::mutate(split_fmb = stringr::str_extract_all(
        .data$feature_metal_binding,
        pattern = "METAL.+?(?=METAL)|METAL.+$"
      )) %>%
      tidyr::unnest(.data$split_fmb) %>%
      dplyr::mutate(fmb_name = str_extract(.data$split_fmb, pattern = '(?<=/note=\\")[^\\";]+(?=[\\";])')) %>%
      dplyr::mutate(metal_identifier = stringr::str_extract(
        .data$fmb_name, 
        pattern = "(?<=[:space:])\\d{1,2}(?=[:space:]|$)")) %>% 
      dplyr::mutate(fmb_name = stringr::str_trim(
        stringr::str_replace_all(
          .data$fmb_name, 
          pattern = "[:space:]\\d{1,2}(?=[:space:]|$)", 
          replacement = ""))) %>% 
      # Deal with cases in which there is no metal name, then just put in "Divalent metal cation"
      dplyr::mutate(no_metal_name = is.na(.data$fmb_name),
        fmb_name = ifelse(is.na(.data$fmb_name),
                                      "Divalent metal cation",
                                      .data$fmb_name)) %>% 
      # Look for special identifiers and extract, but do not remove from the name
      dplyr::mutate(metal_identifier = ifelse(is.na(.data$metal_identifier), 
                                  str_extract(.data$fmb_name,
                                              pattern = "(?<=[:space:])A1$|A2$|Z1$|Z2$|Z3$|Z4$|A$|B$"),
                                  .data$metal_identifier)) %>% 
      # all values still NA get the identifier 1
      dplyr::mutate(metal_identifier = ifelse(is.na(.data$metal_identifier), 
                                              "1",
                                              .data$metal_identifier)) %>% 
      dplyr::mutate(metal_position = stringr::str_extract(.data$split_fmb, pattern = "(?<=METAL )\\d+")) %>% 
      dplyr::mutate(binding_mode = str_extract(
        str_extract(.data$split_fmb, pattern = '(?<=/note=\\")[^\\"]+(?=\\";)'),
        pattern = '(?<=; ).+(?=$)'
      )
      ) %>%
      dplyr::mutate(evidence = stringr::str_extract(.data$split_fmb, pattern = regex('(?<=evidence=).+(?=$)', ignore_case = TRUE))) %>% 
      dplyr::mutate(evidence = str_remove_all(.data$evidence, pattern = '"|;')) %>% 
      dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>% 
      tidyr::unnest(.data$evidence_split) %>% 
      tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>% 
      dplyr::distinct() %>% 
      # The the data that will be joined is a dataset provided with protti
      dplyr::left_join(fmb_annotation_uniprot, by = "fmb_name") %>% 
      # Label cases in which there was no metal name
      dplyr::mutate(comment = dplyr::case_when(.data$no_metal_name & !is.na(.data$comment) ~ paste0(.data$comment, ' There was no metal name therefore "Divalent metal cation" was used.'),
                                               .data$no_metal_name & is.na(.data$comment) ~ 'There was no metal name therefore "Divalent metal cation" was used.',
                                               TRUE ~ .data$comment)) %>% 
      dplyr::select(-c(.data$feature_metal_binding, .data$no_metal_name))
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
    
    # Check if there are any new terms not manually annotated in the annotation data frame
    if(any(!(fmb_uniprot$fmb_name %in% fmb_annotation_uniprot$fmb_name))){
      missing_names <- fmb_uniprot %>% 
        dplyr::distinct(.data$id, .data$fmb_name) %>% 
        dplyr::filter(!(.data$fmb_name %in% fmb_annotation_uniprot$fmb_name))
        
      warning(strwrap("The following names have been found in the feature_metal_binding column and have not yet been 
                      manually annotated in the reference data frame provided by protti. Please contact the package 
                      maintainer to let them be added.",
                      prefix = "\n", initial = ""),
              "\n",
              paste0(utils::capture.output(missing_names), collapse = "\n"))
    }
    
    # Extract comment_cofactor information from UniProt
    if (show_progress == TRUE) {
      message("Extract comment_cofactor information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    cofactor_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$id, .data$comment_cofactor) %>% 
      tidyr::drop_na(.data$comment_cofactor) %>% 
      dplyr::mutate(cofactor_split = stringr::str_extract_all(
                                   .data$comment_cofactor,
                                   pattern = "(?<=COFACTOR:).+?(?=COFACTOR|$)")) %>% 
      tidyr::unnest(.data$cofactor_split) %>% 
      # Extract notes
      dplyr::mutate(note = stringr::str_extract(
                                    .data$cofactor_split,
                                    pattern = "(?<=Note\\=).+?(?=;)"
                                  )) %>% 
      tidyr::unnest(.data$note) %>% 
      # Split names
      dplyr::mutate(name_split = stringr::str_extract_all(
                                   .data$cofactor_split,
                                   pattern = "(?<=Name\\=).+?(?=Name|Note|COFACTOR|$)"
                                 )) %>% 
      # prevent discarding entries without a "name" field e.g. some chlorophyl entries
      dplyr::mutate(
        name_split = ifelse(
          as.character(.data$name_split) == "character(0)",
          NA,
          .data$name_split
        )
      ) %>% 
      tidyr::unnest(.data$name_split) %>% 
      # extract evidence from comment_cofactor
      dplyr::mutate(evidence = stringr::str_extract(.data$name_split, pattern = '(?<=Evidence=).+(?=;|$)')) %>% 
      dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ';|\\{|\\}')) %>% 
      dplyr::mutate(note_evidence = str_extract(.data$note,
                                                pattern = "(?<=\\{).+(?=\\})")) %>% 
      dplyr::mutate(chebi_id = stringr::str_extract(
        .data$name_split,
        pattern = "(?<=CHEBI:)\\d+"
      )) %>% 
      dplyr::left_join(data_chebi_filtered, by = "chebi_id") %>% 
      # Filter for all ChEBI entries that contain a metal or that are in the reference data provided by protti
      # metal_list and metal_chebi_uniprot are provided by protti
      dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) | 
                      chebi_id %in% as.character(metal_chebi_uniprot$id))
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  
    # Extract comment_catalytic_activity information from UniProt
    if (show_progress == TRUE) {
      message("Extract comment_catalytic_activity information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    catalytic_activity_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$id, .data$comment_catalytic_activity) %>% 
      tidyr::drop_na(.data$comment_catalytic_activity) %>% 
      dplyr::mutate(catalytic_activity_split = stringr::str_extract_all(
        .data$comment_catalytic_activity,
        pattern = "(?<=CATALYTIC ACTIVITY:).+?(?=CATALYTIC ACTIVITY|$)")) %>% 
      tidyr::unnest(.data$catalytic_activity_split) %>% 
      dplyr::mutate(catalytic_activity_split = stringr::str_remove(
        stringr::str_trim(.data$catalytic_activity_split),
        pattern = "Reaction=")) %>% 
      dplyr::mutate(evidence = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=Evidence=)[^;]+(?=;)')) %>% 
      dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ';|\\{|\\}')) %>% 
      dplyr::mutate(physiological_direction = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=PhysiologicalDirection=).+(?=;$)')) %>%
      dplyr::mutate(physiological_direction = stringr::str_split(.data$physiological_direction, pattern = 'PhysiologicalDirection=')) %>% 
      tidyr::unnest(.data$physiological_direction) %>%
      dplyr::mutate(ec = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=EC=)[^;]+(?=;)')) %>% 
      dplyr::mutate(rhea = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=RHEA:)\\d+')) %>% 
      dplyr::mutate(chebi_id = stringr::str_extract_all(
        .data$catalytic_activity_split,
        pattern = "(?<=CHEBI:)\\d+"
      )) %>% 
      tidyr::unnest(.data$chebi_id) %>% 
      dplyr::left_join(data_chebi_filtered, by = "chebi_id") %>% 
      # Filter for all ChEBI entries that contain a metal or that are in the reference data provided by protti
      # metal_list and metal_chebi_uniprot are provided by protti
      dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) |
                      chebi_id %in% as.character(metal_chebi_uniprot$id))
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }

    # Check if there are any metal containing entries that are not yet part of the ChEBI dataset provided by protti
    # This could be an indirect indication that some of the manually added ChEBI entries (without formula but metal related)
    # are also missing.
    if(any(!(unique(c(cofactor_uniprot$chebi_id, catalytic_activity_uniprot$chebi_id)) %in% as.character(metal_chebi_uniprot$id)))){
      missing_chebi_ids <- cofactor_uniprot %>% 
        dplyr::distinct(.data$id, .data$chebi_id) %>% 
        dplyr::bind_rows(dplyr::distinct(catalytic_activity_uniprot, .data$id, .data$chebi_id)) %>% 
        dplyr::distinct() %>% 
        dplyr::filter(!(.data$chebi_id %in% as.character(metal_chebi_uniprot$id)))
      
      warning(strwrap("The following ChEBI IDs have been found in the comment_cofactor or comment_catalytic_activity 
                      column and have not yet been manually annotated in the reference data frame provided by protti. 
                      This could be an indicator that there are additional ChEBI IDs missing that do not contain a 
                      formula but are that are metal related IDs.
                      Please contact the package maintainer to let potentially missing ChEBI IDs be added.",
                      prefix = "\n", initial = ""),
              "\n",
              paste0(utils::capture.output(missing_chebi_ids), collapse = "\n"))
    }
    
    
    
    return(fmb_uniprot)
    
    
    
    
  }

# Fix example

# Extract feature metal binding information 
## Decide what to do with catalytic info, and ligand atom info
## Extract evidence correctly and annotate

# Extract comment cofactor and CA 
## if there is only a note make sure that a search for key terms is done to see if there is anything that is missing
## Decide what to do about physiological direction evidence and rhea

# Extract go information

# There might need to be some more manual annotation for go and uniprot chebi

# Extract ECO information
# filter go for "molecular_function" and "UniProtKB"
# fetch eco information and extract IDs based on manual assertion and IDs based on automatic assertion

# Combine information


