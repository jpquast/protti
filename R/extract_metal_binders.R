#' Extract metal-bind protein information from UniProt
#'
#' Information of metal binding proteins is extracted from UniProt data retrieved with
#' \code{fetch_uniprot}. ChEBI IDs, potential sub-IDs for metal cations, binding site locations
#' in the protein and sub-ID evidence level (based on metal presence as cofactor) are extracted.
#'
#' @param data a data frame containing at least the input columns.
#' @param protein_id a character column in the \code{data} data frame that contains the protein
#' identifiers.
#' @param ft_metal a character column in the \code{data} data frame that contains the
#' feature metal binding information from UniProt.
#' @param chebi_cofactor a character column in the \code{data} data frame that contains the ChEBI
#' cofactor information from UniProt.
#' @param chebi_catalytic_activity a character column in the \code{data} data frame that contains
#' the ChEBI catalytic activity information from UniProt.
#' @param cc_cofactor a character cholumn in the \code{data} data frame that contains
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
#' \item{\code{metal_type}: }{Metal name extracted from \code{ft_metal} information.
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
#' \item{\code{source}: }{The source of the information, can be either \code{ft_metal}, \code{cc_cofactor} 
#' or \code{chebi_catalytic_activity}.}
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
#'  # protein_id = accession,
#'  # ft_metal = ft_metal,
#'  # chebi_cofactor = chebi_cofactor,
#'  # chebi_catalytic_activity = chebi_catalytic_activity,
#'  # cc_cofactor = cc_cofactor,
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
    if(!("ft_metal" %in% colnames(data_uniprot) & #TODO change name
         "cc_cofactor" %in% colnames(data_uniprot) & 
         "cc_catalytic_activity" %in% colnames(data_uniprot))){
      stop('Please include at least the columns "ft_metal", "cc_cofactor" and "cc_catalytic_activity" in "data_uniprot"!')
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
      
      data_chebi <- fetch_chebi(stars = c(2, 3))
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
    
    if (show_progress == TRUE) {
      message("Preparing annotation data frames ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    # Prepare ChEBI dataframe for annotation
    data_chebi_filtered <- data_chebi %>% 
      dplyr::filter(.data$type_name == "STANDARD") %>% 
      dplyr::mutate(chebi_id = as.character(.data$id)) %>% 
      dplyr::distinct(chebi_id, definition, name, formula, mass, charge) #TODO check if we should really provide all these columns
    
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
    
    # Create a data frame that maps which complex IDs belong to which metal IDs. This information 
    # can be created using the fmb_annotation_uniprot data that is provided with protti. It should
    # contain all relevant mappings between metal and complex IDs. 
    # TODO this should be overworked and replaced with the new version
    complex_id_metal_mapping <- fmb_annotation_uniprot %>% 
      # remove any character IDs since they are only manual annotations not used by any of the databases we match to.
      dplyr::filter(!stringr::str_detect(.data$complex_id, pattern = "[:alpha:]")) %>% 
      dplyr::distinct(.data$metal_id, .data$complex_id, .data$complex_sub_id) %>% 
      dplyr::rowwise() %>%
      # combine complex_id and complex_sub_ids into one column since they are all metal complexes.
      dplyr::mutate(complex_id = stringr::str_split(paste0(c(complex_id, complex_sub_id), collapse = ","), "\\||,")) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(complex_id) %>%
      dplyr::distinct(metal_id, complex_id) %>% 
      dplyr::filter(complex_id != "") %>% 
      # For complex IDs that have multiple metal IDs only keep the ones with the fewest sub IDs. 
      # This makes the assumption that the one with the fewest is always correct, which might not 
      # hold true anymore in the future if more IDs are added. 
      dplyr::mutate(subs_length = purrr::map(find_all_subs(data_chebi_relation, .data$metal_id),
                                             .f = ~ length(.x))
      ) %>% 
      dplyr::group_by(.data$complex_id) %>% 
      dplyr::filter(.data$subs_length == min(as.numeric(.data$subs_length))) %>% 
      # If there are still multiple metal IDs per complex ID left they are combined.
      # This should usually not be the case.
      dplyr::mutate(metal_id = paste0(.data$metal_id, collapse = "|")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct(.data$metal_id, .data$complex_id)
    
    # Create a data frame that contains all "pure" metal ChEBI IDs from UniProt. IDs that are sub IDs of 
    # other present IDs are filtered out. This data frame will be used in order to annotate complex_ids 
    # that are missing from the current fmb_annotation_uniprot data provided with protti. There should 
    # usually be no IDs that are not already part of fmb_annotation_uniprot so this is an attempt to 
    # rescue IDs in case of an outdated fmb_annotation_uniprot data frame.
    # TODO update the comment 
    chebi_metal_id <- metal_chebi_uniprot %>% 
      dplyr::distinct(.data$formula, .data$id) %>% 
      dplyr::mutate(id = as.character(.data$id)) %>% 
      dplyr::filter(.data$formula %in% unique(metal_list$symbol)) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(metal_sub_id = paste0(find_all_subs(data_chebi_relation, .data$id), collapse = ",")) %>% 
      dplyr::group_by(.data$formula) %>% 
      # We make the assumption that any IDs that do not have sub IDs but that are in a group that contains
      # sub IDs can be removed since they are covered by the sub IDs of that group. This holds true for now
      # but might change in the future.
      dplyr::filter(!(any(.data$metal_sub_id != "") & .data$metal_sub_id == "")) %>% 
      dplyr::mutate(metal_id = paste0(.data$id, collapse = "|")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct(.data$formula, .data$metal_id)
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
    
    # Extract ft_metal information from UniProt
    if (show_progress == TRUE) {
      message("Extract ft_metal information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    fmb_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$accession, .data$ft_metal) %>% 
      tidyr::drop_na(.data$ft_metal) %>% 
      # Extract names to be compatible with manual annotation
      dplyr::mutate(split_fmb = stringr::str_extract_all(
        .data$ft_metal,
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
      dplyr::mutate(eco = stringr::str_trim(.data$eco)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco_type = dplyr::case_when(.data$eco %in% manual_eco ~ "manual_assertion",
                                                .data$eco %in% automatic_eco ~ "automatic_assertion")) %>% 
      dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>% 
      # The data that will be joined is a dataset provided with protti
      dplyr::left_join(fmb_annotation_uniprot, by = "fmb_name") %>% 
      # Label cases in which there was no metal name
      dplyr::mutate(comment = dplyr::case_when(.data$no_metal_name & !is.na(.data$comment) ~ paste0(.data$comment, ' There was no metal name therefore "Divalent metal cation" was used.'),
                                               .data$no_metal_name & is.na(.data$comment) ~ 'There was no metal name therefore "Divalent metal cation" was used.',
                                               TRUE ~ .data$comment)) %>% 
      dplyr::select(-c(.data$ft_metal, .data$no_metal_name, .data$split_fmb, .data$evidence)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco = ifelse(is.na(.data$eco), " ", .data$eco),
                    evidence_source = ifelse(is.na(.data$evidence_source), " ", .data$evidence_source),
                    binding_mode = ifelse(is.na(.data$binding_mode), " ", .data$binding_mode)) %>% 
      # combine the binding_mode column beforehand to prevent duplicates
      # remove the ";" symbol and replace with ":" to not have conflicts with the separator
      dplyr::mutate(binding_mode = stringr::str_replace(.data$binding_mode, pattern = ";", replacement = ":")) %>% 
      dplyr::group_by(.data$accession, .data$fmb_name, .data$metal_identifier, .data$metal_position, .data$eco, .data$metal_id) %>%
      dplyr::mutate(binding_mode = paste0(.data$binding_mode, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # Combine all evidence source information for each ECO ID and metal ID
      dplyr::group_by(.data$accession, .data$fmb_name, .data$metal_identifier, .data$metal_position, .data$eco, .data$metal_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # combine the binding_mode further
      dplyr::group_by(.data$accession, .data$fmb_name, .data$metal_identifier, .data$metal_position, .data$metal_id) %>% 
      dplyr::mutate(binding_mode = paste0(.data$binding_mode, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$fmb_name, .data$metal_identifier, .data$metal_position, .data$metal_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";;"),
             eco = paste0(.data$eco, collapse = ";;"),
             eco_type = paste0(.data$eco_type, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(source = "fmb") %>% 
      dplyr::rename(note = .data$comment) %>% 
      dplyr::distinct() 
      

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
    
    # Check if there are any new terms not manually annotated in the annotation data frame
    if(any(!(fmb_uniprot$fmb_name %in% fmb_annotation_uniprot$fmb_name))){
      missing_names <- fmb_uniprot %>% 
        dplyr::distinct(.data$accession, .data$fmb_name) %>% 
        dplyr::filter(!(.data$fmb_name %in% fmb_annotation_uniprot$fmb_name))
        
      warning(strwrap("The following names have been found in the ft_metal column and have not yet been 
                      manually annotated in the reference data frame provided by protti. Please contact the package 
                      maintainer to let them be added.",
                      prefix = "\n", initial = ""),
              "\n",
              paste0(utils::capture.output(missing_names), collapse = "\n"))
    }
    
    # Extract cc_cofactor information from UniProt
    if (show_progress == TRUE) {
      message("Extract cc_cofactor information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    cofactor_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$accession, .data$cc_cofactor) %>% 
      tidyr::drop_na(.data$cc_cofactor) %>% 
      dplyr::mutate(cofactor_split = stringr::str_extract_all(
                                   .data$cc_cofactor,
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
      # extract evidence from cc_cofactor
      dplyr::mutate(evidence = stringr::str_extract(.data$name_split, pattern = '(?<=Evidence=).+(?=;|$)')) %>% 
      dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ';|\\{|\\}')) %>% 
      dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>% 
      tidyr::unnest(.data$evidence_split) %>% 
      tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>% 
      dplyr::mutate(eco = stringr::str_trim(.data$eco)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco_type = dplyr::case_when(.data$eco %in% manual_eco ~ "manual_assertion",
                                                .data$eco %in% automatic_eco ~ "automatic_assertion")) %>% 
      dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>% 
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
                      chebi_id %in% as.character(metal_chebi_uniprot$id)) %>% 
      # is_complex is TRUE if the corresponding ChEBI ID is a complex. We rescue an incomplete
      # dataset by checking for element symbols in addition to using the complex_id_metal_mapping 
      # data frame. Please note that even though we can rescue some cases entries without a formula 
      # cannot be rescued. 
      dplyr::mutate(is_complex = (!.data$formula %in% unique(metal_list$symbol) & !is.na(.data$formula))|
                      .data$chebi_id %in% as.character(complex_id_metal_mapping$complex_id)) %>% 
      dplyr::left_join(complex_id_metal_mapping, by = c("chebi_id" = "complex_id")) %>% 
      dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(metal_id = ifelse(.data$is_complex & is.na(.data$metal_id), 
                                      paste0(chebi_metal_id[chebi_metal_id$formula %in% extract_formula,]$metal_id, collapse = "|"), 
                                      .data$metal_id)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(complex_id = ifelse(!is.na(.data$metal_id), .data$chebi_id, NA)) %>% 
      dplyr::mutate(metal_id = ifelse(is.na(.data$metal_id), .data$chebi_id, .data$metal_id)) %>% 
      dplyr::mutate(evidence_source = stringr::str_replace(.data$evidence_source, pattern = "\\|", replacement = ",")) %>% 
      dplyr::select(-c(.data$cc_cofactor, .data$cofactor_split, .data$name_split, .data$evidence, .data$definition, .data$mass)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco = ifelse(is.na(.data$eco), " ", .data$eco),
                    evidence_source = ifelse(is.na(.data$evidence_source), " ", .data$evidence_source),
                    note = ifelse(is.na(.data$note), " ", .data$note),
                    note_evidence = ifelse(is.na(.data$note_evidence), " ", .data$note_evidence)) %>% 
      # Combine all evidence source information for each ECO ID and ChEBI ID
      # combine the note and note_evidence column beforehand to prevent duplicates 
      # This is a bit annoying due to the nature of the note column, which can contain multiple notes for the same evidence
      # or the same note for different evidences.
      # remove the ";" symbol and replace with ":" to not have conflicts with the separator
      dplyr::mutate(note = stringr::str_replace(.data$note, pattern = ";", replacement = ":")) %>%
      # First concatenate different notes for the same evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(note = paste0(.data$note, collapse = ","),
                    note_evidence = paste0(.data$note_evidence, collapse = ",")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # run condensation of notes again to concatenate different notes for the same condensed evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(note = paste0(.data$note, collapse = ";"),
                    note_evidence = paste0(.data$note_evidence, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";;"),
                    eco = paste0(.data$eco, collapse = ";;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # combine the note column further
      dplyr::group_by(.data$accession, .data$chebi_id) %>%
      dplyr::mutate(note = paste0(.data$note, collapse = ";;"),
                    note_evidence = paste0(.data$note_evidence, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(source = "cofactor")

    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }

    # Extract cc_catalytic_activity information from UniProt
    if (show_progress == TRUE) {
      message("Extract cc_catalytic_activity information from UniProt ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    catalytic_activity_uniprot <- data_uniprot %>% 
      dplyr::distinct(.data$accession, .data$cc_catalytic_activity) %>% 
      tidyr::drop_na(.data$cc_catalytic_activity) %>% 
      dplyr::mutate(catalytic_activity_split = stringr::str_extract_all(
        .data$cc_catalytic_activity,
        pattern = "(?<=CATALYTIC ACTIVITY:).+?(?=CATALYTIC ACTIVITY|$)")) %>% 
      tidyr::unnest(.data$catalytic_activity_split) %>% 
      dplyr::mutate(catalytic_activity_split = stringr::str_remove(
        stringr::str_trim(.data$catalytic_activity_split),
        pattern = "Reaction=")) %>% 
      dplyr::mutate(evidence = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=Evidence=)[^;]+(?=;)')) %>% 
      dplyr::mutate(evidence = stringr::str_remove_all(.data$evidence, pattern = ';|\\{|\\}')) %>% 
      dplyr::mutate(evidence_split = stringr::str_split(.data$evidence, pattern = ", ")) %>% 
      tidyr::unnest(.data$evidence_split) %>% 
      tidyr::separate(.data$evidence_split, into = c("eco", "evidence_source"), sep = "\\|", fill = "right") %>% 
      dplyr::mutate(eco = stringr::str_trim(.data$eco)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco_type = dplyr::case_when(.data$eco %in% manual_eco ~ "manual_assertion",
                                                .data$eco %in% automatic_eco ~ "automatic_assertion")) %>% 
      dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>% 
      dplyr::mutate(physiological_direction = stringr::str_extract(.data$catalytic_activity_split, pattern = '(?<=PhysiologicalDirection=).+(?=;$)')) %>%
      dplyr::mutate(physiological_direction = stringr::str_split(.data$physiological_direction, pattern = 'PhysiologicalDirection=')) %>% 
      tidyr::unnest(.data$physiological_direction) %>%
      dplyr::mutate(physiological_direction_evidence = stringr::str_extract(.data$physiological_direction, pattern = '(?<=Evidence=\\{)[^\\}]+(?=\\})')) %>%
      dplyr::mutate(physiological_direction = stringr::str_extract(.data$physiological_direction, pattern = '^[^;]+(?=;)')) %>%
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
                      chebi_id %in% as.character(metal_chebi_uniprot$id)) %>% 
      # is_complex is TRUE if the corresponding ChEBI ID is a complex. We rescue an incomplete
      # dataset by checking for element symbols in addition to using the complex_id_metal_mapping 
      # data frame. Please note that even though we can rescue some cases entries without a formula 
      # cannot be rescued. 
      dplyr::mutate(is_complex = (!.data$formula %in% unique(metal_list$symbol) & !is.na(.data$formula))|
                      .data$chebi_id %in% as.character(complex_id_metal_mapping$complex_id)) %>% 
      dplyr::left_join(complex_id_metal_mapping, by = c("chebi_id" = "complex_id")) %>% 
      dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(metal_id = ifelse(.data$is_complex & is.na(.data$metal_id), 
                                      paste0(chebi_metal_id[chebi_metal_id$formula %in% extract_formula,]$metal_id, collapse = "|"), 
                                      .data$metal_id)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(complex_id = ifelse(!is.na(.data$metal_id), .data$chebi_id, NA)) %>% 
      dplyr::mutate(metal_id = ifelse(is.na(.data$metal_id), .data$chebi_id, .data$metal_id)) %>% 
      dplyr::mutate(evidence_source = stringr::str_replace(.data$evidence_source, pattern = "\\|", replacement = ",")) %>% 
      dplyr::select(-c(.data$cc_catalytic_activity, .data$catalytic_activity_split, .data$evidence, .data$definition, .data$mass)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(eco = ifelse(is.na(.data$eco), " ", .data$eco),
                    evidence_source = ifelse(is.na(.data$evidence_source), " ", .data$evidence_source),
                    physiological_direction = ifelse(is.na(.data$physiological_direction), " ", .data$physiological_direction),
                    physiological_direction_evidence = ifelse(is.na(.data$physiological_direction_evidence), " ", .data$physiological_direction_evidence),
                    ec = ifelse(is.na(.data$ec), " ", .data$ec),
                    rhea = ifelse(is.na(.data$rhea), " ", .data$rhea)) %>% 
      # Combine all evidence source information for each ECO ID and Chebi ID
      # combine the physiological_direction, physiological_direction_evidence, rhea and ec column beforehand to prevent duplicates 
      # This is a bit annoying due to the nature of these column, which can contain multiple different entries for the same evidence
      # or the same entries for for different evidences.
      # First concatenate different entries for the same evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(physiological_direction = paste0(.data$physiological_direction, collapse = ","),
                    physiological_direction_evidence = paste0(.data$physiological_direction_evidence, collapse = ","),
                    rhea = paste0(.data$rhea, collapse = ","),
                    ec = paste0(.data$ec, collapse = ",")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # run condensation of entries again to concatenate different entries for the same condensed evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(physiological_direction = paste0(.data$physiological_direction, collapse = ";"),
                    physiological_direction_evidence = paste0(.data$physiological_direction_evidence, collapse = ";"),
                    rhea = paste0(.data$rhea, collapse = ";"),
                    ec = paste0(.data$ec, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";;"),
                    eco = paste0(.data$eco, collapse = ";;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # combine the additional columns further
      dplyr::group_by(.data$accession, .data$chebi_id) %>%
      dplyr::mutate(physiological_direction = paste0(.data$physiological_direction, collapse = ";;"),
                    physiological_direction_evidence = paste0(.data$physiological_direction_evidence, collapse = ";;"),
                    rhea = paste0(.data$rhea, collapse = ";;"),
                    ec = paste0(.data$ec, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(source = "catalytic_activity")
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }

    # Check if there are any metal containing entries that are not yet part of the ChEBI dataset provided by protti
    # This could be an indirect indication that some of the manually added ChEBI entries (without formula but metal related)
    # are also missing.
    if(any(!(unique(c(cofactor_uniprot$chebi_id, catalytic_activity_uniprot$chebi_id)) %in% as.character(metal_chebi_uniprot$id)))){
      missing_chebi_ids <- cofactor_uniprot %>% 
        dplyr::distinct(.data$accession, .data$chebi_id) %>% 
        dplyr::bind_rows(dplyr::distinct(catalytic_activity_uniprot, .data$accession, .data$chebi_id)) %>% 
        dplyr::distinct() %>% 
        dplyr::filter(!(.data$chebi_id %in% as.character(metal_chebi_uniprot$id)))
      
      warning(strwrap("The following ChEBI IDs have been found in the cc_cofactor or cc_catalytic_activity 
                      column and have not yet been manually annotated in the reference data frame provided by protti. 
                      This could be an indicator that there are additional ChEBI IDs missing that do not contain a 
                      formula but are that are metal related IDs.
                      Please contact the package maintainer to let potentially missing ChEBI IDs be added.",
                      prefix = "\n", initial = ""),
              "\n",
              paste0(utils::capture.output(missing_chebi_ids), collapse = "\n"))
    }
    
    # Extract information from QuickGO
    if (show_progress == TRUE) {
      message("Extract molecular_function information from QuickGO ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    mf_quickgo <- data_quickgo %>% 
      dplyr::filter(.data$gene_product_db == "UniProtKB" &
                      .data$go_aspect == "molecular_function") %>% 
      # The data that will be used for filtering is a dataset provided with protti
      # It contains all MF GO IDs that are associated with metals
      dplyr::filter(.data$go_term %in% metal_go_slim_subset$slims_from_id) %>% 
      # Filter GO data to only contain protein IDs also in the UniProt input
      dplyr::filter(.data$gene_product_id %in% data_uniprot$accession) %>% 
      dplyr::select(.data$gene_product_id,
                    .data$go_term,
                    .data$go_name,
                    .data$eco_id,
                    .data$reference,
                    .data$with_from,
                    .data$assigned_by) %>% 
      dplyr::distinct() %>% 
      dplyr::rename(eco = .data$eco_id,
                    accession = .data$gene_product_id) %>% 
      # join ChEBI annotations to data
      dplyr::left_join(dplyr::distinct(metal_go_slim_subset, .data$slims_from_id, .data$chebi_id, .data$database),
                       by = c("go_term" = "slims_from_id")) %>% 
      dplyr::mutate(eco_type = dplyr::case_when(.data$eco %in% manual_eco ~ "manual_assertion",
                                                .data$eco %in% automatic_eco ~ "automatic_assertion")) %>% 
      dplyr::mutate(eco_type = ifelse(is.na(.data$eco_type), "automatic_assertion", .data$eco_type)) %>% 
      # format ChEBI ID and annotate table properly
      dplyr::mutate(chebi_id = stringr::str_extract(.data$chebi_id, pattern = "\\d+")) %>% 
      dplyr::left_join(data_chebi_filtered, by = "chebi_id") %>% 
      # is_complex is TRUE if the corresponding ChEBI ID is a complex. We rescue an incomplete
      # dataset by checking for element symbols in addition to using the complex_id_metal_mapping 
      # data frame. Please note that even though we can rescue some cases entries without a formula 
      # cannot be rescued. 
      dplyr::mutate(is_complex = (!.data$formula %in% unique(metal_list$symbol) & !is.na(.data$formula))|
                      .data$chebi_id %in% as.character(complex_id_metal_mapping$complex_id)) %>% 
      dplyr::left_join(complex_id_metal_mapping, by = c("chebi_id" = "complex_id")) %>% 
      dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(metal_id = ifelse(.data$is_complex & is.na(.data$metal_id), 
                                      paste0(chebi_metal_id[chebi_metal_id$formula %in% extract_formula,]$metal_id, collapse = "|"), 
                                      .data$metal_id)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(complex_id = ifelse(!is.na(.data$metal_id), .data$chebi_id, NA)) %>% 
      dplyr::mutate(metal_id = ifelse(is.na(.data$metal_id), .data$chebi_id, .data$metal_id)) %>% 
      tidyr::unite(.data$reference, .data$with_from, col = "evidence_source", na.rm = TRUE) %>% 
      dplyr::mutate(evidence_source = stringr::str_replace_all(.data$evidence_source, pattern = "\\|", replacement = ",")) %>% 
      dplyr::select(-c(.data$definition, .data$mass)) %>% 
      dplyr::distinct() %>% 
      # Combine all evidence source information for each ECO ID and Chebi ID
      # combine the go_term, go_name, assigned_by and database beforehand to prevent duplicates 
      # This is a bit annoying due to the nature of these column, which can contain multiple different entries for the same evidence
      # or the same entries for for different evidences.
      # First concatenate different entries for the same evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(go_term = ifelse(paste0(.data$go_term, collapse = ",") == "NA",
                                     NA,
                                     paste0(.data$go_term, collapse = ",")),
                    go_name = ifelse(paste0(.data$go_name, collapse = ",") == "NA",
                                     NA,
                                     paste0(.data$go_name, collapse = ",")),
                    assigned_by = ifelse(paste0(.data$assigned_by, collapse = ",") == "NA",
                                         NA,
                                         paste0(.data$assigned_by, collapse = ",")),
                    database = ifelse(paste0(.data$database, collapse = ",") == "NA",
                                      NA,
                                      paste0(.data$database, collapse = ","))) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$eco, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # run condensation of entries again to concatenate different entries for the same condensed evidence_source
      dplyr::group_by(.data$accession, .data$eco, .data$evidence_source, .data$chebi_id) %>%
      dplyr::mutate(go_term = paste0(.data$go_term, collapse = ";"),
                    go_name = paste0(.data$go_name, collapse = ";"),
                    assigned_by = paste0(.data$assigned_by, collapse = ";"),
                    database = paste0(.data$database, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$chebi_id) %>% 
      dplyr::mutate(evidence_source = paste0(.data$evidence_source, collapse = ";;"),
                    eco = paste0(.data$eco, collapse = ";;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      # combine the additional columns further
      dplyr::group_by(.data$accession, .data$chebi_id) %>%
      dplyr::mutate(go_term = paste0(.data$go_term, collapse = ";;"),
                    go_name = paste0(.data$go_name, collapse = ";;"),
                    assigned_by = paste0(.data$assigned_by, collapse = ";;"),
                    database = paste0(.data$database, collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(source = "go_term")
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
    
    # Check if there are any complexed metals that are not yet part of the fmb_annotation_uniprot data frame 
    # provided by protti. This could be an indirect indication that some of the entries without formula but
    # that are metal related are also missing.
    cofactor_complex_ids <- cofactor_uniprot %>% 
      tidyr::drop_na(.data$complex_id)
    catalytic_activity_complex_ids <- catalytic_activity_uniprot %>% 
      tidyr::drop_na(.data$complex_id)
    go_complex_ids <- mf_quickgo %>% 
      tidyr::drop_na(.data$complex_id)

    if(any(!(unique(c(cofactor_complex_ids$complex_id, catalytic_activity_complex_ids$complex_id, go_complex_ids$complex_id)) %in% as.character(complex_id_metal_mapping$complex_id)))){
      missing_complex_ids <- cofactor_complex_ids %>% 
        dplyr::distinct(.data$accession, .data$complex_id) %>% 
        dplyr::bind_rows(dplyr::distinct(catalytic_activity_complex_ids, .data$accession, .data$complex_id)) %>% 
        dplyr::bind_rows(dplyr::distinct(go_complex_ids, .data$accession, .data$complex_id)) %>% 
        dplyr::distinct() %>% 
        dplyr::filter(!(.data$complex_id %in% as.character(complex_id_metal_mapping$complex_id)))
      
      warning(strwrap("The following ChEBI IDs contain a metal complex and have been found in the cc_cofactor, 
                  cc_catalytic_activity column or the GO molecular function data and have not yet been 
                  manually annotated as a complex in the reference data frame provided by protti. This 
                  could be an indicator that there are additional metal complex ChEBI IDs missing that 
                  do not contain a formula but are that are metal related complex IDs. Please contact 
                  the package maintainer to let potentially missing ChEBI IDs be added.",
                  prefix = "\n", initial = ""),
              "\n",
              paste0(utils::capture.output(missing_complex_ids), collapse = "\n"))
    }
    
    # Find ChEBI sub IDs for IDs from cofactor, catalytic activity and GO.
    if (show_progress == TRUE) {
      message("Find ChEBI sub IDs ... ", appendLF = FALSE)
      start_time <- Sys.time()
    }
    
    # Metal sub IDs
    metal_ids <- unique(c(cofactor_clean$metal_id, ca_clean$metal_id, go_clean$metal_id))
    metal_sub_id_mapping <- tibble::tibble(metal_id = metal_ids) %>% 
      dplyr::mutate(split_id = stringr::str_split(.data$metal_id, pattern = "\\|")) %>% 
      tidyr::unnest(.data$split_id) %>% 
      dplyr::mutate(metal_sub_id = purrr::map_chr(find_all_subs(data_chebi_relation, .data$split_id),
                                                  .f = ~ paste0(.x, collapse = ","))) %>% 
      dplyr::group_by(.data$metal_id) %>% 
      dplyr::mutate(metal_sub_id = paste0(metal_sub_id, collapse = "|")) %>% 
      dplyr::distinct(.data$metal_id, .data$metal_sub_id)
    
    # Complex sub IDs
    complex_ids <- unique(c(cofactor_clean$complex_id, ca_clean$complex_id, go_clean$complex_id))
    complex_sub_id_mapping <- tibble::tibble(complex_id = complex_ids) %>% 
      tidyr::drop_na(complex_id) %>% 
      dplyr::mutate(split_id = stringr::str_split(.data$complex_id, pattern = "\\|")) %>% 
      tidyr::unnest(.data$split_id) %>% 
      dplyr::mutate(complex_sub_id = purrr::map_chr(find_all_subs(data_chebi_relation, .data$split_id),
                                                    .f = ~ paste0(.x, collapse = ","))) %>% 
      dplyr::group_by(.data$complex_id) %>% 
      dplyr::mutate(complex_sub_id = paste0(complex_sub_id, collapse = "|")) %>% 
      dplyr::distinct(.data$complex_id, .data$complex_sub_id)
    
    if (show_progress == TRUE) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
    
    # Prepare data for combining
    
    ## Binding
    
    prep_fmb <- fmb_uniprot %>% 
      dplyr::mutate(complex_id = ifelse(.data$complex_id == "", " ", .data$complex_id),
                    complex_sub_id = ifelse(.data$complex_sub_id == "", " ", .data$complex_sub_id)) %>% 
      dplyr::group_by(.data$accession, .data$fmb_name, .data$metal_identifier) %>% 
      dplyr::mutate(metal_position =  paste0(.data$metal_position, collapse = ";"),
                    eco = paste0(.data$eco, collapse = ";"),
                    eco_type = paste0(.data$eco_type, collapse = ";"),
                    evidence_source = paste0(.data$evidence_source, collapse = ";"),
                    binding_mode = paste0(.data$binding_mode, collapse = ";")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(fmb_name = paste0(.data$fmb_name, collapse = ";;"),
                    metal_identifier =  paste0(.data$metal_identifier, collapse = ";;"),
                    metal_position =  paste0(.data$metal_position, collapse = ";;"),
                    eco = paste0(.data$eco, collapse = ";;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;"),
                    evidence_source = paste0(.data$evidence_source, collapse = ";;"),
                    binding_mode = paste0(.data$binding_mode, collapse = ";;"),
                    complex_id = paste0(.data$complex_id, collapse = ";;"),
                    complex_sub_id = paste0(.data$complex_sub_id, collapse = ";;"),
                    note = paste0(.data$note[.data$note != ""], collapse = ";;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(.data$accession,
                    .data$fmb_name, 
                    .data$metal_identifier,
                    .data$metal_position,
                    .data$eco, 
                    .data$eco_type, 
                    .data$evidence_source,
                    .data$metal_id,
                    .data$metal_sub_id,
                    .data$complex_id,
                    .data$complex_sub_id,
                    .data$note,
                    .data$source) %>% 
      dplyr::distinct()
    
    ## Cofactor
    
    prep_cofactor <- cofactor_uniprot %>% 
      dplyr::left_join(metal_sub_id_mapping, by = "metal_id") %>% 
      dplyr::left_join(complex_sub_id_mapping, by = "complex_id") %>% 
      dplyr::mutate(complex_id = ifelse(is.na(.data$complex_id), " ", .data$complex_id),
                    complex_sub_id = ifelse(is.na(.data$complex_sub_id), " ", .data$complex_sub_id)) %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(complex_id = paste0(.data$complex_id, collapse = ";;;"),
                    complex_sub_id = paste0(.data$complex_sub_id, collapse = ";;;"),
                    charge = paste0(.data$charge, collapse = ";;;"),
                    formula =  paste0(.data$formula, collapse = ";;;"),
                    extract_formula =  paste0(.data$extract_formula, collapse = ";;;"),
                    is_complex = paste0(.data$is_complex, collapse = ";;;"),
                    name = paste0(.data$name, collapse = ";;;"),
                    chebi_id = paste0(.data$chebi_id, collapse = ";;;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(eco = paste0(.data$eco, collapse = ";;;"),
                    evidence_source = paste0(.data$evidence_source, collapse = ";;;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;;"),
                    note = paste0(.data$note, collapse = ";;;"),
                    note_evidence = paste0(.data$note_evidence, collapse = ";;;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::select(.data$accession,
                    .data$eco,
                    .data$eco_type,
                    .data$evidence_source,
                    .data$metal_id,
                    .data$metal_sub_id,
                    .data$complex_id,
                    .data$complex_sub_id,
                    .data$note,
                    .data$source) %>% 
      dplyr::distinct()
    
    ## Catalytic activity
    
    prep_ca <- catalytic_activity_uniprot %>% 
      dplyr::left_join(metal_sub_id_mapping, by = "metal_id") %>% 
      dplyr::left_join(complex_sub_id_mapping, by = "complex_id") %>% 
      dplyr::mutate(complex_id = ifelse(is.na(.data$complex_id), " ", .data$complex_id),
                    complex_sub_id = ifelse(is.na(.data$complex_sub_id), " ", .data$complex_sub_id)) %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(complex_id = paste0(.data$complex_id, collapse = ";;;"),
                    complex_sub_id = paste0(.data$complex_sub_id, collapse = ";;;"),
                    charge = paste0(.data$charge, collapse = ";;;"),
                    formula =  paste0(.data$formula, collapse = ";;;"),
                    is_complex = paste0(.data$is_complex, collapse = ";;;"),
                    name = paste0(.data$name, collapse = ";;;"),
                    chebi_id = paste0(.data$chebi_id, collapse = ";;;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(eco = paste0(.data$eco, collapse = ";;;"),
                    evidence_source = paste0(.data$evidence_source, collapse = ";;;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;;"),
                    physiological_direction = paste0(.data$physiological_direction, collapse = ";;;"),
                    physiological_direction_evidence = paste0(.data$physiological_direction_evidence, collapse = ";;;"),
                    ec = paste0(.data$ec, collapse = ";;;"),
                    rhea = paste0(.data$rhea, collapse = ";;;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct() %>% 
      dplyr::select(.data$accession,
                    .data$eco,
                    .data$eco_type,
                    .data$evidence_source,
                    .data$metal_id,
                    .data$metal_sub_id,
                    .data$complex_id,
                    .data$complex_sub_id,
                    .data$source) %>% 
      dplyr::distinct()
    
    ## Molecular function gene ontology
    
    prep_go <- mf_quickgo %>% 
      dplyr::left_join(metal_sub_id_mapping, by = "metal_id") %>% 
      dplyr::left_join(complex_sub_id_mapping, by = "complex_id") %>% 
      dplyr::mutate(complex_id = ifelse(is.na(.data$complex_id), " ", .data$complex_id),
                    complex_sub_id = ifelse(is.na(.data$complex_sub_id), " ", .data$complex_sub_id)) %>% 
      dplyr::select(.data$accession,
                    .data$go_term,
                    .data$go_name,
                    .data$eco,
                    .data$eco_type,
                    .data$evidence_source,
                    .data$metal_id,
                    .data$metal_sub_id,
                    .data$complex_id,
                    .data$complex_sub_id,
                    .data$source) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(.data$accession, .data$metal_id) %>% 
      dplyr::mutate(go_term = paste0(.data$go_term, collapse = ";;;"),
                    go_name = paste0(.data$go_name, collapse = ";;;"),
                    eco = paste0(.data$eco, collapse = ";;;"),
                    evidence_source = paste0(.data$evidence_source, collapse = ";;;"),
                    eco_type = paste0(.data$eco_type, collapse = ";;;"),
                    complex_id = paste0(.data$complex_id, collapse = ";;;"),
                    complex_sub_id = paste0(.data$complex_sub_id, collapse = ";;;")) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct()
    
    
    return(prep_go)
    
    
  }

# Fix example
# Update required r packages in function if not required anymore

# Extract feature metal binding information 
## Decide what to do with catalytic info, and ligand atom info (keep for now)
## Extract evidence correctly and annotate

# Extract comment cofactor and CA 
## if there is only a note make sure that a search for key terms is done to see if there is anything that is missing
## Decide what to do about physiological direction evidence and rhea (keep for now)

# Part 1
# Retrieve and check information

# Part 2
# Extract information

# Part 3 
# Combine information



# Combine information
