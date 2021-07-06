


map_peptides_on_structure <- function(peptide_data, uniprot_id, pdb_id, chain, start_in_pdb, end_in_pdb, map_value, file_format = ".cif", export_location, structure_file = NULL, show_progress = TRUE) {
  if (!file_format %in% c(".cif", ".pdb")){
    stop('Please provide either ".cif" or ".pdb" to the file_format argument.')
  }
  
  if (missing(structure_file)){
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }
    
    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }
    
    peptide_data_filter <- peptide_data %>% 
      tidyr::drop_na({{ pdb_id }}) %>% 
      dplyr::mutate({{map_value }} := round(scale_protti(c({{ map_value }}), method = "01")*50+50, digits = 2)) %>% # Scale values between 0 and 100 
      group_by({{ pdb_id }}, {{ chain }}) %>% 
      dplyr::mutate(residue = list(seq({{ start_in_pdb }}, {{ end_in_pdb }}))) %>% 
      dplyr::ungroup() %>% 
      tidyr::unnest(residue)
    
    # Extract the length of the B-factor replacement value. It is used later to ensure that the spacing stays constant in the structure file.
    length_replacement <- peptide_data_filter %>% 
      dplyr::distinct({{ map_value }}) %>% 
      dplyr::mutate(map_value = as.character({{ map_value }})) %>% 
      dplyr::mutate(length = nchar(.data$map_value)) %>% 
      dplyr::pull(.data$length) %>% 
      max()
    
    # Extract PDB IDs for fetching
    pdb_preparation <- peptide_data_filter %>% 
      dplyr::distinct({{ pdb_id }}, {{ uniprot_id }}) %>% 
      dplyr::group_by({{ pdb_id }}) %>% 
      dplyr::mutate(name = paste({{ uniprot_id }}, collapse = "_")) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(name = paste0({{ pdb_id }}, "_", .data$name)) %>% 
      dplyr::select(-{{ uniprot_id }}) %>% 
      dplyr::distinct()
    
    pdb_ids <- dplyr::pull(pdb_preparation, {{ pdb_id }}) 
    pdb_names <- dplyr::pull(pdb_preparation, .data$name) 
    
    if (file_format == ".cif"){
      batches <- purrr::map(
        .x = pdb_ids,
        .f = ~ paste0("https://files.rcsb.org/download/", .x, ".cif")
      )
      
      names(batches) <- pdb_names
      
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(batches), format = "  Fetching structures [:bar] :current/:total (:percent) :eta")
      }
      
      query_result <- purrr::map(
        .x = batches,
        .f = ~ {
          # query information from database
          if (!is.null(batches)) {
            query <- try_query(.x, as = "text", col_names = FALSE, quote = "")
          }
          if (show_progress == TRUE & "character" %in% class(query)) {
            pb$tick()
          }
          
          # if previous batch had a connection problem change batches to NULL, which breaks the mapping.
          if (!"character" %in% class(query)) {
            batches <<- NULL
          }
          
          # only proceed with data if it was correctly retrieved
          if ("character" %in% class(query)) {
            query %>%
              readr::read_tsv(col_names = FALSE, show_col_types = FALSE)
          }
        }
      )
      
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(query_result), format = "  Replacing B-factors [:bar] :current/:total (:percent) :eta")
      }
      
      purrr::map2(
        .x = query_result,
        .y = names(query_result),
        .f = ~ {
          result <- .x %>% 
            dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), X1, NA)) %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(b_factor = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+[:space:]+")[[1]][15], # extract b-factor values based on positions
                          chain = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][7],
                          residue = suppressWarnings(as.numeric(stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][9]))) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(pdb_id = stringr::str_extract(.y, pattern = "^\\w+(?=_)")) %>% 
            dplyr::mutate(b_factor = stringr::str_replace(.data$b_factor, pattern = "\\.", replacement = "\\\\\\.")) %>% 
            dplyr::left_join(peptide_data_filter, by = c("pdb_id" = rlang::as_name(rlang::enquo(pdb_id)),
                                                         "chain" = rlang::as_name(rlang::enquo(chain)),
                                                         "residue")) %>%
            dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$b_factor) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
            dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = length_replacement + 1)) %>%
            dplyr::mutate(atoms_mod = stringr::str_replace(.data$atoms, pattern = paste0("(?<= )", .data$b_factor), replacement = {{ map_value }})) %>%
            dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
            dplyr::select(.data$X1) %>%
            readr::write_tsv(file = paste0(export_location, .y, file_format), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)
            
          if (show_progress == TRUE) {
            pb$tick()
          }
          
          result
        }
      )    
      
      return(NULL)
      
    }
    if (file_format == ".pdb"){
      # same as above but for .pdb files. Indexing is different.
      batches <- purrr::map(
        .x = pdb_ids,
        .f = ~ paste0("https://files.rcsb.org/download/", .x, ".pdb")
      )
      
      names(batches) <- pdb_names
      
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(batches), format = "  Fetching structures [:bar] :current/:total (:percent) :eta")
      }
      
      query_result <- purrr::map(
        .x = batches,
        .f = ~ {
          # query information from database
          if (!is.null(batches)) {
            query <- suppressMessages(try_query(.x, as = "text", col_names = FALSE, quote = ""))
          }
          if (show_progress == TRUE & "character" %in% class(query)) {
            pb$tick()
          }
          
          # only proceed with data if it was correctly retrieved
          if ("character" %in% class(query)) {
            query %>%
              readr::read_tsv(col_names = FALSE, show_col_types = FALSE)
          }
        }
      )
      
      not_fetched <- query_result %>% 
        purrr::keep(is.null) %>% 
        names()
      
      if(length(not_fetched) != 0) {
        message('The following structures were not fetched, likely because no ".pdb" file is available. Try using the ".cif" format for these. ', paste(not_fetched, collapse = ", "))
      }
      
      query_result <- query_result %>% 
        purrr::discard(is.null)
      
      if (show_progress == TRUE) {
        pb <- progress::progress_bar$new(total = length(query_result), format = "  Replacing B-factors [:bar] :current/:total (:percent) :eta")
      }
      
      result <- purrr::map2(
        .x = query_result,
        .y = names(query_result),
        .f = ~ {
          result <- .x %>% 
            dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), X1, NA)) %>% 
            dplyr::mutate(chain = stringr::str_sub(.data$atoms, start = 22, end = 22),
                          residue = suppressWarnings(as.numeric(stringr::str_sub(.data$atoms, start = 23, end = 26)))) %>% 
            dplyr::mutate(pdb_id = stringr::str_extract(.y, pattern = "^\\w+(?=_)")) %>% 
            dplyr::left_join(peptide_data_filter, by = c("pdb_id" = rlang::as_name(rlang::enquo(pdb_id)),
                                                         "chain" = rlang::as_name(rlang::enquo(chain)),
                                                         "residue")) %>%
            dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$residue) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
            dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = 6, justify = "right")) %>%
            dplyr::mutate({{ map_value }} := ifelse(str_detect({{ map_value }}, pattern = "NA"), NA, {{ map_value }})) %>% 
            dplyr::mutate(atoms_mod = `str_sub<-`(.data$atoms, 61, 66, value = {{ map_value }})) %>% 
            dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
            dplyr::select(.data$X1) %>%
            readr::write_tsv(file = paste0(export_location, .y, file_format), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)
          
          if (show_progress == TRUE) {
            pb$tick()
          }
          
          result
        }
      ) 

      return(NULL)
      
    }
  }
  
  if (!missing(structure_file)){
    
    file_format <- stringr::str_sub(structure_file, start = -4, end = -1)
    
    if (!file_format %in% c(".cif", ".pdb")){
      stop('Please either provide a ".cif" or ".pdb" structure file.')
    }
    
    file_name <- stringr::str_extract(structure_file, pattern = "\\w+[:punct:]\\w+$")
    
    pdb_file <- readr::read_tsv(structure_file, col_names = FALSE, show_col_types = FALSE)
    
    peptide_data_filter <- peptide_data %>% 
      dplyr::mutate({{map_value }} := round(scale_protti(c({{ map_value }}), method = "01")*50+50, digits = 2)) %>% # Scale values between 0 and 100 
      group_by({{ chain }}) %>% 
      dplyr::mutate(residue = list(seq({{ start_in_pdb }}, {{ end_in_pdb }}))) %>% 
      dplyr::ungroup() %>% 
      tidyr::unnest(residue)
    
    length_replacement <- peptide_data_filter %>% 
      dplyr::distinct({{ map_value }}) %>% 
      dplyr::mutate(map_value = as.character({{ map_value }})) %>% 
      dplyr::mutate(length = nchar(.data$map_value)) %>% 
      dplyr::pull(.data$length) %>% 
      max()
    
    if(file_format == ".cif") {
      pdb_file %>% 
        dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), X1, NA)) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(b_factor = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+[:space:]+")[[1]][15], # extract b-factor values based on positions
                      chain = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][7],
                      residue = suppressWarnings(as.numeric(stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][9]))) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(b_factor = stringr::str_replace(.data$b_factor, pattern = "\\.", replacement = "\\\\\\.")) %>% 
        dplyr::left_join(peptide_data_filter, by = c("chain" = rlang::as_name(rlang::enquo(chain)),
                                                     "residue")) %>%
        dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$b_factor) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
        dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = length_replacement + 1)) %>%
        dplyr::mutate(atoms_mod = stringr::str_replace(.data$atoms, pattern = paste0("(?<= )", .data$b_factor), replacement = {{ map_value }})) %>%
        dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
        dplyr::select(.data$X1) %>%
        readr::write_tsv(file = paste0(export_location, "modified_", file_name), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)
      
      return(NULL)
    }
    
    if(file_format == ".pdb") {
      pdb_file %>% 
        dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), X1, NA)) %>% 
        dplyr::mutate(chain = stringr::str_sub(.data$atoms, start = 22, end = 22),
                      residue = suppressWarnings(as.numeric(stringr::str_sub(.data$atoms, start = 23, end = 26)))) %>% 
        dplyr::left_join(peptide_data_filter, by = c("chain" = rlang::as_name(rlang::enquo(chain)),
                                                     "residue")) %>%
        dplyr::mutate({{ map_value }} := ifelse(!is.na(.data$residue) & is.na({{ map_value }}), 0, {{ map_value }})) %>%
        dplyr::mutate({{ map_value }} := format(as.character({{ map_value }}), width = 6, justify = "right")) %>%
        dplyr::mutate({{ map_value }} := ifelse(str_detect({{ map_value }}, pattern = "NA"), NA, {{ map_value }})) %>% 
        dplyr::mutate(atoms_mod = `str_sub<-`(.data$atoms, 61, 66, value = {{ map_value }})) %>% 
        dplyr::mutate(X1 = ifelse(!is.na(.data$atoms_mod), .data$atoms_mod, .data$X1)) %>%
        dplyr::select(.data$X1) %>%
        readr::write_tsv(file = paste0(export_location, "modified_", file_name), quote = "none", escape = "none", col_names = FALSE, progress = FALSE)
        
      return(NULL)
    }

  }

}