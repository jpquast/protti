#' Creates a contact map of all atoms from a structure file (using parallel processing)
#'
#' This function is a wrapper around `create_structure_contact_map()` that allows the use of all
#' system cores for the creation of contact maps. Alternatively, it can be used for sequential
#' processing of large datasets. The benefit of this function over `create_structure_contact_map()`
#' is that it processes contact maps in batches, which is recommended for large datasets. If used
#' for parallel processing it should only be used on systems that have enough memory available.
#' Workers can either be set up manually before running the function with
#' `future::plan(multisession)` or automatically by the function (maximum number of workers
#' is 12 in this case). If workers are set up manually the `processing_type` argument should
#' be set to "parallel manual". In this case workers can be terminated after completion with 
#' `future::plan(sequential)`.
#'
#' @param data a data frame containing at least a column with PDB ID information of which the name
#' can be provided to the \code{id} argument. If only this column is provided, all atom or residue
#' distances are calculated. Additionally, a chain column can be present in the data frame of which
#' the name can be provided to the \code{chain} argument. If chains are provided, only distances
#' of this chain relative to the rest of the structure are calculated. Multiple chains can be
#' provided in multiple rows. If chains are provided for one structure but not for another, the
#' rows should contain NAs. Furthermore, specific residue positions can be provided in start and
#' end columns if the selection should be further reduced. It is not recommended to create full
#' contact maps for more than a few structures due to time and memory limitations. If contact maps are
#' created only for small regions it is possible to create multiple maps at once. By default distances
#' of regions provided in this data frame to the complete structure are computed. If distances of regions
#' from this data frame to another specific subset of regions should be computed, the second subset
#' of regions can be provided through the optional \code{data2} argument.
#' @param data2 optional, a data frame that contains a subset of regions for which distances to regions
#' provided in the \code{data} data frame should be computed. If regions from the \code{data} data
#' frame should be compared to the whole structure, data2 does not need to be provided.
#' This data frame should have the same structure and column names as the \code{data} data frame.
#' @param id a character column in the \code{data} data frame that contains PDB or UniProt IDs for
#' structures or AlphaFold predictions of which contact maps should be created. If a structure not
#' downloaded directly from PDB is provided (i.e. a locally stored structure file) to the
#' \code{structure_file} argument, this column should contain "my_structure" as content.
#' @param chain optional, a character column in the \code{data} data frame that contains chain
#' identifiers for the structure file. Identifiers defined by the structure author should be used.
#' Distances will be only calculated between the provided chains and the rest of the structure.
#' @param start_in_pdb optional, a numeric column in the \code{data} data frame that contains
#' start positions of regions which for distances should be calculated. This needs to be always
#' provided in combination with a corresponding end position in \code{end_in_pdb} and chain in
#' \code{chain}. The position should match the positioning defined by the structure author. For
#' PDB structures this information can be obtained from the \code{find_peptide_in_structure}
#' function. The corresponding column in the output is called \code{auth_seq_id_start}. If an
#' AlphaFold prediction is provided, UniProt positions should be used.
#' @param end_in_pdb optional, a numeric column in the \code{data} data frame that contains end
#' positions of regions which for distances should be calculated. This needs to be always provided
#' in combination with a corresponding start position in \code{start_in_pdb} and chain in
#' \code{chain}. The position should match the positioning defined by the structure author. For
#' PDB structures this information can be obtained from the \code{find_peptide_in_structure}
#' function. The corresponding column in the output is called \code{auth_seq_id_end}. If an
#' AlphaFold prediction is provided, UniProt positions should be used.
#' @param distance_cutoff a numeric value specifying the distance cutoff in Angstrom. All values
#' for pairwise comparisons are calculated but only values smaller than this cutoff will be
#' returned in the output. If a cutoff of e.g. 5 is selected then only residues with a distance of
#' 5 Angstrom and less are returned. Using a small value can reduce the size of the contact map
#' drastically and is therefore recommended. The default value is 10.
#' @param pdb_model_number_selection a numeric vector specifying which models from the structure
#' files should be considered for contact maps. E.g. NMR models often have many models in one file.
#' The default for this argument is c(0, 1). This means the first model of each structure file is
#' selected for contact map calculations. For AlphaFold predictions the model number is 0
#' (only .pdb files), therefore this case is also included here.
#' @param return_min_residue_distance a logical value that specifies if the contact map should be
#' returned for all atom distances or the minimum residue distances. Minimum residue distances are
#' smaller in size. If atom distances are not strictly needed it is recommended to set this
#' argument to TRUE. The default is TRUE.
#' @param export a logical value that indicates if contact maps should be exported as ".csv". The
#' name of the file will be the structure ID. Default is \code{export = FALSE}.
#' @param export_location optional, a character value that specifies the path to the location in
#' which the contact map should be saved if \code{export = TRUE}. If left empty, they will be
#' saved in the current working directory. The location should be provided in the following format
#' "folderA/folderB".
#' @param split_n a numeric value that specifies the number of structures that should be included
#' in each batch. Default is 40.
#' @param processing_type a character value that is either "parallel" for parallel processing or
#' "sequential" for sequential processing. Alternatively it can also be "parallel manual" in this
#' case you have to set up the number of cores on your own using the `future::plan(multisession)`
#' function. The default is "parallel".
#'
#' @return A list of contact maps for each PDB or UniProt ID provided in the input is returned.
#' If the \code{export} argument is TRUE, each contact map will be saved as a ".csv" file in the
#' current working directory or the location provided to the \code{export_location} argument.
#' @import dplyr
#' @importFrom purrr map2 flatten
#' @importFrom readr write_csv
#' @importFrom magrittr %>%
#' @importFrom rlang .data as_name enquo
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data <- data.frame(
#'   pdb_id = c("6NPF", "1C14", "3NIR"),
#'   chain = c("A", "A", NA),
#'   start = c(1, NA, NA),
#'   end = c(10, NA, NA)
#' )
#'
#' # Create contact map
#' contact_maps <- parallel_create_structure_contact_map(
#'   data = data,
#'   id = pdb_id,
#'   chain = chain,
#'   start_in_pdb = start,
#'   end_in_pdb = end,
#'   split_n = 1,
#' )
#'
#' str(contact_maps[["3NIR"]])
#'
#' contact_maps
#' }
parallel_create_structure_contact_map <- function(data,
                                                  data2 = NULL,
                                                  id,
                                                  chain = NULL,
                                                  start_in_pdb = NULL,
                                                  end_in_pdb = NULL,
                                                  distance_cutoff = 10,
                                                  pdb_model_number_selection = c(0, 1),
                                                  return_min_residue_distance = TRUE,
                                                  export = FALSE,
                                                  export_location = NULL,
                                                  split_n = 40,
                                                  processing_type = "parallel") {
  . <- NULL
  dependency_test <- c(
    furrr = !requireNamespace("furrr", quietly = TRUE),
    future = !requireNamespace("future", quietly = TRUE),
    parallel = !requireNamespace("parallel", quietly = TRUE),
    drc = !requireNamespace("drc", quietly = TRUE)
  )
  if (any(dependency_test)) {
    dependency_name <- names(dependency_test[dependency_test == TRUE])
    if (length(dependency_name) == 1) {
      stop("Package \"",
        paste(dependency_name),
        "\" is needed for this function to work. Please install it.",
        call. = FALSE
      )
    } else {
      stop("Packages \"",
        paste(dependency_name, collapse = "\" and \""),
        "\" are needed for this function to work. Please install them.",
        call. = FALSE
      )
    }
  }

  terminate <- FALSE

  if (processing_type == "parallel") {
    message("Setting up workers ... ", appendLF = FALSE)
    terminate <- TRUE
    n_cores <- parallel::detectCores()
    n_cores <- ifelse(n_cores > 12, 12, n_cores)
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
    message("DONE", appendLF = TRUE)
  }

  # if data2 was provided make sure that only IDs that overlap are retained
  if (!missing(data2)) {
    data <- data %>%
      dplyr::filter({{ id }} %in% unique(dplyr::pull(data2, {{ id }})))

    data2 <- data2 %>%
      dplyr::filter({{ id }} %in% unique(dplyr::pull(data, {{ id }})))
  }

  unique_ids <- unique(dplyr::pull(data, {{ id }}))

  pieces <- rep(1:ceiling(length(unique_ids) / split_n), split_n)[1:length(unique_ids)]

  pieces_mapping <- tibble::tibble({{ id }} := unique(dplyr::pull(data, {{ id }})), piece = pieces)

  input <- data %>%
    dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(id))) %>%
    split(dplyr::pull(., .data$piece))

  if (!missing(data2)) {
    input2 <- data2 %>%
      dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(id))) %>%
      split(dplyr::pull(., .data$piece))
  }

  message("Performing contact map calculations ... ", appendLF = FALSE)

  if (missing(data2)) {
    result <- furrr::future_map(
      .x = input,
      .f = ~ protti::create_structure_contact_map(
        data = .x,
        id = {{ id }},
        chain = {{ chain }},
        start_in_pdb = {{ start_in_pdb }},
        end_in_pdb = {{ end_in_pdb }},
        distance_cutoff = distance_cutoff,
        pdb_model_number_selection = pdb_model_number_selection,
        return_min_residue_distance = return_min_residue_distance,
        show_progress = FALSE
      ),
      .options = furrr::furrr_options(globals = FALSE)
    )

    message("DONE", appendLF = TRUE)
  }

  if (!missing(data2)) {
    result <- furrr::future_map2(
      .x = input,
      .y = input2,
      .f = ~ protti::create_structure_contact_map(
        data = .x,
        data2 = .y,
        id = {{ id }},
        chain = {{ chain }},
        start_in_pdb = {{ start_in_pdb }},
        end_in_pdb = {{ end_in_pdb }},
        distance_cutoff = distance_cutoff,
        pdb_model_number_selection = pdb_model_number_selection,
        return_min_residue_distance = return_min_residue_distance,
        show_progress = FALSE
      ),
      .options = furrr::furrr_options(globals = FALSE)
    )

    message("DONE", appendLF = TRUE)
  }

  if (terminate == TRUE) {
    future::plan(future::sequential)
  }

  result <- result %>%
    purrr::flatten()

  if (export == FALSE) {
    return(result)
  } else {
    message("Exporting contact maps ... ", appendLF = FALSE)
    # make sure export location is correct if or if not provided.
    if (missing(export_location)) {
      export_location <- ""
    } else {
      export_location <- paste0(export_location, "/")
    }
    purrr::map2(
      .x = result,
      .y = names(result),
      .f = ~ {
        readr::write_csv(x = .x, file = paste0(export_location, .y, "_contact_map.csv"), progress = FALSE)
      }
    )
    message("DONE", appendLF = TRUE)
    return(invisible(NULL))
  }
}