#' Fetch structure information from RCSB
#'
#' Fetches structure metadata from RCSB. If you want to retrieve atom data such as positions, use
#' the function \code{fetch_pdb_structure()}.
#'
#' @param pdb_ids a character vector of PDB identifiers.
#' @param batchsize a numeric value that specifies the number of structures to be processed in a
#' single query. Default is 100.
#' @param show_progress a logical value that indicates if a progress bar will be shown. Default is
#' TRUE.
#'
#' @return A data frame that contains structure metadata for the PDB IDs provided. The data frame
#' contains some columns that might not be self explanatory.
#'
#' * auth_asym_id: Chain identifier provided by the author of the structure in order to
#' match the identification used in the publication that describes the structure.
#' * label_asym_id: Chain identifier following the standardised convention for mmCIF files.
#' * entity_beg_seq_id, ref_beg_seq_id, length, pdb_sequence: \code{entity_beg_seq_id} is a
#' position in the structure sequence (\code{pdb_sequence}) that matches the position given in
#' \code{ref_beg_seq_id}, which is a position within the protein sequence (not included in the
#' data frame). \code{length} identifies the stretch of sequence for which positions match
#' accordingly between structure and protein sequence. \code{entity_beg_seq_id} is a residue ID
#' based on the standardised convention for mmCIF files.
#' * auth_seq_id: Residue identifier provided by the author of the structure in order to
#' match the identification used in the publication that describes the structure. This character
#' vector has the same length as the \code{pdb_sequence} and each position is the identifier for
#' the matching amino acid position in \code{pdb_sequence}. The contained values are not
#' necessarily numbers and the values do not have to be positive.
#' * modified_monomer: Is composed of first the composition ID of the modification, followed
#' by the \code{label_seq_id} position. In parenthesis are the parent monomer identifiers as
#' they appear in the sequence.
#' * ligand_*: Any column starting with the \code{ligand_*} prefix contains information about
#' the position, identity and donors for ligand binding sites. If there are multiple entities of
#' ligands they are separated by "|". Specific donor level information is separated by ";".
#' * secondar_structure: Contains information about helix and sheet secondary structure elements.
#' Individual regions are separated by ";".
#' * unmodeled_structure: Contains information about unmodeled or partially modeled regions in
#' the model. Individual regions are separated by ";".
#' * auth_seq_id_original: In some cases the sequence positions do not match the number of residues
#' in the sequence either because positions are missing or duplicated. This always coincides with modified
#' residues, however does not always occur when there is a modified residue in the sequence. This column
#' contains the original \code{auth_seq_id} information that does not have these positions corrected.
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom R.utils insert
#' @importFrom httr modify_url
#' @importFrom stringr str_replace_all
#' @importFrom curl has_internet
#' @importFrom utils URLencode
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \donttest{
#' pdb <- fetch_pdb(pdb_ids = c("6HG1", "1E9I", "6D3Q", "4JHW"))
#'
#' head(pdb)
#' }
fetch_pdb <- function(pdb_ids, batchsize = 100, show_progress = TRUE) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  . <- NULL
  # query that is used for fetching information
  query <- 'query={
  entries(entry_ids: ["pdb_ids"]) {
    rcsb_id
    struct_keywords {
      pdbx_keywords
    }
    exptl {
      method
    }
    exptl_crystal_grow {
      pH
      temp
      method
    }
    rcsb_binding_affinity {
      comp_id
      value
    }
    rcsb_entry_info {
      experimental_method
      assembly_count
      resolution_combined
      inter_mol_metalic_bond_count
    }
    pdbx_nmr_exptl_sample_conditions {
      ionic_strength
      pH
      temperature
    }
    pdbx_nmr_refine {
      method
    }
    pdbx_nmr_exptl {
      type
    }
    polymer_entities {
      polymer_entity_instances{
        rcsb_ligand_neighbors{
            atom_id
            auth_seq_id
            comp_id
            ligand_asym_id
            ligand_atom_id
            ligand_comp_id
            ligand_entity_id
            ligand_is_bound
            seq_id
        }
        rcsb_polymer_instance_feature{
          name
          type
          feature_positions{
            beg_comp_id
            beg_seq_id
            end_seq_id
          }
        }
      	rcsb_polymer_entity_instance_container_identifiers {
          asym_id
          auth_asym_id
          entry_id
      		auth_to_entity_poly_seq_mapping
    		}
      }
      rcsb_polymer_entity_feature {
      type
      name
        feature_positions {
          beg_comp_id
          beg_seq_id
          end_seq_id
        }
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        rcsb_artifact_monomer_count
        rcsb_conflict_count
        rcsb_deletion_count
        rcsb_insertion_count
        rcsb_mutation_count
        rcsb_non_std_monomer_count
        rcsb_non_std_monomers
      }
      rcsb_polymer_entity{
        pdbx_mutation
      }
      rcsb_entity_source_organism{
        ncbi_scientific_name
        ncbi_taxonomy_id
      }
      rcsb_polymer_entity_container_identifiers {
        entry_id
        auth_asym_ids
      }
      rcsb_polymer_entity_align {
      aligned_regions {
        entity_beg_seq_id
        ref_beg_seq_id
        length
      }
      reference_database_accession
      reference_database_isoform
      reference_database_name
    }
      uniprots {
      rcsb_uniprot_container_identifiers {
        uniprot_id
      }
      rcsb_uniprot_protein {
        name {
          value
        }
      }
    }
    }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers{
        auth_asym_ids
        entry_id
      }
      nonpolymer_comp {
        chem_comp {
          id
          type
          formula_weight
          name
          formula
        }
      }
    }
  }
}'

  # remove NA values
  pdb_ids <- pdb_ids[!is.na(pdb_ids)]
  if (length(pdb_ids) == 0) {
    stop("No PDB IDs were provided.")
  }
  # split pdb_ids into batches
  batches <- split(pdb_ids, ceiling(seq_along(pdb_ids) / batchsize))
  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(batches),
      format = "[1/6] Fetching PDB entries [:bar] (:percent) :eta"
    )
  }

  # query information from database

  query_result <- purrr::map(batches, function(x) {
    pdb_ids <- paste0(x, collapse = '", "')
    full_query <- stringr::str_replace_all(query, pattern = "pdb_ids", replacement = pdb_ids)
    url_encode_query <- utils::URLencode(full_query) %>%
      stringr::str_replace_all(pattern = "\\[", replacement = "%5B") %>%
      stringr::str_replace_all(pattern = "\\]", replacement = "%5D")

    query <- try_query(
      httr::modify_url("https://data.rcsb.org/graphql",
        query = url_encode_query
      ),
      type = "application/json",
      simplifyDataFrame = TRUE
    )

    if (show_progress == TRUE) {
      pb$tick()
    }
    # only proceed with data if it was correctly retrieved
    if ("list" %in% class(query)) {
      query <- query %>%
        purrr::flatten() %>%
        as.data.frame(stringsAsFactors = FALSE)
    }
    query
  })

  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>%
    purrr::keep(.p = ~ is.character(.x))

  if (length(error_list) != 0) {
    error_table <- tibble::tibble(
      id = paste0("IDs: ", as.numeric(names(error_list)) * batchsize - batchsize + 1, " to ", as.numeric(names(error_list)) * batchsize),
      error = unlist(error_list)
    ) %>%
      dplyr::distinct()

    message("The following IDs have not been retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output

  query_result <- query_result %>%
    purrr::keep(.p = ~ !is.character(.x))

  if (length(query_result) == 0) {
    message("No valid information could be retrieved!")
    return(invisible(NULL))
  }

  query_result <- purrr::map_dfr(
    .x = query_result,
    .f = ~.x
  ) %>%
    distinct() # make sure entries are unique otherwise there will be problems with unnesting

  if (nrow(query_result) == 0) {
    stop("None of the provided IDs could be retrieved!")
  }

  queried_ids <- unique(query_result$entries.rcsb_id)
  not_retrieved <- setdiff(stringr::str_to_lower(pdb_ids), stringr::str_to_lower(queried_ids))
  if (length(not_retrieved) > 0) {
    message("The following IDs have not been retrieved:")
    message(paste0(utils::capture.output(not_retrieved), collapse = "\n"))
  }

  # process information from database
  if (show_progress == TRUE) {
    message("[2/6] Extract experimental conditions ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  query_result_clean <- query_result %>%
    dplyr::bind_cols(
      pdb_ids = .$entries.rcsb_id,
      structure_keywords = .$entries.struct_keywords,
      .$entries.rcsb_entry_info
    ) %>%
    dplyr::select(-c(
      "entries.rcsb_id",
      "entries.struct_keywords",
      "entries.rcsb_entry_info"
    )) %>%
    tidyr::unnest("entries.exptl") %>%
    dplyr::rename(structure_method = "method")

  crystal_growth_info <- query_result_clean %>%
    dplyr::select("pdb_ids", "entries.exptl_crystal_grow") %>%
    tidyr::unnest("entries.exptl_crystal_grow")

  # make sure that the data is complete even if there is no crystal structure
  should_not_be_here <- colnames(crystal_growth_info)[!colnames(crystal_growth_info) %in% c(
    "pdb_ids",
    "method",
    "pH",
    "temp"
  )]
  should_be_here <- c(
    "pdb_ids",
    "method",
    "pH",
    "temp"
  )[!c(
    "pdb_ids",
    "method",
    "pH",
    "temp"
  ) %in% colnames(crystal_growth_info)]

  crystal_growth_info <- crystal_growth_info %>%
    dplyr::select(-all_of(should_not_be_here)) %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here),
        nrow = nrow(crystal_growth_info)
      )),
      should_be_here
    )) %>%
    dplyr::rename(
      pH_crystallisation = "pH",
      method_crystallisation = "method",
      temp_crystallisation = "temp"
    )

  resolution_info <- query_result_clean %>%
    dplyr::select("pdb_ids", "resolution_combined") %>%
    tidyr::unnest("resolution_combined")

  nmr_info <- query_result_clean %>%
    dplyr::select(
      "pdb_ids",
      "entries.pdbx_nmr_exptl",
      "entries.pdbx_nmr_exptl_sample_conditions",
      "entries.pdbx_nmr_refine"
    ) %>%
    tidyr::unnest("entries.pdbx_nmr_exptl") %>%
    tidyr::unnest("entries.pdbx_nmr_exptl_sample_conditions") %>%
    tidyr::unnest("entries.pdbx_nmr_refine")

  # make sure that the data is complete even if there is no NMR structure
  should_not_be_here <- colnames(nmr_info)[!colnames(nmr_info) %in% c(
    "pdb_ids",
    "type",
    "pH",
    "temperature",
    "method",
    "ionic_strength"
  )]
  should_be_here <- c(
    "pdb_ids",
    "type",
    "pH",
    "temperature",
    "method",
    "ionic_strength"
  )[!c(
    "pdb_ids",
    "type",
    "pH",
    "temperature",
    "method",
    "ionic_strength"
  ) %in% colnames(nmr_info)]

  nmr_info <- nmr_info %>%
    dplyr::select(-all_of(should_not_be_here)) %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here),
        nrow = nrow(nmr_info)
      )),
      should_be_here
    )) %>%
    dplyr::rename(
      type_nmr = "type",
      pH_nmr = "pH",
      temp_nmr = "temperature",
      method_nmr = "method",
      ionic_strength_nmr = "ionic_strength"
    )

  rcsb_binding_affinity <- query_result_clean %>%
    dplyr::select("pdb_ids", "entries.rcsb_binding_affinity") %>%
    tidyr::unnest("entries.rcsb_binding_affinity")

  # make sure that the data is complete even if there is no affinity information
  should_not_be_here <- colnames(rcsb_binding_affinity)[!colnames(rcsb_binding_affinity) %in% c(
    "pdb_ids",
    "comp_id",
    "value"
  )]
  should_be_here <- c(
    "pdb_ids",
    "comp_id",
    "value"
  )[!c(
    "pdb_ids",
    "comp_id",
    "value"
  ) %in% colnames(rcsb_binding_affinity)]

  rcsb_binding_affinity <- rcsb_binding_affinity %>%
    dplyr::select(-all_of(should_not_be_here)) %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here),
        nrow = nrow(rcsb_binding_affinity)
      )),
      should_be_here
    )) %>%
    dplyr::rename(
      affinity_comp_id = "comp_id",
      affinity_value = "value"
    )

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("[3/6] Extracting polymer information: ")
    message("-> 1/7 UniProt IDs ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  polymer_entities <- query_result_clean %>%
    dplyr::select("pdb_ids", "entries.polymer_entities") %>%
    tidyr::unnest("entries.polymer_entities") %>%
    dplyr::bind_cols(
      .$entity_poly,
      .$rcsb_polymer_entity_container_identifiers
    ) %>%
    dplyr::select(-c(
      "entity_poly",
      "rcsb_polymer_entity_container_identifiers",
      "rcsb_entity_source_organism"
    )) %>%
    tidyr::unnest("rcsb_polymer_entity") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(rcsb_non_std_monomers = ifelse(!is.null(unlist(.data$rcsb_non_std_monomers)),
      paste0(.data$rcsb_non_std_monomers, collapse = ";"),
      NA
    )) %>%
    dplyr::ungroup()

  # Deal with cases in which UniProt information of some entries is available but not for others
  polymer_entities_no_uniprots <- polymer_entities %>%
    dplyr::rowwise() %>%
    dplyr::mutate(no_uniprots = is.null(unlist(.data$uniprots))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$no_uniprots) %>%
    dplyr::select(-c("uniprots", "no_uniprots"))

  if (nrow(polymer_entities_no_uniprots) > 0) {
    polymer_entities <- polymer_entities %>%
      tidyr::unnest(c("uniprots")) %>%
      dplyr::bind_rows(polymer_entities_no_uniprots)
  } else {
    polymer_entities <- polymer_entities %>%
      tidyr::unnest(c("uniprots"))
  }

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 2/7 UniProt alignment ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  # Deal with cases in which polymer entity alignment information of some entries is available but not for others
  polymer_entities_no_rcsb_polymer_entity_align <- polymer_entities %>%
    dplyr::rowwise() %>%
    dplyr::mutate(no_rcsb_polymer_entity_align = is.null(unlist(.data$rcsb_polymer_entity_align))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$no_rcsb_polymer_entity_align) %>%
    dplyr::select(-c("rcsb_polymer_entity_align", "no_rcsb_polymer_entity_align"))

  if (nrow(polymer_entities_no_rcsb_polymer_entity_align) > 0) {
    polymer_entities <- polymer_entities %>%
      tidyr::unnest(c("rcsb_polymer_entity_align")) %>%
      dplyr::bind_rows(polymer_entities_no_rcsb_polymer_entity_align)
  } else {
    polymer_entities <- polymer_entities %>%
      tidyr::unnest(c("rcsb_polymer_entity_align"))
  }

  # some proteins do not contain UniProt information therefore data needs to be extracted differently
  if ("rcsb_uniprot_container_identifiers" %in% colnames(polymer_entities)) {
    polymer_entities <- polymer_entities %>%
      dplyr::bind_cols(
        uniprot_container_identifiers = .$rcsb_uniprot_container_identifiers,
        uniprot_protein = .$rcsb_uniprot_protein
      ) %>%
      dplyr::select(-c("rcsb_uniprot_container_identifiers", "rcsb_uniprot_protein"))
  } else {
    polymer_entities <- polymer_entities %>%
      dplyr::select(-c("uniprots")) %>%
      dplyr::mutate(
        uniprot_id = NA,
        name = data.frame(value = NA)
      )
  }

  # The alignment can also be missing independent of the UniProt information
  if ("aligned_regions" %in% colnames(polymer_entities)) {
    # if there are entries that do not have aligned regions these need to be processed separately to not lose them to an unnest
    polymer_entities_no_aligned_regions <- polymer_entities %>%
      dplyr::rowwise() %>%
      dplyr::mutate(no_aligned_regions = is.null(unlist(.data$aligned_regions))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.data$no_aligned_regions)

    polymer_entities <- polymer_entities %>%
      tidyr::unnest(c("aligned_regions")) %>%
      dplyr::bind_cols(.$name) %>%
      dplyr::select(-c("name", "entry_id")) %>%
      dplyr::rename(name_protein = "value") %>%
      tidyr::unnest(c("auth_asym_ids", "polymer_entity_instances")) %>%
      dplyr::bind_cols(
        rcsb_polymer_entity_instance_container_identifiers = .$rcsb_polymer_entity_instance_container_identifiers
      ) %>%
      dplyr::select(-c("rcsb_polymer_entity_instance_container_identifiers"))

    if (nrow(polymer_entities_no_aligned_regions) > 0) {
      polymer_entities_no_aligned_regions <- polymer_entities_no_aligned_regions %>%
        dplyr::select(-c("aligned_regions", "no_aligned_regions")) %>%
        dplyr::bind_cols(.$name) %>%
        dplyr::select(-c("name", "entry_id")) %>%
        dplyr::rename(name_protein = "value") %>%
        tidyr::unnest(c("auth_asym_ids", "polymer_entity_instances")) %>%
        dplyr::bind_cols(
          rcsb_polymer_entity_instance_container_identifiers = .$rcsb_polymer_entity_instance_container_identifiers
        ) %>%
        dplyr::select(-c("rcsb_polymer_entity_instance_container_identifiers")) %>%
        dplyr::mutate(
          entity_beg_seq_id = NA,
          ref_beg_seq_id = NA,
          length = NA
        )

      # Join columns back into main data.frame
      polymer_entities <- polymer_entities %>%
        dplyr::bind_rows(polymer_entities_no_aligned_regions)
    }
  } else {
    polymer_entities <- polymer_entities %>%
      dplyr::select(-c("rcsb_polymer_entity_align")) %>%
      dplyr::bind_cols(.$name) %>%
      dplyr::select(-c("name", "entry_id")) %>%
      dplyr::rename(name_protein = "value") %>%
      tidyr::unnest(c("auth_asym_ids", "polymer_entity_instances")) %>%
      dplyr::bind_cols(
        rcsb_polymer_entity_instance_container_identifiers = .$rcsb_polymer_entity_instance_container_identifiers
      ) %>%
      dplyr::select(-c("rcsb_polymer_entity_instance_container_identifiers")) %>%
      dplyr::mutate(
        entity_beg_seq_id = NA,
        ref_beg_seq_id = NA,
        length = NA,
        reference_database_accession = NA,
        reference_database_isoform = NA,
        reference_database_name = NA
      )
  }

  # Extract ligand information
  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 3/7 Ligand binding sites ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  polymer_entities_no_ligands <- polymer_entities %>%
    dplyr::rowwise() %>%
    dplyr::mutate(no_ligands = is.null(unlist(.data$rcsb_ligand_neighbors))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$no_ligands) %>%
    dplyr::select(-c("rcsb_ligand_neighbors", "no_ligands")) %>%
    dplyr::mutate(
      atom_id = as.character(NA),
      auth_seq_id = as.integer(NA),
      comp_id = as.character(NA),
      ligand_asym_id = as.character(NA),
      ligand_atom_id = as.character(NA),
      ligand_comp_id = as.character(NA),
      ligand_entity_id = as.character(NA),
      ligand_is_bound = as.character(NA),
      seq_id = as.integer(NA)
    )

  polymer_entities <- polymer_entities %>%
    tidyr::unnest(c("rcsb_ligand_neighbors")) %>%
    dplyr::bind_rows(polymer_entities_no_ligands) %>%
    dplyr::mutate(ligand_is_bound = ifelse(.data$ligand_is_bound == "Y", "TRUE", "FALSE")) %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id, .data$ligand_entity_id) %>%
    dplyr::mutate(dplyr::across(
      .cols = c(
        "atom_id",
        "auth_seq_id",
        "comp_id",
        "ligand_asym_id",
        "ligand_atom_id",
        "ligand_comp_id",
        "ligand_is_bound",
        "seq_id"
      ),
      .fns = ~ paste0(.x, collapse = ";")
    )) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id) %>%
    dplyr::mutate(dplyr::across(
      .cols = c(
        "atom_id",
        "auth_seq_id",
        "comp_id",
        "ligand_asym_id",
        "ligand_atom_id",
        "ligand_comp_id",
        "ligand_entity_id",
        "ligand_is_bound",
        "seq_id"
      ),
      .fns = ~ paste0(.x, collapse = "|")
    )) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(
      .cols = c(
        "atom_id",
        "auth_seq_id",
        "comp_id",
        "ligand_asym_id",
        "ligand_atom_id",
        "ligand_comp_id",
        "ligand_entity_id",
        "ligand_is_bound",
        "seq_id"
      ),
      .fns = ~ ifelse(str_detect(.data$atom_id, pattern = "NA"), NA, .x)
    )) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      ligand_donor_atom_id = "atom_id",
      ligand_donor_auth_seq_id = "auth_seq_id",
      ligand_donor_id = "comp_id",
      ligand_id = "ligand_comp_id",
      ligand_label_asym_id = "ligand_asym_id",
      ligand_bond_is_covalent_or_coordinating = "ligand_is_bound",
      ligand_donor_label_seq_id = "seq_id"
    )

  if ("rcsb_ligand_neighbors" %in% colnames(polymer_entities)) {
    # if none of the retrieved entries contains any ligands then this column needs to be removed manually
    polymer_entities <- polymer_entities %>%
      select(-"rcsb_ligand_neighbors")
  }
  # extract modified monomer information

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 4/7 Modified monomers ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  rcsb_polymer_entity_feature <- polymer_entities %>%
    dplyr::select("pdb_ids", "auth_asym_id", "rcsb_polymer_entity_feature") %>%
    tidyr::unnest("rcsb_polymer_entity_feature") %>%
    tidyr::unnest("feature_positions")

  modified_monomer <- rcsb_polymer_entity_feature %>%
    dplyr::filter(.data$type %in% c("modified_monomer")) %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id) %>%
    dplyr::mutate(residue = paste0(
      .data$beg_comp_id,
      ":",
      .data$beg_seq_id,
      "(",
      stringr::str_replace(.data$name, pattern = "Parent monomer ", replacement = ""),
      ")"
    )) %>%
    dplyr::mutate(modified_monomer = paste0(.data$residue, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$modified_monomer, .data$pdb_ids, .data$auth_asym_id)

  # extract secondary structure information
  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 5/7 Secondary structure ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  rcsb_polymer_instance_feature_data <- polymer_entities %>%
    dplyr::select("pdb_ids", "auth_asym_id", "rcsb_polymer_instance_feature") %>%
    tidyr::unnest("rcsb_polymer_instance_feature") %>%
    tidyr::unnest("feature_positions")

  secondary_structures <- rcsb_polymer_instance_feature_data %>%
    dplyr::filter(.data$name %in% c("helix", "sheet")) %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id) %>%
    dplyr::mutate(from_to = paste0(.data$name, ":", .data$beg_seq_id, "-", .data$end_seq_id)) %>%
    dplyr::mutate(secondary_structure = paste0(.data$from_to, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$secondary_structure, .data$pdb_ids, .data$auth_asym_id)

  # extract info about unmodeled residues
  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 6/7 Unmodeled residues ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  unmodeled_residues <- rcsb_polymer_instance_feature_data %>%
    dplyr::filter(.data$name %in% c("partially modeled residue", "unmodeled residue")) %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id) %>%
    dplyr::mutate(from_to = paste0(.data$name, ":", .data$beg_seq_id, "-", .data$end_seq_id)) %>%
    dplyr::mutate(unmodeled_structure = paste0(.data$from_to, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$unmodeled_structure, .data$pdb_ids, .data$auth_asym_id)

  # join secondary structure, unmodeled info and modified monomer

  polymer_entities <- polymer_entities %>%
    left_join(modified_monomer, by = c("pdb_ids", "auth_asym_id")) %>%
    left_join(secondary_structures, by = c("pdb_ids", "auth_asym_id")) %>%
    left_join(unmodeled_residues, by = c("pdb_ids", "auth_asym_id")) %>%
    select(-c("rcsb_polymer_instance_feature", "rcsb_polymer_entity_feature"))

  # extract info about binding-sites
  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("-> 7/7 Ligand Binding ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  # Ligand-binding information is joined to the non-polymer ligand info later on based on the ligands name
  ligand_binding <- rcsb_polymer_instance_feature_data %>%
    dplyr::filter(stringr::str_detect(.data$name, pattern = "ligand")) %>%
    dplyr::mutate(ligand = str_extract(.data$name, pattern = "(?<=ligand ).+")) %>%
    dplyr::group_by(.data$pdb_ids, .data$auth_asym_id, .data$ligand) %>%
    dplyr::mutate(nonpolymer_donor_label_seq_id = paste0(.data$beg_seq_id, collapse = ";")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$nonpolymer_donor_label_seq_id, .data$pdb_ids, .data$auth_asym_id, .data$ligand)

  # Modify auth_seq_id positions that are either duplicated or missing.
  # Missing or duplicated entries are identified by comparing the length of auth_seq_id to the length of the sequence.
  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("[4/6] Correct author sequence positions for some PDB IDs ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  fix_auth_seq_positions <- polymer_entities %>%
    dplyr::rowwise() %>%
    dplyr::filter(length(unlist(.data$auth_to_entity_poly_seq_mapping)) != nchar(.data$pdbx_seq_one_letter_code_can))

  # Remove duplicated positions
  polymer_entities_seq_id_duplications <- fix_auth_seq_positions %>%
    dplyr::filter(length(unlist(.data$auth_to_entity_poly_seq_mapping)) > nchar(.data$pdbx_seq_one_letter_code_can)) %>%
    dplyr::mutate(auth_seq_id = list(unique(unlist(.data$auth_to_entity_poly_seq_mapping)))) %>%
    dplyr::ungroup()

  # Add missing positions
  polymer_entities_seq_id_missing <- fix_auth_seq_positions %>%
    dplyr::filter(length(unlist(.data$auth_to_entity_poly_seq_mapping)) < nchar(.data$pdbx_seq_one_letter_code_can))

  if (nrow(polymer_entities_seq_id_missing) > 0) { # do not run the code below if the data frame is empty, this causes error
    polymer_entities_seq_id_missing <- polymer_entities_seq_id_missing %>%
      dplyr::mutate(auth_seq_id_pdb_numeric = list(stringr::str_replace(.data$auth_to_entity_poly_seq_mapping, pattern = "[:alpha:]", replacement = ""))) %>%
      dplyr::mutate(non_consecutive = list(as.numeric(.data$auth_seq_id_pdb_numeric)[c(diff(as.numeric(.data$auth_seq_id_pdb_numeric)) > 1)])) %>%
      dplyr::mutate(n_missing = list(diff(as.numeric(.data$auth_seq_id_pdb_numeric))[which(diff(as.numeric(.data$auth_seq_id_pdb_numeric)) > 1)] - 1)) %>%
      dplyr::mutate(replacement_values_addition = list(purrr::map(
        .x = .data$n_missing,
        .f = ~ 1:.x
      ))) %>%
      dplyr::mutate(replacement_positions = list(rep(which(.data$auth_seq_id_pdb_numeric %in% as.character(.data$non_consecutive)) + 1, .data$n_missing))) %>%
      dplyr::mutate(auth_seq_id = list(R.utils::insert(.data$auth_to_entity_poly_seq_mapping,
        .data$replacement_positions,
        values = as.character(rep(.data$non_consecutive, .data$n_missing) + unlist(.data$replacement_values_addition))
      ))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(
        "auth_seq_id_pdb_numeric",
        "non_consecutive",
        "n_missing",
        "replacement_values_addition",
        "replacement_positions"
      ))
  }
  # Join corrected entries back
  polymer_entities <- polymer_entities %>%
    dplyr::rowwise() %>%
    dplyr::filter(length(unlist(.data$auth_to_entity_poly_seq_mapping)) == nchar(.data$pdbx_seq_one_letter_code_can) |
      is.null(.data$auth_to_entity_poly_seq_mapping) |
      is.na(.data$pdbx_seq_one_letter_code_can)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(auth_seq_id = .data$auth_to_entity_poly_seq_mapping) %>%
    dplyr::bind_rows(polymer_entities_seq_id_duplications) %>%
    dplyr::bind_rows(polymer_entities_seq_id_missing)

  if (show_progress == TRUE) {
    if (nrow(fix_auth_seq_positions) > 0) {
      message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
      message("Corrected entries: ", paste0(unique(fix_auth_seq_positions$pdb_ids), collapse = ", "))
    } else {
      message("None to correct", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    }
  }

  entity_instance_info <- polymer_entities %>%
    dplyr::distinct(
      .data$asym_id,
      .data$auth_asym_id,
      .data$entry_id,
      .data$auth_to_entity_poly_seq_mapping
    ) %>%
    dplyr::rename(
      auth_asym_ids = "auth_asym_id",
      pdb_ids = "entry_id"
    )

  uniprot_info <- polymer_entities %>%
    dplyr::distinct(.data$uniprot_id, .data$name_protein) %>%
    dplyr::rename(
      reference_database_accession = "uniprot_id",
      protein_name = "name_protein"
    )

  polymer_entities <- polymer_entities %>%
    select(-c(
      "uniprot_id",
      "name_protein",
      "asym_id",
      "auth_asym_id",
      "entry_id",
      "auth_to_entity_poly_seq_mapping"
    )) %>%
    distinct()

  if (show_progress == TRUE) {
    message("[5/6] Extract non-polymer information ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  if (!all(is.na(query_result_clean$entries.nonpolymer_entities))) {
    nonpolymer_entities <- query_result_clean %>%
      dplyr::select("pdb_ids", "entries.nonpolymer_entities") %>%
      tidyr::unnest("entries.nonpolymer_entities") %>%
      dplyr::bind_cols(
        .$rcsb_nonpolymer_entity_container_identifiers,
        .$nonpolymer_comp
      ) %>%
      dplyr::bind_cols(.$chem_comp) %>%
      dplyr::select(-c(
        "nonpolymer_comp",
        "rcsb_nonpolymer_entity_container_identifiers",
        "chem_comp",
        "entry_id"
      )) %>%
      tidyr::unnest("auth_asym_ids") %>%
      dplyr::rename(
        name_nonpolymer = "name",
        formula_nonpolymer = "formula",
        formula_weight_nonpolymer = "formula_weight",
        type_nonpolymer = "type",
        id_nonpolymer = "id"
      ) %>% 
      # Fix "SODIUM ION" IDs
      dplyr::mutate(id_nonpolymer = ifelse(.data$name_nonpolymer == "SODIUM ION", "NA", .data$id_nonpolymer)) 
    
    # Create an annotation data frame for ligand_binding data
    np_annotation <- nonpolymer_entities %>% 
      dplyr::select(-c("pdb_ids", "auth_asym_ids")) %>% 
      dplyr::distinct()
    
    # annotate ligand_binding
    ligand_binding <- ligand_binding %>% 
      left_join(np_annotation, by = c("ligand" = "id_nonpolymer"))
    
    nonpolymer_entities <- nonpolymer_entities %>%
      dplyr::full_join(ligand_binding, by = c("id_nonpolymer" = "ligand", 
                                              "auth_asym_ids" = "auth_asym_id", 
                                              "pdb_ids",
                                              "formula_nonpolymer",
                                              "formula_weight_nonpolymer",
                                              "name_nonpolymer",
                                              "type_nonpolymer"))
      
  } else {
    nonpolymer_entities <- polymer_entities %>%
      dplyr::select("pdb_ids", "auth_asym_ids") %>%
      dplyr::mutate(
        name_nonpolymer = NA,
        formula_nonpolymer = NA,
        formula_weight_nonpolymer = NA,
        type_nonpolymer = NA,
        id_nonpolymer = NA,
        nonpolymer_donor_label_seq_id = NA
      )
  }

  additional_info <- query_result_clean %>%
    dplyr::select(-c(
      "entries.nonpolymer_entities",
      "entries.polymer_entities",
      "entries.rcsb_binding_affinity",
      "entries.pdbx_nmr_exptl",
      "entries.pdbx_nmr_exptl_sample_conditions",
      "entries.pdbx_nmr_refine",
      "entries.exptl_crystal_grow",
      "resolution_combined"
    ))

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
    message("[6/6] Combine information ... ", appendLF = FALSE)
    start_time <- Sys.time()
  }

  combined <- polymer_entities %>%
    dplyr::full_join(nonpolymer_entities, by = c("pdb_ids", "auth_asym_ids"), relationship = "many-to-many") %>%
    dplyr::left_join(rcsb_binding_affinity, by = "pdb_ids", relationship = "many-to-many") %>%
    dplyr::left_join(additional_info, by = "pdb_ids") %>%
    dplyr::left_join(crystal_growth_info, by = "pdb_ids") %>%
    dplyr::left_join(nmr_info, by = "pdb_ids") %>%
    dplyr::left_join(resolution_info, by = "pdb_ids") %>%
    dplyr::left_join(uniprot_info, by = "reference_database_accession") %>%
    dplyr::left_join(entity_instance_info, by = c("pdb_ids", "auth_asym_ids")) %>%
    dplyr::rename(
      auth_asym_id = "auth_asym_ids",
      label_asym_id = "asym_id",
      pdb_sequence = "pdbx_seq_one_letter_code_can",
      auth_seq_id_original = "auth_to_entity_poly_seq_mapping",
      engineered_mutation = "pdbx_mutation",
      non_std_monomer = "rcsb_non_std_monomers"
    ) %>%
    dplyr::rowwise() %>%
    # make character string out of list column
    dplyr::mutate(
      auth_seq_id = paste0(.data$auth_seq_id, collapse = ";"),
      auth_seq_id_original = paste0(.data$auth_seq_id_original, collapse = ";")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      auth_seq_id = ifelse(.data$auth_seq_id == "",
        NA,
        .data$auth_seq_id
      ),
      auth_seq_id_original = ifelse(.data$auth_seq_id_original == "",
        NA,
        .data$auth_seq_id_original
      )
    ) %>%
    dplyr::select(
      "pdb_ids",
      "auth_asym_id",
      "label_asym_id",
      "reference_database_accession",
      "protein_name",
      "reference_database_name",
      "entity_beg_seq_id",
      "ref_beg_seq_id",
      "length",
      "pdb_sequence",
      "auth_seq_id",
      "auth_seq_id_original",
      "engineered_mutation",
      "modified_monomer",
      "ligand_donor_atom_id",
      "ligand_donor_auth_seq_id",
      "ligand_donor_label_seq_id",
      "ligand_donor_id",
      "ligand_label_asym_id",
      "ligand_atom_id",
      "ligand_id",
      "ligand_entity_id",
      "ligand_bond_is_covalent_or_coordinating",
      "secondary_structure",
      "unmodeled_structure",
      "id_nonpolymer",
      "nonpolymer_donor_label_seq_id",
      "type_nonpolymer",
      "formula_weight_nonpolymer",
      "name_nonpolymer",
      "formula_nonpolymer",
      "experimental_method",
      "structure_method",
      "affinity_comp_id",
      "affinity_value",
      "pdbx_keywords",
      "assembly_count",
      "inter_mol_metalic_bond_count",
      "pH_crystallisation",
      "temp_crystallisation",
      "method_crystallisation",
      "type_nmr",
      "ionic_strength_nmr",
      "pH_nmr",
      "temp_nmr",
      "method_nmr",
      "resolution_combined"
    )

  if (show_progress == TRUE) {
    message("DONE ", paste0("(", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), digits = 2), "s)"))
  }

  combined
}
