#' Fetch structure information from RCSB
#'
#' Fetches structure metadata from RCSB. If you want to retrieve atom data such as positions, use the function \code{fetch_pdb_structure()}.
#'
#' @param pdb_ids a character vector of PDB identifiers.
#' @param batchsize numeric, specifying the number of structures to be processed in a single query. Default is 2000.
#' @param show_progress logical, if true, a progress bar will be shown. Default is TRUE.
#'
#' @return A data frame that contains all structure metadata for the PDB IDs provided.
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom stringr str_replace_all
#' @importFrom curl has_internet
#' @importFrom utils URLencode
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_pdb(c("6HG1", "1E9I", "6D3Q", "4JHW")))
#' }
fetch_pdb <- function(pdb_ids, batchsize = 200, show_progress = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
  }

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
      	rcsb_polymer_entity_instance_container_identifiers {
          asym_id
          auth_asym_id
          entry_id
      		auth_to_entity_poly_seq_mapping
    		}
      }
      entity_poly {
        pdbx_seq_one_letter_code_can
        rcsb_artifact_monomer_count
        rcsb_conflict_count
        rcsb_deletion_count
        rcsb_insertion_count
        rcsb_mutation_count
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

  # split pdb_ids into batches
  batches <- split(pdb_ids, ceiling(seq_along(pdb_ids) / batchsize))
  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(batches))
  }

  # query information from database

  query_result <- purrr::map_dfr(batches, function(x) {
    pdb_ids <- paste0(x, collapse = '", "')
    full_query <- stringr::str_replace_all(query, pattern = "pdb_ids", replacement = pdb_ids)
    url_encode_query <- utils::URLencode(full_query) %>%
      stringr::str_replace_all(pattern = "\\[", replacement = "%5B") %>%
      stringr::str_replace_all(pattern = "\\]", replacement = "%5D")
    # only try to fetch more batches if previous cycle did not encounter a connection problem.
    if (!is.null(batches)) {
      query <- try_query(httr::modify_url("https://data.rcsb.org/graphql", query = url_encode_query), type = "application/json", simplifyDataFrame = TRUE)
    }
    if (show_progress == TRUE & "list" %in% class(query)) {
      pb$tick()
    }
    # if previous batch had a connection problem change batches to NULL, which breaks the mapping.
    if (!"list" %in% class(query)) {
      batches <<- NULL
    }
    # only proceed with data if it was correctly retrieved
    if ("list" %in% class(query)) {
      query %>%
        purrr::flatten() %>%
        as.data.frame(stringsAsFactors = FALSE)
    }
  })

  # process information from database

  query_result_clean <- query_result %>%
    dplyr::bind_cols(
      pdb_ids = .$entries.rcsb_id,
      structure_keywords = .$entries.struct_keywords,
      .$entries.rcsb_entry_info
    ) %>%
    dplyr::select(-c(
      .data$entries.rcsb_id,
      .data$entries.struct_keywords,
      .data$entries.rcsb_entry_info
    )) %>%
    tidyr::unnest(.data$entries.exptl) %>%
    dplyr::rename(structure_method = .data$method)

  crystal_growth_info <- query_result_clean %>%
    dplyr::select(.data$pdb_ids, .data$entries.exptl_crystal_grow) %>%
    tidyr::unnest(.data$entries.exptl_crystal_grow)

  # make sure that the data is complete even if there is no crystal structure
  should_not_be_here <- colnames(crystal_growth_info)[!colnames(crystal_growth_info) %in% c("pdb_ids", "method", "pH", "temp")]
  should_be_here <- c("pdb_ids", "method", "pH", "temp")[!c("pdb_ids", "method", "pH", "temp") %in% colnames(crystal_growth_info)]

  crystal_growth_info <- crystal_growth_info %>%
    dplyr::select(-should_not_be_here) %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here), nrow = nrow(crystal_growth_info))), should_be_here)) %>%
    dplyr::rename(
      pH_crystallisation = .data$pH,
      method_crystallisation = .data$method,
      temp_crystallisation = .data$temp
    )

  resolution_info <- query_result_clean %>%
    dplyr::select(.data$pdb_ids, .data$resolution_combined) %>%
    tidyr::unnest(.data$resolution_combined)

  nmr_info <- query_result_clean %>%
    dplyr::select(.data$pdb_ids, .data$entries.pdbx_nmr_exptl, .data$entries.pdbx_nmr_exptl_sample_conditions, .data$entries.pdbx_nmr_refine) %>%
    tidyr::unnest(.data$entries.pdbx_nmr_exptl) %>%
    tidyr::unnest(.data$entries.pdbx_nmr_exptl_sample_conditions) %>%
    tidyr::unnest(.data$entries.pdbx_nmr_refine)

  # make sure that the data is complete even if there is no NMR structure
  should_not_be_here <- colnames(nmr_info)[!colnames(nmr_info) %in% c("pdb_ids", "type", "pH", "temperature", "method", "ionic_strength")]
  should_be_here <- c("pdb_ids", "type", "pH", "temperature", "method", "ionic_strength")[!c("pdb_ids", "type", "pH", "temperature", "method", "ionic_strength") %in% colnames(nmr_info)]

  nmr_info <- nmr_info %>%
    dplyr::select(-should_not_be_here) %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here), nrow = nrow(nmr_info))), should_be_here)) %>%
    dplyr::rename(
      type_nmr = .data$type,
      pH_nmr = .data$pH,
      temp_nmr = .data$temperature,
      method_nmr = .data$method,
      ionic_strength_nmr = .data$ionic_strength
    )

  rcsb_binding_affinity <- query_result_clean %>%
    dplyr::select(.data$pdb_ids, .data$entries.rcsb_binding_affinity) %>%
    tidyr::unnest(.data$entries.rcsb_binding_affinity)

  # make sure that the data is complete even if there is no affinity information
  should_not_be_here <- colnames(rcsb_binding_affinity)[!colnames(rcsb_binding_affinity) %in% c("pdb_ids", "comp_id", "value")]
  should_be_here <- c("pdb_ids", "comp_id", "value")[!c("pdb_ids", "comp_id", "value") %in% colnames(rcsb_binding_affinity)]

  rcsb_binding_affinity <- rcsb_binding_affinity %>%
    dplyr::select(-should_not_be_here) %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here), nrow = nrow(rcsb_binding_affinity))), should_be_here)) %>%
    dplyr::rename(
      affinity_comp_id = .data$comp_id,
      affinity_value = .data$value
    )

  polymer_entities <- query_result_clean %>%
    dplyr::select(.data$pdb_ids, .data$entries.polymer_entities) %>%
    tidyr::unnest(.data$entries.polymer_entities) %>%
    dplyr::bind_cols(
      .$entity_poly,
      .$rcsb_polymer_entity_container_identifiers
    ) %>%
    dplyr::select(-c(.data$entity_poly, .data$rcsb_polymer_entity_container_identifiers, .data$rcsb_entity_source_organism)) %>%
    tidyr::unnest(c(.data$uniprots, .data$rcsb_polymer_entity_align)) %>%
    dplyr::bind_cols(
      uniprot_container_identifiers = .$rcsb_uniprot_container_identifiers,
      uniprot_protein = .$rcsb_uniprot_protein
    ) %>%
    dplyr::select(-c(.data$rcsb_uniprot_container_identifiers, .data$rcsb_uniprot_protein)) %>%
    tidyr::unnest(c(.data$aligned_regions)) %>%
    dplyr::bind_cols(.$name) %>%
    dplyr::select(-c(.data$name, .data$entry_id)) %>%
    dplyr::rename(name_protein = .data$value) %>%
    tidyr::unnest(c(.data$auth_asym_ids, .data$polymer_entity_instances)) %>% 
    dplyr::bind_cols(
      rcsb_polymer_entity_instance_container_identifiers = .$rcsb_polymer_entity_instance_container_identifiers
    ) %>% 
    dplyr::select(-c(.data$rcsb_polymer_entity_instance_container_identifiers))
  
  entity_instance_info <- polymer_entities %>% 
    dplyr::distinct(.data$asym_id, .data$auth_asym_id, .data$entry_id, .data$auth_to_entity_poly_seq_mapping) %>%
    dplyr::rename(
      auth_asym_ids = .data$auth_asym_id,
      pdb_ids = .data$entry_id
    )

  uniprot_info <- polymer_entities %>%
    dplyr::distinct(.data$uniprot_id, .data$name_protein) %>%
    dplyr::rename(
      reference_database_accession = .data$uniprot_id,
      protein_name = .data$name_protein
    )

  polymer_entities <- polymer_entities %>%
    select(-c(.data$uniprot_id, .data$name_protein, .data$asym_id, .data$auth_asym_id, .data$entry_id, .data$auth_to_entity_poly_seq_mapping)) %>%
    distinct()

  if (!all(is.na(query_result_clean$entries.nonpolymer_entities))) {
    nonpolymer_entities <- query_result_clean %>%
      dplyr::select(.data$pdb_ids, .data$entries.nonpolymer_entities) %>%
      tidyr::unnest(.data$entries.nonpolymer_entities) %>%
      dplyr::bind_cols(
        .$rcsb_nonpolymer_entity_container_identifiers,
        .$nonpolymer_comp
      ) %>%
      dplyr::bind_cols(.$chem_comp) %>%
      dplyr::select(-c(.data$nonpolymer_comp, .data$rcsb_nonpolymer_entity_container_identifiers, .data$chem_comp, .data$entry_id)) %>%
      tidyr::unnest(.data$auth_asym_ids) %>%
      dplyr::rename(
        name_nonpolymer = .data$name,
        formula_nonpolymer = .data$formula,
        formula_weight_nonpolymer = .data$formula_weight,
        type_nonpolymer = .data$type,
        id_nonpolymer = .data$id
      )
  } else {
    nonpolymer_entities <- polymer_entities %>%
      dplyr::select(.data$pdb_ids, .data$auth_asym_ids) %>%
      dplyr::mutate(
        name_nonpolymer = NA,
        formula_nonpolymer = NA,
        formula_weight_nonpolymer = NA,
        type_nonpolymer = NA,
        id_nonpolymer = NA
      )
  }

  additional_info <- query_result_clean %>%
    dplyr::select(-c(
      .data$entries.nonpolymer_entities,
      .data$entries.polymer_entities,
      .data$entries.rcsb_binding_affinity,
      .data$entries.pdbx_nmr_exptl,
      .data$entries.pdbx_nmr_exptl_sample_conditions,
      .data$entries.pdbx_nmr_refine,
      .data$entries.exptl_crystal_grow,
      .data$resolution_combined
    ))

  combined <- polymer_entities %>%
    dplyr::full_join(nonpolymer_entities, by = c("pdb_ids", "auth_asym_ids")) %>%
    dplyr::left_join(rcsb_binding_affinity, by = "pdb_ids") %>%
    dplyr::left_join(additional_info, by = "pdb_ids") %>%
    dplyr::left_join(crystal_growth_info, by = "pdb_ids") %>%
    dplyr::left_join(nmr_info, by = "pdb_ids") %>%
    dplyr::left_join(resolution_info, by = "pdb_ids") %>%
    dplyr::left_join(uniprot_info, by = "reference_database_accession") %>%
    dplyr::left_join(entity_instance_info, by = c("pdb_ids", "auth_asym_ids")) %>% 
    dplyr::rename(
      chain = .data$auth_asym_ids,
      chain_alternative = .data$asym_id,
      pdb_sequence = .data$pdbx_seq_one_letter_code_can
    ) %>%
    dplyr::select(
      .data$pdb_ids,
      .data$chain,
      .data$chain_alternative,
      .data$reference_database_accession,
      .data$protein_name,
      .data$reference_database_name,
      .data$entity_beg_seq_id,
      .data$ref_beg_seq_id,
      .data$length,
      .data$pdb_sequence,
      .data$auth_to_entity_poly_seq_mapping,
      .data$id_nonpolymer,
      .data$type_nonpolymer,
      .data$formula_weight_nonpolymer,
      .data$name_nonpolymer,
      .data$formula_nonpolymer,
      .data$experimental_method,
      .data$structure_method,
      .data$affinity_comp_id,
      .data$affinity_value,
      .data$pdbx_keywords,
      .data$assembly_count,
      .data$inter_mol_metalic_bond_count,
      .data$pH_crystallisation,
      .data$temp_crystallisation,
      .data$method_crystallisation,
      .data$type_nmr,
      .data$ionic_strength_nmr,
      .data$pH_nmr,
      .data$temp_nmr,
      .data$method_nmr,
      .data$resolution_combined
    )

  combined
}
