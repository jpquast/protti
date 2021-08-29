context("test-structure_functions")

if (Sys.getenv("TEST_PROTTI") == "true") {
  peptide_data <- tibble::tibble(
    uniprot_id = c("P0A8T7", "P0A8T7", "P60906"),
    peptide_sequence = c("SGIVSFGKETKGKRRLVITPVDGSDPYEEMIPKWRQLNV", "NVFEGERVER", "AIGEVTDVVEKE"),
    start = c(1160, 1197, 55),
    end = c(1198, 1206, 66),
    map_value = c(70, 100, 100)
  )

  positions_structure <- find_peptide_in_structure(peptide_data,
    peptide = peptide_sequence,
    start = start,
    end = end,
    uniprot_id = uniprot_id,
    retain_columns = c(map_value)
  )

  test_that("find_peptide_in_structure works", {
    expect_is(positions_structure, "data.frame")
    expect_equal(nrow(positions_structure), 312)
    expect_equal(ncol(positions_structure), 16)
  })

  positions_structure_filter <- positions_structure %>%
    filter(pdb_ids %in% c("6UU2", "2EL9"))

  test_that("map_peptide_on_structure works", {
    map_peptides_on_structure(
      peptide_data = positions_structure_filter,
      uniprot_id = uniprot_id,
      pdb_id = pdb_ids,
      chain = auth_asym_id,
      start_in_pdb = auth_seq_id_start,
      end_in_pdb = auth_seq_id_end,
      map_value = map_value,
      file_format = ".cif",
      export_location = tempdir()
    )

    expect_gt(file.info(paste0(tempdir(), "/6UU2_P0A8T7.cif"))$size, 5000000)
    expect_gt(file.info(paste0(tempdir(), "/2EL9_P60906.cif"))$size, 1000000)

    file_cif_6UU2_P0A8T7 <- readr::read_tsv(paste0(tempdir(), "/6UU2_P0A8T7.cif"), col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        b_factor = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+[:space:]+")[[1]][15], # extract b-factor values based on positions
        chain = stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][19],
        residue = suppressWarnings(as.numeric(stringr::str_extract_all(.data$atoms, "[\\w[:punct:]]+")[[1]][17]))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(residue %in% c(1160, 1197, 1206) & chain == "DDD") %>%
      dplyr::distinct(b_factor, residue, chain)

    expect_equal(file_cif_6UU2_P0A8T7$b_factor, c("50  ", "100 ", "100 "))

    expect_message(map_peptides_on_structure(
      peptide_data = positions_structure_filter,
      uniprot_id = uniprot_id,
      pdb_id = pdb_ids,
      chain = auth_asym_id,
      start_in_pdb = auth_seq_id_start,
      end_in_pdb = auth_seq_id_end,
      map_value = map_value,
      file_format = ".pdb",
      export_location = tempdir()
    ))

    expect_false(file.exists(paste0(tempdir(), "/6UU2_P0A8T7.pdb")))
    expect_gt(file.info(paste0(tempdir(), "/2EL9_P60906.pdb"))$size, 1000000)

    file_pdb_6UU2_P0A8T7 <- readr::read_tsv(paste0(tempdir(), "/2EL9_P60906.pdb"), col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::mutate(atoms = ifelse(stringr::str_detect(.data$X1, pattern = "^ATOM|^HETATM"), .data$X1, NA)) %>%
      dplyr::mutate(
        chain = stringr::str_sub(.data$atoms, start = 22, end = 22),
        residue = suppressWarnings(as.numeric(stringr::str_sub(.data$atoms, start = 23, end = 26))),
        b_factor = stringr::str_sub(.data$atoms, start = 61, end = 66)
      ) %>%
      dplyr::filter(residue %in% c(50, 55) & chain == "A") %>%
      dplyr::distinct(b_factor, residue, chain)

    expect_equal(file_pdb_6UU2_P0A8T7$b_factor, c("     0", "   100"))
  })

  test_that("create_structure_contact_map works", {
    data_input <- tibble::tibble(
      pdb_id = c("6NPF", "1C14"),
      chain = c("A", "A"),
      start = c(1, NA),
      end = c(10, NA)
    )
    contact_maps <- create_structure_contact_map(
      data = data_input,
      id = pdb_id,
      chain = chain,
      start_in_pdb = start,
      end_in_pdb = end,
      return_min_residue_distance = TRUE
    )

    expect_is(contact_maps, "list")
    expect_equal(length(contact_maps), 2)
    expect_equal(ncol(contact_maps[["6NPF"]]), 14)
    expect_equal(nrow(contact_maps[["6NPF"]]), 504)
    expect_equal(nrow(contact_maps[["1C14"]]), 18553)
  })
}
