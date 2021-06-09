context("test-fetch_extract_and_enrichment_functions")

# make tests conditional. They onlyr run if environmental variable "TEST_PROTTI" is set to TRUE
if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("fetch_uniprot works", {
    unis <- c("iRT", "P36578", "O43324", "Q00796", "P0CX31;P0CX32", "P00163;P03873;P03879", "P06873_1-100")
    expect_warning(uniprot <- fetch_uniprot(unis))
    expect_is(uniprot, "data.frame")
    expect_equal(nrow(uniprot), 9)
    expect_equal(ncol(uniprot), 18)
  })

  proteome <- fetch_uniprot_proteome(organism_id = "83333", columns = c("id", "go(molecular function)", "database(String)"))
  test_that("fetch_uniprot_proteome works", {
    expect_is(proteome, "data.frame")
    expect_equal(ncol(proteome), 3)
    expect_gt(nrow(proteome), 10)
  })

  test_that("fetch_mobidb works", {
    unis <- c("iRT", "P25437", "P30870", "P0A6P9")
    expect_warning(mobidb <- fetch_mobidb("83333", unis))
    expect_is(mobidb, "data.frame")
    expect_equal(nrow(mobidb), 2)
    expect_equal(ncol(mobidb), 7)
  })

  database <- fetch_chebi()
  relations <- fetch_chebi(relation = TRUE)
  test_that("fetch_chebi works", {
    expect_is(database, "data.frame")
    expect_is(relations, "data.frame")
    expect_equal(ncol(database), 13)
    expect_equal(ncol(relations), 3)
    expect_gt(nrow(database), 10)
    expect_gt(nrow(relations), 10)
  })

  kegg <- fetch_kegg(species = "eco")
  test_that("fetch_kegg works", {
    expect_is(kegg, "data.frame")
    expect_equal(ncol(kegg), 4)
    expect_gt(nrow(kegg), 10)
  })

  go_eco <- fetch_go(organism_id = "83333")

  test_that("fetch_go works", {
    go_sac <- fetch_go(organism_id = "559292")
    go_hs <- fetch_go(organism_id = "9606")
    expect_is(go_eco, "data.frame")
    expect_equal(ncol(go_eco), 17)
    expect_gt(nrow(go_eco), 10)
    expect_is(go_sac, "data.frame")
    expect_equal(ncol(go_sac), 17)
    expect_gt(nrow(go_sac), 10)
    expect_is(go_hs, "data.frame")
    expect_equal(ncol(go_hs), 17)
    expect_gt(nrow(go_hs), 10)
  })

  test_that("fetch_pdb works", {
    pdb_ids <- c("6HG1", "1E9I", "6D3Q", "4JHW")
    pdb <- fetch_pdb(pdb_ids)
    expect_is(pdb, "data.frame")
    expect_equal(nrow(pdb), 34)
    expect_equal(ncol(pdb), 30)
  })

  test_that("extract_metal_binders works", {
    data_uniprot <- fetch_uniprot(c("Q03640", "Q03778", "P22276"))
    metal_info <- extract_metal_binders(data = data_uniprot, chebi_data = database, chebi_relation_data = relations)

    expect_is(metal_info, "data.frame")
    expect_equal(ncol(metal_info), 9)
    expect_gt(nrow(metal_info), 40)
  })

  test_that("kegg_enrichment works", {
    # first fake significances are generated based on the first 10 rows of every group
    kegg_input <- kegg %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$kegg_id, .data$kegg_id) <= 10), TRUE, FALSE)) %>%
      dplyr::group_by(uniprot_id) %>%
      dplyr::mutate(is_significant = rep(.data$is_significant[1], dplyr::n()))

    kegg_enriched <- kegg_enrichment(
      data = kegg_input,
      protein_id = uniprot_id,
      is_significant = is_significant,
      pathway_id = pathway_id,
      pathway_name = pathway_name,
      plot = FALSE
    )

    expect_is(kegg_enriched, "data.frame")
    expect_equal(ncol(kegg_enriched), 10)
    expect_gt(nrow(kegg_enriched), 100)

    p <- kegg_enrichment(
      data = kegg_input,
      protein_id = uniprot_id,
      is_significant = is_significant,
      pathway_id = pathway_id,
      pathway_name = pathway_name,
      plot = TRUE,
      plot_cutoff = "adj_pval 0.01"
    )
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("go_enrichment works", {
    go_input <- proteome %>%
      dplyr::distinct(.data$id, .data$go_molecular_function) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$id, .data$id) <= 500), TRUE, FALSE)) %>%
      dplyr::slice(1:3000)

    go_enriched <- go_enrichment(
      data = go_input,
      protein_id = id,
      is_significant = is_significant,
      ontology_type = "MF",
      organism_id = "83333",
      plot = FALSE
    )

    expect_is(go_enriched, "data.frame")
    expect_equal(ncol(go_enriched), 9)
    expect_gt(nrow(go_enriched), 1000)

    go_enriched_data <- go_enrichment(
      data = go_input,
      protein_id = id,
      is_significant = is_significant,
      ontology_type = "MF",
      go_data = go_eco,
      plot = FALSE
    )

    expect_is(go_enriched_data, "data.frame")
    expect_equal(ncol(go_enriched_data), 9)
    expect_gt(nrow(go_enriched_data), 1000)

    go_enriched_uniprot <- go_enrichment(
      data = go_input,
      protein_id = id,
      is_significant = is_significant,
      go_annotations_uniprot = go_molecular_function,
      plot = FALSE
    )

    expect_is(go_enriched_uniprot, "data.frame")
    expect_equal(ncol(go_enriched_uniprot), 10)
    expect_gt(nrow(go_enriched_uniprot), 1000)

    p <- go_enrichment(
      data = go_input,
      protein_id = id,
      is_significant = is_significant,
      ontology_type = "MF",
      organism_id = "83333",
      plot = TRUE,
      plot_cutoff = "adj_pval 0.05"
    )
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("treatment_enrichment works", {
    enrichment_input <- go_eco %>%
      dplyr::distinct(.data$db_id) %>%
      dplyr::mutate(binds_treatment = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE)) %>%
      dplyr::arrange(.data$db_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE))

    treatment_enriched <- treatment_enrichment(
      data = enrichment_input,
      protein_id = db_id,
      is_significant = is_significant,
      binds_treatment = binds_treatment,
      treatment_name = "test treatment",
      plot = FALSE
    )
    expect_is(treatment_enriched, "data.frame")
    expect_equal(ncol(treatment_enriched), 4)
    expect_equal(nrow(treatment_enriched), 4)

    p <- treatment_enrichment(
      data = enrichment_input,
      protein_id = db_id,
      is_significant = is_significant,
      binds_treatment = binds_treatment,
      treatment_name = "test treatment",
      plot = TRUE
    )
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("network_analysis works", {
    # does not check halo_color argument. value of score_threshold is not changed. Only E. coli protein interactions are checked.
    input <- proteome %>%
      dplyr::slice(1:400) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 300)))

    input_many <- proteome %>%
      dplyr::slice(1:1000) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 900)))

    network <- network_analysis(
      data = input_many,
      protein_id = id,
      string_id = database_string,
      organism_id = 511145,
      score_threshold = 900,
      binds_treatment = is_known,
      plot = FALSE
    )
    expect_is(network, "data.frame")
    expect_equal(ncol(network), 5)
    expect_gt(nrow(network), 100)

    expect_error(network_analysis(
      data = input,
      protein_id = id,
      string_id = database_string,
      organism_id = 511145,
      score_threshold = 900,
      binds_treatment = is_known,
      plot = TRUE
    ), NA)
  })
}
