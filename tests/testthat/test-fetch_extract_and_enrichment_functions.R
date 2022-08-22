context("test-fetch_extract_and_enrichment_functions")

# make tests conditional. They onlyr run if environmental variable "TEST_PROTTI" is set to TRUE
if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("fetch_uniprot works", {
    unis <- c("iRT", "P36578", "O43324", "Q00796", "P0CX31;P0CX32", "P00163;P03873;P03879", "P06873_1-100")
    expect_warning(uniprot <- fetch_uniprot(unis))
    expect_is(uniprot, "data.frame")
    expect_equal(nrow(uniprot), 9)
    expect_equal(ncol(uniprot), 17)
  })

  proteome <- fetch_uniprot_proteome(organism_id = "83333", columns = c("accession", "go_f", "xref_string"))
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

  database <- fetch_chebi(stars = c(2, 3))
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

  pdb_ids <- c("6HG1", "1E9I", "6D3Q", "4JHW")
  test_that("fetch_pdb works", {
    pdb <- fetch_pdb(pdb_ids)
    expect_is(pdb, "data.frame")
    expect_equal(nrow(pdb), 36)
    expect_equal(ncol(pdb), 46)
  })

  test_that("fetch_pdb_structure works", {
    pdb_structure <- fetch_pdb_structure(pdb_ids, return_data_frame = TRUE)
    expect_is(pdb_structure, "data.frame")
    expect_equal(nrow(pdb_structure), 45731)
    expect_equal(ncol(pdb_structure), 18)
  })

  test_that("fetch_alphafold_prediction works", {
    af_prediction <- fetch_alphafold_prediction(uniprot_ids = c("F4HVG8", "O15552"), return_data_frame = TRUE)
    expect_is(af_prediction, "data.frame")
    expect_equal(nrow(af_prediction), 7310)
    expect_equal(ncol(af_prediction), 15)
  })

  if (.Platform$OS.type != "windows") {
    test_that("fetch_alphafold_prediction organism fetching works", {
      af_prediction_organism <- fetch_alphafold_prediction(organism_name = "Helicobacter pylori", return_data_frame = FALSE)
      expect_is(af_prediction_organism, "list")
      expect_equal(length(af_prediction_organism), 1538)
      expect_equal(ncol(af_prediction_organism[["O24860"]]), 15)
      expect_equal(nrow(af_prediction_organism[["O24860"]]), 746)
    })
  }

  test_that("fetch_metal_pdb works", {
    metal_pdb <- fetch_metal_pdb(id_type = "pdb", id_value = c("1g54"), metal = "Zn")
    expect_is(metal_pdb, "data.frame")
    expect_equal(nrow(metal_pdb), 5)
    expect_equal(ncol(metal_pdb), 25)
  })

  test_that("deprecated kegg_enrichment works", {
    # first fake significances are generated based on the first 10 rows of every group
    kegg_input <- kegg %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$kegg_id, .data$kegg_id) <= 10), TRUE, FALSE)) %>%
      dplyr::group_by(uniprot_id) %>%
      dplyr::mutate(is_significant = rep(.data$is_significant[1], dplyr::n()))

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(kegg_enriched <- kegg_enrichment(
        data = kegg_input,
        protein_id = uniprot_id,
        is_significant = is_significant,
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        plot = FALSE
      ))
    })

    expect_is(kegg_enriched, "data.frame")
    expect_equal(ncol(kegg_enriched), 10)
    expect_gt(nrow(kegg_enriched), 100)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- kegg_enrichment(
        data = kegg_input,
        protein_id = uniprot_id,
        is_significant = is_significant,
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        plot = TRUE,
        plot_cutoff = "adj_pval 0.01"
      ))
    })
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("calculate_kegg_enrichment works", {
    # first fake significances are generated based on the first 10 rows of every group
    kegg_input <- kegg %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$kegg_id, .data$kegg_id) <= 10), TRUE, FALSE)) %>%
      dplyr::group_by(uniprot_id) %>%
      dplyr::mutate(is_significant = rep(.data$is_significant[1], dplyr::n()))

    kegg_enriched <- calculate_kegg_enrichment(
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

    p <- calculate_kegg_enrichment(
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

  test_that("deprecated go_enrichment works", {
    go_input <- proteome %>%
      dplyr::distinct(.data$accession, .data$go_f) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$accession, .data$accession) <= 500), TRUE, FALSE)) %>%
      dplyr::slice(1:3000)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(go_enriched <- go_enrichment(
        data = go_input,
        protein_id = accession,
        is_significant = is_significant,
        ontology_type = "MF",
        organism_id = "83333",
        plot = FALSE
      ))
    })

    expect_is(go_enriched, "data.frame")
    expect_equal(ncol(go_enriched), 9)
    expect_gt(nrow(go_enriched), 1000)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(go_enriched_data <- go_enrichment(
        data = go_input,
        protein_id = accession,
        is_significant = is_significant,
        ontology_type = "MF",
        go_data = go_eco,
        plot = FALSE
      ))
    })

    expect_is(go_enriched_data, "data.frame")
    expect_equal(ncol(go_enriched_data), 9)
    expect_gt(nrow(go_enriched_data), 1000)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(go_enriched_uniprot <- go_enrichment(
        data = go_input,
        protein_id = accession,
        is_significant = is_significant,
        go_annotations_uniprot = go_f,
        plot = FALSE
      ))
    })

    expect_is(go_enriched_uniprot, "data.frame")
    expect_equal(ncol(go_enriched_uniprot), 10)
    expect_gt(nrow(go_enriched_uniprot), 1000)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- go_enrichment(
        data = go_input,
        protein_id = accession,
        is_significant = is_significant,
        ontology_type = "MF",
        organism_id = "83333",
        plot = TRUE,
        plot_cutoff = "adj_pval 0.05"
      ))
    })
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("calculate_go_enrichment works", {
    go_input <- proteome %>%
      dplyr::distinct(.data$accession, .data$go_f) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$accession, .data$accession) <= 500), TRUE, FALSE)) %>%
      dplyr::slice(1:3000)

    go_enriched <- calculate_go_enrichment(
      data = go_input,
      protein_id = accession,
      is_significant = is_significant,
      ontology_type = "MF",
      organism_id = "83333",
      plot = FALSE
    )

    expect_is(go_enriched, "data.frame")
    expect_equal(ncol(go_enriched), 9)
    expect_gt(nrow(go_enriched), 1000)

    go_enriched_data <- calculate_go_enrichment(
      data = go_input,
      protein_id = accession,
      is_significant = is_significant,
      ontology_type = "MF",
      go_data = go_eco,
      plot = FALSE
    )

    expect_is(go_enriched_data, "data.frame")
    expect_equal(ncol(go_enriched_data), 9)
    expect_gt(nrow(go_enriched_data), 1000)

    go_enriched_uniprot <- calculate_go_enrichment(
      data = go_input,
      protein_id = accession,
      is_significant = is_significant,
      go_annotations_uniprot = go_f,
      plot = FALSE
    )

    expect_is(go_enriched_uniprot, "data.frame")
    expect_equal(ncol(go_enriched_uniprot), 10)
    expect_gt(nrow(go_enriched_uniprot), 1000)

    p <- calculate_go_enrichment(
      data = go_input,
      protein_id = accession,
      is_significant = is_significant,
      ontology_type = "MF",
      organism_id = "83333",
      plot = TRUE,
      plot_cutoff = "adj_pval 0.05"
    )
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("deprecated treatment_enrichment works", {
    enrichment_input <- go_eco %>%
      dplyr::distinct(.data$db_id) %>%
      dplyr::mutate(binds_treatment = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE)) %>%
      dplyr::arrange(.data$db_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE))

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(treatment_enriched <- treatment_enrichment(
        data = enrichment_input,
        protein_id = db_id,
        is_significant = is_significant,
        binds_treatment = binds_treatment,
        treatment_name = "test treatment",
        plot = FALSE
      ))
    })
    expect_is(treatment_enriched, "data.frame")
    expect_equal(ncol(treatment_enriched), 4)
    expect_equal(nrow(treatment_enriched), 4)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- treatment_enrichment(
        data = enrichment_input,
        protein_id = db_id,
        is_significant = is_significant,
        binds_treatment = binds_treatment,
        treatment_name = "test treatment",
        plot = TRUE
      ))
    })
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })

  test_that("calculate_treatment_enrichment works", {
    enrichment_input <- go_eco %>%
      dplyr::distinct(.data$db_id) %>%
      dplyr::mutate(binds_treatment = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE)) %>%
      dplyr::arrange(.data$db_id) %>%
      dplyr::mutate(is_significant = ifelse((match(.data$db_id, .data$db_id) <= 500), TRUE, FALSE))

    treatment_enriched <- calculate_treatment_enrichment(
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

    p <- calculate_treatment_enrichment(
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

  test_that("deprecated network_analysis works", {
    # does not check halo_color argument. value of score_threshold is not changed. Only E. coli protein interactions are checked.
    input <- proteome %>%
      dplyr::slice(1:400) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 300)))

    input_many <- proteome %>%
      dplyr::slice(1:1000) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 900)))

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(network <- network_analysis(
        data = input_many,
        protein_id = accession,
        string_id = xref_string,
        organism_id = 511145,
        score_threshold = 900,
        binds_treatment = is_known,
        plot = FALSE
      ))
    })
    expect_is(network, "data.frame")
    expect_equal(ncol(network), 5)
    expect_gt(nrow(network), 100)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(expect_error(network_analysis(
        data = input,
        protein_id = accession,
        string_id = xref_string,
        organism_id = 511145,
        score_threshold = 900,
        binds_treatment = is_known,
        plot = TRUE
      ), NA))
    })
  })

  test_that("analyse_functional_network works", {
    # does not check halo_color argument. value of score_threshold is not changed. Only E. coli protein interactions are checked.
    input <- proteome %>%
      dplyr::slice(1:400) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 300)))

    input_many <- proteome %>%
      dplyr::slice(1:1000) %>%
      dplyr::mutate(is_known = c(rep(TRUE, 100), rep(FALSE, 900)))

    network <- analyse_functional_network(
      data = input_many,
      protein_id = accession,
      string_id = xref_string,
      organism_id = 511145,
      score_threshold = 900,
      binds_treatment = is_known,
      plot = FALSE
    )
    expect_is(network, "data.frame")
    expect_equal(ncol(network), 5)
    expect_gt(nrow(network), 100)

    expect_error(analyse_functional_network(
      data = input,
      protein_id = accession,
      string_id = xref_string,
      organism_id = 511145,
      score_threshold = 900,
      binds_treatment = is_known,
      plot = TRUE
    ), NA)
  })

  eco <- fetch_eco()
  eco_relation <- fetch_eco(return_relation = TRUE)
  test_that("fetch_eco works", {
    expect_is(eco, "data.frame")
    expect_gt(nrow(eco), 2300)
    expect_equal(ncol(eco), 10)

    expect_is(eco_relation, "data.frame")
    expect_gt(nrow(eco_relation), 3200)
    expect_equal(ncol(eco_relation), 3)

    eco_history <- fetch_eco(return_history = TRUE)
    expect_is(eco_history, "data.frame")
    expect_gt(nrow(eco_history), 15500)
    expect_equal(ncol(eco_history), 5)

    expect_error(fetch_eco(return_relation = TRUE, return_history = TRUE))
  })
  
  uniprot_ids <- c("P00393", "P06129", "A0A0C5Q309", "A0A0C9VD04")
  annotations <- fetch_quickgo(type = "annotations", id = uniprot_ids, ontology = "molecular_function")
  test_that("fetch_quickgo works", {
    expect_is(annotations, "data.frame")
    expect_equal(nrow(annotations), 22)
    expect_equal(ncol(annotations), 15)
    
    terms <- fetch_quickgo(type = "terms")
    expect_is(terms, "data.frame")
    expect_equal(nrow(terms), 47906)
    expect_equal(ncol(terms), 13)

    slims <- fetch_quickgo(type = "slims", go_id_slims = c("GO:0046872", "GO:0051540"))
    expect_is(slims, "data.frame")
    expect_equal(nrow(slims), 43)
    expect_equal(ncol(slims), 2)
    
    expect_warning(fetch_quickgo(type = "annotations", 
                               id = c("P63328","Q4FFP4"), 
                               ontology = "molecular_function", 
                               go_id_slims = c("GO:0046872", "GO:0051540")))
  })
  
  test_that("extract_metal_binders works", {
    data_uniprot <- fetch_uniprot(uniprot_ids,
                                  columns = c(
                                    "ft_binding",
                                    "cc_cofactor",
                                    "cc_catalytic_activity"
                                  ))
    
    metal_info <- extract_metal_binders(data_uniprot = data_uniprot,
                                        data_quickgo = annotations,
                                        data_chebi = database,
                                        data_chebi_relation = relations,
                                        data_eco = eco,
                                        data_eco_relation = eco_relation)
    
    expect_is(metal_info, "data.frame")
    expect_equal(ncol(metal_info), 20)
    expect_equal(nrow(metal_info), 34)
  })
}