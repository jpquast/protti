context("test-auxiliary_functions")

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("read_protti works", {
    test_data <- read_protti(test_path("test_import.csv"))
    expect_is(test_data, "data.frame")
    expect_equal(ncol(test_data), 3)
    expect_equal(nrow(test_data), 3)
    expect_equal(names(test_data), c("test_column", "test_column_2", "test_column_3"))
    expect_equal(unname(unlist(lapply(test_data, class))), c("numeric", "character", "integer"))
  })

  protein <- fetch_uniprot(uniprot_ids = "P36578")
  if (!is.null(protein)) {
    data <- tibble::tibble(
      protein_id = rep("P36578", 3),
      protein_sequence = rep(protein$sequence, 3),
      peptide = c(
        stringr::str_sub(protein$sequence, start = 87, end = 97),
        stringr::str_sub(protein$sequence, start = 59, end = 71),
        stringr::str_sub(protein$sequence, start = 10, end = 18)
      )
    )

    assigned_types <- data %>%
      find_peptide(
        protein_sequence = protein_sequence,
        peptide_sequence = peptide
      ) %>%
      peptide_type(aa_before = aa_before, last_aa = last_aa)

    test_that("find_peptide and peptide_type work", {
      expect_is(assigned_types, "data.frame")
      expect_equal(nrow(assigned_types), 3)
      expect_equal(ncol(assigned_types), 9)
      expect_equal(assigned_types$pep_type, c("fully-tryptic", "semi-tryptic", "non-tryptic"))
    })

    coverage <- sequence_coverage(data = assigned_types, protein_sequence = protein_sequence, peptides = peptide)

    test_that("sequence_coverage works", {
      expect_is(coverage, "data.frame")
      expect_equal(nrow(coverage), 3)
      expect_equal(ncol(coverage), 10)
      expect_equal(unique(round(coverage$coverage, digits = 1)), 7.7)
    })

    plot_data <- coverage %>%
      dplyr::mutate(
        fold_change = c(3, -0.4, 2.1),
        protein_length = nchar(protein_sequence)
      )

    test_that("woods_plot works", {
      p <- woods_plot(
        data = plot_data,
        fold_change = fold_change,
        start_position = start,
        end_position = end,
        protein_length = protein_length,
        coverage = coverage,
        protein_id = protein_id,
        colouring = pep_type
      )
      expect_is(p, "ggplot")
      expect_error(print(p), NA)
    })

    test_that("barcode_plot works", {
      p <- barcode_plot(
        data = plot_data,
        start_position = start,
        end_position = end,
        protein_length = protein_length,
        coverage = coverage,
        protein_id = protein_id,
        colouring = pep_type
      )
      expect_is(p, "ggplot")
      expect_error(print(p), NA)
    })
  }
}

test_that("scale_protti works", {
  scale01 <- round(scale_protti(c(1, 2, 1, 4, 6, 8), method = "01"), digits = 1)
  expect_equal(scale01, c(0.0, 0.1, 0.0, 0.4, 0.7, 1.0))

  scale_center <- round(scale_protti(c(1, 2, 1, 4, 6, 8), method = "center"), digits = 1)
  expect_equal(scale_center, c(-0.9, -0.6, -0.9, 0.1, 0.8, 1.5))
})
