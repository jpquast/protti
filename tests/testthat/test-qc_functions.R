context("test-qc_functions")

set.seed(123)
data <- create_synthetic_data(
  n_proteins = 100,
  frac_change = 0.05,
  n_replicates = 3,
  n_conditions = 2,
  method = "random_effect"
)

random_proteins <- sample(unique(data$protein), size = 6)

test_that("create_synthetic_data works", {
  expect_equal(nrow(data), 8190)
  expect_equal(ncol(data), 14)
  expect_equal(sum(is.na(data$peptide_intensity_missing)), 633)
})

test_that("qc_data_completeness works", {
  completeness <- qc_data_completeness(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = FALSE)
  expect_is(completeness, "data.frame")
  expect_equal(round(completeness$completeness, digits = 5), c(92.89377, 91.64835, 91.72161, 92.16117, 91.86813, 93.33333))

  p <- qc_data_completeness(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_data_completeness(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }
})

test_that("qc_intensity_distribution works", {
  p_facet <- qc_intensity_distribution(data = data, sample = sample, grouping = peptide, intensity_log2 = peptide_intensity_missing, plot_style = "histogram")
  expect_is(p_facet, "ggplot")
  expect_error(print(p_facet), NA)

  p <- qc_intensity_distribution(data = data, grouping = peptide, intensity_log2 = peptide_intensity_missing, plot_style = "histogram")
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  p_boxplot <- qc_intensity_distribution(data = data, sample = sample, grouping = peptide, intensity_log2 = peptide_intensity_missing, plot_style = "boxplot")
  expect_is(p_boxplot, "ggplot")
  expect_error(print(p_boxplot), NA)

  p_violin <- qc_intensity_distribution(data = data, sample = sample, grouping = peptide, intensity_log2 = peptide_intensity_missing, plot_style = "violin")
  expect_is(p_violin, "ggplot")
  expect_error(print(p_violin), NA)
})

test_that("qc_median_intensities works", {
  medians <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = FALSE)
  expect_is(medians, "data.frame")
  expect_equal(round(medians$median_intensity, digits = 2), c(17.34, 17.33, 17.41, 17.36, 17.40, 17.39))

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }

  p <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)
})

test_that("qc_ids works", {
  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = TRUE, interactive = TRUE, title = "Test Title")
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }

  p <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = TRUE, interactive = FALSE, title = "Test Title")
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  ids <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = FALSE)
  expect_is(ids, "data.frame")
  expect_equal(ids$count, c(1252, 1258, 1254, 1274, 1251, 1268))
})

test_that("qc_cvs works", {
  data_non_log2 <- data %>%
    dplyr::mutate(peptide_intensity_missing = 2^peptide_intensity_missing)

  p_density <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = TRUE)
  expect_is(p_density, "ggplot")
  expect_error(print(p_density), NA)

  p_violin <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = TRUE, plot_style = "violin")
  expect_is(p_violin, "ggplot")
  expect_error(print(p_violin), NA)

  p_boxplot <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = TRUE, plot_style = "boxplot")
  expect_is(p_boxplot, "ggplot")
  expect_error(print(p_boxplot), NA)

  cvs <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = FALSE)
  expect_is(cvs, "data.frame")
  expect_equal(round(cvs$median_cv, digits = 2), c(6.06, 6.07))
})

test_that("qc_pca works", {
  p <- qc_pca(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)
})

test_that("qc_sample_correlation works", {
  p <- qc_sample_correlation(data = data, sample = sample, grouping = peptide, intensity_log2 = peptide_intensity_missing, condition = condition, interactive = FALSE)
  expect_is(p, "pheatmap")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_sample_correlation(data = data, sample = sample, grouping = peptide, intensity_log2 = peptide_intensity_missing, condition = condition, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }
})

test_that("qc_proteome_coverage works", {
  proteome <- tibble::tibble(id = 1:4518)
  data_proteins <- tibble::tibble(
    sample = c(rep("A", 101), rep("B", 1000), rep("C", 1000)),
    protein_id = c(proteome$id[1:100], proteome$id[1:1000], proteome$id[1000:2000])
  )

  coverage <- qc_proteome_coverage(data = data_proteins, sample = sample, protein_id = protein_id, organism_id = "83333", plot = FALSE)
  expect_is(coverage, "data.frame")
  expect_equal(round(coverage$percentage, digits = 0), c(2, 98, 22, 78, 22, 78, 44, 56))

  p <- qc_proteome_coverage(data = data_proteins, sample = sample, protein_id = protein_id, organism_id = "83333", plot = TRUE, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_proteome_coverage(data = data_proteins, sample = sample, protein_id = protein_id, organism_id = "83333", plot = TRUE, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }
})

test_that("qc_sequence_coverage works", {
  p <- qc_sequence_coverage(data = data, protein_identifier = protein, coverage = coverage, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  p_facet <- qc_sequence_coverage(data = data, protein_identifier = protein, coverage = coverage, sample = sample, interactive = FALSE)
  expect_is(p_facet, "ggplot")
  expect_error(print(p_facet), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_sequence_coverage(data = data, protein_identifier = protein, coverage = coverage, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)

    p_interactive_facet <- qc_sequence_coverage(data = data, protein_identifier = protein, coverage = coverage, sample = sample, interactive = TRUE)
    expect_is(p_interactive_facet, "plotly")
    expect_error(print(p_interactive_facet), NA)
  }
})

test_that("qc_charge_states works", {
  data_no_na <- data %>%
    tidyr::drop_na(peptide_intensity_missing)
  data_charge_count <- qc_charge_states(data = data_no_na, sample = sample, grouping = peptide, charge_states = charge, intensity = peptide_intensity_missing, method = "count", plot = FALSE)
  expect_is(data_charge_count, "data.frame")
  expect_equal(round(data_charge_count$charge_per, digits = 0), c(8, 54, 33, 4, 0, 8, 55, 33, 4, 0, 8, 54, 33, 4, 0, 8, 54, 33, 4, 0, 9, 54, 33, 4, 0, 8, 55, 33, 4, 0))

  data_charge_intensity <- qc_charge_states(data = data, sample = sample, grouping = peptide, charge_states = charge, intensity = peptide_intensity_missing, method = "intensity", plot = FALSE)
  expect_is(data_charge_intensity, "data.frame")
  expect_equal(round(data_charge_intensity$charge_per, digits = 0), c(8, 54, 33, 4, 0, 8, 55, 33, 4, 0, 8, 54, 33, 4, 0, 8, 55, 33, 4, 0, 9, 54, 33, 4, 0, 8, 55, 33, 4, 0))

  p_count <- qc_charge_states(data = data_no_na, sample = sample, grouping = peptide, charge_states = charge, intensity = peptide_intensity_missing, method = "count", plot = TRUE)
  expect_is(p_count, "ggplot")
  expect_error(print(p_count), NA)

  p_intensity <- qc_charge_states(data = data, sample = sample, grouping = peptide, charge_states = charge, intensity = peptide_intensity_missing, method = "intensity", plot = TRUE)
  expect_is(p_intensity, "ggplot")
  expect_error(print(p_intensity), NA)
})

test_that("qc_missed_cleavages works", {
  data_no_na <- data %>%
    tidyr::drop_na(peptide_intensity_missing)
  data_mc_count <- qc_missed_cleavages(data = data_no_na, sample = sample, grouping = peptide, missed_cleavages = n_missed_cleavage, intensity = peptide_intensity_missing, method = "count", plot = FALSE)
  expect_is(data_mc_count, "data.frame")
  expect_equal(round(data_mc_count$mc_percent, digits = 1), c(76.8, 20.9, 2.3, 76.9, 20.9, 2.2, 77.2, 20.4, 2.4, 77.0, 20.5, 2.5, 76.8, 21.0, 2.2, 77.2, 20.6, 2.3))

  data_mc_intensity <- qc_missed_cleavages(data = data, sample = sample, grouping = peptide, missed_cleavages = n_missed_cleavage, intensity = peptide_intensity_missing, method = "intensity", plot = FALSE)
  expect_is(data_mc_intensity, "data.frame")
  expect_equal(round(data_mc_intensity$mc_percent, digits = 1), c(76.9, 20.8, 2.3, 76.9, 20.9, 2.2, 77.3, 20.4, 2.4, 77.1, 20.5, 2.4, 76.9, 20.9, 2.2, 77.2, 20.5, 2.3))

  p_count <- qc_missed_cleavages(data = data, sample = sample, grouping = peptide, missed_cleavages = n_missed_cleavage, intensity = peptide_intensity_missing, method = "count", plot = TRUE)
  expect_is(p_count, "ggplot")
  expect_error(print(p_count), NA)

  p_intensity <- qc_missed_cleavages(data = data, sample = sample, grouping = peptide, missed_cleavages = n_missed_cleavage, intensity = peptide_intensity_missing, method = "intensity", plot = TRUE)
  expect_is(p_intensity, "ggplot")
  expect_error(print(p_intensity), NA)
})

test_that("qc_contaminants works", {
  contaminant_data <- data %>%
    dplyr::mutate(contaminant = ifelse(protein %in% random_proteins, TRUE, FALSE))

  contaminants <- qc_contaminants(data = contaminant_data, sample = sample, protein = protein, is_contaminant = contaminant, intensity = peptide_intensity_missing, plot = FALSE)
  expect_is(contaminants, "data.frame")
  expect_equal(round(contaminants$contaminant_percentage, digits = 1)[1:12], c(0.2, 0.3, 1.0, 1.2, 1.3, 2.3, 0.2, 0.3, 1.1, 1.2, 1.3, 2.1))
  expect_equal(as.character(unique(contaminants$protein)), c("other", "protein_75", "protein_62", "protein_73", "protein_37", "protein_83", "protein_55"))

  p <- qc_contaminants(data = contaminant_data, sample = sample, protein = protein, is_contaminant = contaminant, intensity = peptide_intensity_missing, plot = TRUE, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_contaminants(data = contaminant_data, sample = sample, protein = protein, is_contaminant = contaminant, intensity = peptide_intensity_missing, plot = TRUE, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }
})

test_that("qc_peptide_type works", {
  data_no_na <- data %>%
    tidyr::drop_na(peptide_intensity_missing)
  peptide_types_count <- qc_peptide_type(data = data_no_na, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "count", plot = FALSE)
  expect_is(peptide_types_count, "data.frame")
  expect_equal(round(peptide_types_count$peptide_type_percent, digits = 1), c(55.8, 5.9, 38.2, 56.1, 6.2, 37.7, 56.2, 5.8, 38.0, 56.4, 5.9, 37.8, 56.0, 5.9, 38.1, 55.8, 6.2, 38.0))

  peptide_types_intenisty <- qc_peptide_type(data = data, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "intensity", plot = FALSE)
  expect_is(peptide_types_intenisty, "data.frame")
  expect_equal(round(peptide_types_intenisty$peptide_type_percent, digits = 1), c(38.0, 37.9, 38.2, 38.0, 6.0, 5.8, 5.8, 5.8, 6.1, 55.9, 56.1, 56.2, 56.3, 56.0, 55.9, 38.2, 37.8, 5.9))

  p_count <- qc_peptide_type(data = data_no_na, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "count", plot = TRUE, interactive = FALSE)
  expect_is(p_count, "ggplot")
  expect_error(print(p_count), NA)

  p_intensity <- qc_peptide_type(data = data, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "intensity", plot = TRUE, interactive = FALSE)
  expect_is(p_intensity, "ggplot")
  expect_error(print(p_intensity), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_count_interactive <- qc_peptide_type(data = data_no_na, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "count", plot = TRUE, interactive = TRUE)
    expect_is(p_count_interactive, "plotly")
    expect_error(print(p_count_interactive), NA)

    p_intensity_interactive <- qc_peptide_type(data = data_no_na, sample = sample, peptide = peptide, pep_type = pep_type, intensity = peptide_intensity_missing, method = "intensity", plot = TRUE, interactive = TRUE)
    expect_is(p_intensity_interactive, "plotly")
    expect_error(print(p_intensity_interactive), NA)
  }
})

test_that("qc_peak_width works", {
  # not testing retention_time_start and retention_time_end arguments
  data_no_na <- data %>%
    tidyr::drop_na(peptide_intensity_missing)
  p <- qc_peak_width(data_no_na, sample = sample, intensity = peptide_intensity_missing, retention_time = retention_time, peak_width = peak_width, interactive = FALSE)
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == TRUE) {
    p_interactive <- qc_peak_width(data_no_na, sample = sample, intensity = peptide_intensity_missing, retention_time = retention_time, peak_width = peak_width, interactive = TRUE)
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }
})
