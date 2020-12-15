context("test-qc_functions")

test_that("create_synthetic_data works", {
  set.seed(123)
  data <- create_synthetic_data(n_proteins = 1000,
                                frac_change = 0.05,
                                n_replicates = 3,
                                n_conditions = 2,
                                method = "random_effect")
  expect_equal(nrow(data), 81336)
  expect_equal(ncol(data), 8)
  expect_equal(sum(is.na(data$peptide_intensity_missing)), 10663)
})

test_that("test_function works", {
  test <- test_function(data = data, sample = sample)
  expect_is(test, "data.frame")
  expect_equal(nrow(test), 6)
})

# test_that("qc_data_completeness works", {
#   completeness <- qc_data_completeness(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = FALSE)
#   expect_is(completeness, "data.frame")
#   expect_equal(round(completeness$completeness, digits = 5), c(86.85453, 86.96518, 86.69962, 87.12009, 86.83240, 87.23075))
# })

# test_that("qc_log2_intensity_distribution works", {
#   p <- qc_log2_intensity_distribution(data = data, sample = sample, grouping = peptide, log2_intensity = peptide_intensity_missing)
#   expect_is(p,"ggplot")
#   expect_error(print(p), NA)
# })
# 
# test_that("qc_run_intensity works", {
#   p <- qc_run_intensity(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing)
#   expect_is(p,"ggplot")
#   expect_error(print(p), NA)
# })
# 
# test_that("qc_median_intensities works", {
#   p <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = TRUE)
#   expect_is(p,"plotly")
#   expect_error(print(p), NA)
# 
#   p <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = TRUE, interactive = FALSE)
#   expect_is(p,"ggplot")
#   expect_error(print(p), NA)
# 
#   medians <- qc_median_intensities(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, plot = FALSE)
#   expect_is(ids,"data.frame")
#   expect_equal(round(medians$median_intensity, digits = 2), c(17.19, 17.19, 17.20, 17.20, 17.20, 17.18))
# })
# 
# test_that("qc_ids works", {
#   p <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = TRUE, interactive = TRUE, title = "Test Title")
#   expect_is(p,"plotly")
#   expect_error(print(p), NA)
# 
#   p <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = TRUE, interactive = FALSE, title = "Test Title")
#   expect_is(p,"ggplot")
#   expect_error(print(p), NA)
# 
#   ids <- qc_ids(data = data, sample = sample, grouping = peptide, intensity = peptide_intensity_missing, condition = condition, plot = FALSE)
#   expect_is(ids,"data.frame")
#   expect_equal(ids$count, c(11774, 11789, 11753, 11810, 11771, 11825))
# })
# 
# test_that("qc_cvs works", {
#   data_non_log2 <- data %>%
#     dplyr::mutate(peptide_intensity_missing = 2^peptide_intensity_missing)
# 
#   expect_warning(p <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = TRUE))
#   expect_is(p,"ggplot")
#   expect_error(print(p), NA)
# 
#   cvs <- qc_cvs(data = data_non_log2, grouping = peptide, condition = condition, intensity = peptide_intensity_missing, plot = FALSE)
#   expect_is(cvs, "data.frame")
#   expect_equal(round(cvs$median_cv, digits = 2), c(5.70, 5.72))
# })
# 
# 
