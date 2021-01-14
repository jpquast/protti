context("test-workflow")

set.seed(123)
data <- create_synthetic_data(
  n_proteins = 500,
  frac_change = 0.05,
  n_replicates = 3,
  n_conditions = 2,
  method = "random_effect",
  additional_metadata = FALSE
)

data_drc <- create_synthetic_data(
  n_proteins = 100,
  frac_change = 0.05,
  n_replicates = 3,
  n_conditions = 8,
  method = "dose_response",
  concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
  additional_metadata = FALSE
)

normalised_data <- data %>%
  median_normalisation(sample = sample, intensity_log2 = peptide_intensity_missing)

normalised_data_drc <- data_drc %>%
  median_normalisation(sample = sample, intensity_log2 = peptide_intensity_missing)

test_that("median_normalization works", {
  non_normalised <- normalised_data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median = median(peptide_intensity_missing, na.rm = TRUE), .groups = "drop")
  normalised <- normalised_data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median = median(normalised_intensity_log2, na.rm = TRUE), .groups = "drop")

  all_equal_non_normalised <- range(non_normalised$median) / mean(non_normalised$median)
  expect_false(isTRUE(all.equal(all_equal_non_normalised[1], all_equal_non_normalised[2]))) # test that medians are unequal before normalizing

  all_equal_normalised <- range(normalised$median) / mean(normalised$median)
  expect_equal(all_equal_normalised[1], all_equal_normalised[2]) # test that medians are equal after normalizing
  
  ## test drc data
  non_normalised_drc <- normalised_data_drc %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median = median(peptide_intensity_missing, na.rm = TRUE), .groups = "drop")
  normalised_drc <- normalised_data_drc %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median = median(normalised_intensity_log2, na.rm = TRUE), .groups = "drop")
  
  all_equal_non_normalised_drc <- range(non_normalised_drc$median) / mean(non_normalised_drc$median)
  expect_false(isTRUE(all.equal(all_equal_non_normalised_drc[1], all_equal_non_normalised_drc[2]))) # test that medians are unequal before normalizing
  
  all_equal_normalised_drc <- range(normalised_drc$median) / mean(normalised_drc$median)
  expect_equal(all_equal_normalised_drc[1], all_equal_normalised_drc[2]) # test that medians are equal after normalizing
})

missing_data <- normalised_data %>%
  assign_missingness(sample = sample, condition = condition, grouping = peptide, intensity = normalised_intensity_log2, ref_condition = "condition_1", retain_columns = c(protein))

missing_data_drc <- normalised_data_drc %>%
  assign_missingness(sample = sample, condition = concentration, grouping = peptide, intensity = normalised_intensity_log2, ref_condition = "0", retain_columns = c(protein))

test_that("assign_missingness works", {
  # not testing noise argument. Also no change of default values for completeness_MAR and completeness_MNAR
  expect_equal(nrow(data), nrow(missing_data))
  expect_true("missingness" %in% colnames(missing_data))

  missingness_count <- missing_data %>%
    dplyr::count(missingness)

  expect_equal(sort(missingness_count$n), c(222, 4134, 5550, 28926))
  
  ## test drc data
  
  expect_equal(nrow(missing_data_drc), 49812)
  expect_true("missingness" %in% colnames(missing_data_drc))
  
  missingness_count_drc <- missing_data_drc %>%
    dplyr::count(missingness)
  
  expect_equal(sort(missingness_count_drc$n), c(198, 4590, 5304, 39720))
})

test_that("impute works", {
  # only test method = "ludovic" and not method = "noise". Does not test switching off log2 transformation error.
  imputed_data <- impute(missing_data, sample = sample, grouping = peptide, intensity = normalised_intensity_log2, condition = condition, comparison = comparison, missingness = missingness, method = "ludovic")

  expect_is(imputed_data, "data.frame")
  expect_equal(sum(imputed_data$imputed), 1259)
  arranged_data <- imputed_data %>%
    dplyr::filter(peptide == "peptide_1_1")
  expect_equal(round(arranged_data$imputed_intensity, digits = 1), c(15.9, 15.6, 15.6, 12.6, 12.6, 12.9))
})

protein_abundance <- calculate_protein_abundance(data = missing_data, sample = sample, protein_id = protein, precursor = peptide, intensity = normalised_intensity_log2, method = "iq", retain_columns = condition)

test_that("calculate_protein_abundance works", {
  expect_is(protein_abundance, "data.frame")
  expect_equal(round(protein_abundance$normalised_intensity_log2[1:6], digits = 2), c(16.24, 16.15, 16.23, 16.25, 16.06, 16.06))

  protein_abundance_sum <- calculate_protein_abundance(data = missing_data, sample = sample, protein_id = protein, precursor = peptide, intensity = normalised_intensity_log2, method = "sum", retain_columns = condition)
  expect_is(protein_abundance_sum, "data.frame")
  expect_equal(round(protein_abundance_sum$normalised_intensity_log2[1:6], digits = 2), c(18.21, 18.31, 18.24, 18.01, 17.82, 17.77))
})

diff <- diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity = normalised_intensity_log2, missingness = missingness, comparison = comparison, ref_condition = "condition_1", method = "t-test", retain_columns = c(protein))

test_that("diff_abundance works", {
  data_mean_sd <- missing_data %>%
    tidyr::drop_na() %>%
    dplyr::group_by(condition, peptide, protein) %>%
    dplyr::summarise(mean = mean(normalised_intensity_log2, na.rm = TRUE), sd = sd(normalised_intensity_log2, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

  diff_mean_sd <- diff_abundance(data = data_mean_sd, condition = condition, grouping = peptide, mean = mean, sd = sd, n_samples = n, ref_condition = "condition_1", method = "t-test_mean_sd", retain_columns = c(protein))
  diff_moderated <- diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity = normalised_intensity_log2, missingness = missingness, comparison = comparison, ref_condition = "condition_1", method = "moderated_t-test", retain_columns = c(protein))
  diff_proDA <- diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity = normalised_intensity_log2, missingness = missingness, comparison = comparison, ref_condition = "condition_1", method = "proDA", retain_columns = c(protein))
  
  expect_is(diff, "data.frame")
  expect_is(diff_mean_sd, "data.frame")
  expect_is(diff_moderated, "data.frame")
  expect_is(diff_proDA, "data.frame")
  expect_equal(nrow(diff), 5783)
  expect_equal(nrow(diff_mean_sd), 5746)
  expect_equal(nrow(diff_moderated), 5783)
  expect_equal(nrow(diff_proDA), 5783)
  expect_equal(ncol(diff), 9)
  expect_equal(ncol(diff_mean_sd), 14)
  expect_equal(ncol(diff_moderated), 13)
  expect_equal(ncol(diff_proDA), 12)
  expect_equal(round(min(diff$adj_pval, na.rm = TRUE), digits = 9), 0.007930396)
  expect_equal(round(min(diff_mean_sd$adj_pval, na.rm = TRUE), digits = 9), 0.007930396)
  expect_equal(round(min(diff_moderated$adj_pval, na.rm = TRUE), digits = 9), 1.3834e-05)
  expect_equal(round(min(diff_proDA$adj_pval, na.rm = TRUE), digits = 9), 0.000897949)
})

test_that("volcano_protti works", {
  sig_prots <- paste0("protein_", 1:25)
  p <- volcano_protti(data = diff, 
                      grouping = peptide, 
                      log2FC = diff, 
                      significance = adj_pval, 
                      method = "significant", 
                      protein_identifier = protein, 
                      title = "Test tile", 
                      x_axis_label = "test x-Axis", 
                      y_axis_label = "test y-Axis", 
                      log2FC_cutoff = 1, 
                      significance_cutoff = 0.05, 
                      interactive = FALSE)
  expect_is(p,"ggplot")
  expect_warning(expect_error(print(p), NA))
  
  p_interactive <- volcano_protti(data = diff, 
                      grouping = peptide, 
                      log2FC = diff, 
                      significance = adj_pval, 
                      method = "significant", 
                      protein_identifier = protein, 
                      title = "Test tile", 
                      x_axis_label = "test x-Axis", 
                      y_axis_label = "test y-Axis", 
                      log2FC_cutoff = 1, 
                      significance_cutoff = 0.05, 
                      interactive = TRUE)
  expect_is(p_interactive,"plotly")
  expect_error(print(p_interactive), NA)
  
  p_target <- volcano_protti(data = diff, 
                      grouping = peptide, 
                      log2FC = diff, 
                      significance = adj_pval, 
                      method = "target", 
                      target = "protein_3",
                      protein_identifier = protein, 
                      title = "Test tile", 
                      x_axis_label = "test x-Axis", 
                      y_axis_label = "test y-Axis", 
                      log2FC_cutoff = 1, 
                      significance_cutoff = 0.05, 
                      interactive = FALSE)
  expect_is(p_target,"ggplot")
  expect_warning(expect_error(print(p_target), NA))
})

test_that("protein_abundance_normalisation works", {
  # does not test the paired argument. No change of the completeness argument is tested
  protein_abundance <- protein_abundance %>% 
    dplyr::rename(protein_intensity = normalised_intensity_log2) %>% 
    assign_missingness(sample = sample,
                       condition = condition,
                       grouping = protein,
                       intensity = protein_intensity,
                       ref_condition = "condition_1")
  
  abundance_normalised_ref <- protein_abundance_normalisation(lip_peptides = missing_data, 
                                                          tc_proteins = protein_abundance, 
                                                          sample = sample,
                                                          grouping = peptide,
                                                          condition = condition, 
                                                          protein_id = protein,
                                                          peptide_intensity = normalised_intensity_log2,
                                                          protein_intensity = protein_intensity,
                                                          method = "ref_vs_rest",
                                                          comparison = comparison)
  expect_is(abundance_normalised_ref, "data.frame")
  expect_equal(nrow(abundance_normalised_ref), 5648)
  expect_equal(ncol(abundance_normalised_ref), 14)
  expect_equal(round(min(abundance_normalised_ref$adj_pval, na.rm = TRUE), digits = 9), 0.01394052)
  
  abundance_normalised_all <- protein_abundance_normalisation(lip_peptides = missing_data, 
                                                              tc_proteins = protein_abundance, 
                                                              sample = sample,
                                                              grouping = peptide,
                                                              condition = condition, 
                                                              protein_id = protein,
                                                              peptide_intensity = normalised_intensity_log2,
                                                              protein_intensity = protein_intensity,
                                                              method = "all_vs_all",
                                                              comparison = comparison)
  expect_is(abundance_normalised_all, "data.frame")
  expect_equal(nrow(abundance_normalised_all), 5648)
  expect_equal(ncol(abundance_normalised_all), 15)
  expect_equal(round(min(abundance_normalised_all$adj_pval, na.rm = TRUE), digits = 9), 0.002556753)
  expect_equal(round(min(abundance_normalised_all$pval_tukey, na.rm = TRUE), digits = 9), 0)
})

drc_fit <- fit_drc_4p(data = missing_data_drc, sample = sample, grouping = peptide, response = normalised_intensity_log2, dose = concentration)

test_that("fit_drc_4p works", {
  # did not test the argument include_models = TRUE
  expect_is(drc_fit, "data.frame")
  expect_equal(nrow(drc_fit), 1067)
  expect_equal(ncol(drc_fit), 9)
  expect_equal(round(max(drc_fit$correlation, na.rm = TRUE), digits = 3), 0.995)
  expect_equal(round(min(drc_fit$p_value, na.rm = TRUE), digits = 40), 2.318708e-21)
})

test_that("plot_drc_4p works", {
  p <- plot_drc_4p(data = drc_fit, 
                   grouping = peptide,
                   response = normalised_intensity_log2,
                   dose = concentration,
                   targets = c("peptide_4_1"),
                   unit = "uM",
                   y_axis_name = "test y-Axis")
  
  expect_is(p,"ggplot")
  expect_warning(expect_error(print(p), NA))
  
  p_facet <- plot_drc_4p(data = drc_fit, 
                   grouping = peptide,
                   response = normalised_intensity_log2,
                   dose = concentration,
                   targets = c("peptide_2_1", "peptide_2_4"),
                   unit = "uM")
  
  expect_is(p_facet,"ggplot")
  expect_warning(expect_error(print(p_facet), NA))
})