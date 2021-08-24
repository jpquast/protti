context("test-workflow")

set.seed(123)
data <- create_synthetic_data(
  n_proteins = 50,
  frac_change = 0.05,
  n_replicates = 3,
  n_conditions = 2,
  method = "effect_random",
  additional_metadata = FALSE
)

data_drc <- create_synthetic_data(
  n_proteins = 20,
  frac_change = 0.05,
  n_replicates = 3,
  n_conditions = 8,
  method = "dose_response",
  concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
  additional_metadata = FALSE
)

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("deprecated median_normalization works", {
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(normalised_data <- data %>%
        median_normalisation(sample = sample, intensity_log2 = peptide_intensity_missing))
    })
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(normalised_data_drc <- data_drc %>%
        median_normalisation(sample = sample, intensity_log2 = peptide_intensity_missing))
    })

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
}

normalised_data <- data %>%
  normalise(sample = sample, intensity_log2 = peptide_intensity_missing, method = "median")

normalised_data_drc <- data_drc %>%
  normalise(sample = sample, intensity_log2 = peptide_intensity_missing, method = "median")

test_that("normalise works", {
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

  if (Sys.getenv("TEST_PROTTI") == "true") {
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
  }
})

missing_data <- normalised_data %>%
  assign_missingness(sample = sample, condition = condition, grouping = peptide, intensity = normalised_intensity_log2, ref_condition = "condition_1", retain_columns = c(protein))

test_that("assign_missingness works", {
  # not testing noise argument. Also no change of default values for completeness_MAR and completeness_MNAR
  expect_equal(nrow(data), nrow(missing_data))
  expect_true("missingness" %in% colnames(missing_data))

  missingness_count <- missing_data %>%
    dplyr::count(missingness)

  expect_equal(sort(missingness_count$n), c(12, 516, 618, 3078))
})

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("impute works", {
    # only test method = "ludovic" and not method = "noise". Does not test switching off log2 transformation error.
    imputed_data <- impute(missing_data, sample = sample, grouping = peptide, intensity = normalised_intensity_log2, condition = condition, comparison = comparison, missingness = missingness, method = "ludovic", retain_columns = protein)

    expect_is(imputed_data, "data.frame")
    expect_equal(sum(imputed_data$imputed), 115)
    arranged_data <- imputed_data %>%
      dplyr::filter(peptide == "peptide_1_4")
    expect_false(is.na(arranged_data$imputed_intensity[4]))
  })
}

protein_abundance <- calculate_protein_abundance(data = missing_data, sample = sample, protein_id = protein, precursor = peptide, peptide = peptide, intensity_log2 = normalised_intensity_log2, method = "iq", retain_columns = condition)
protein_abundance_all <- calculate_protein_abundance(data = missing_data, sample = sample, protein_id = protein, precursor = peptide, peptide = peptide, intensity_log2 = normalised_intensity_log2, method = "sum", for_plot = TRUE)

test_that("calculate_protein_abundance works", {
  arranged_data <- protein_abundance %>%
    dplyr::filter(protein == "protein_1")
  expect_is(protein_abundance, "data.frame")
  expect_equal(round(arranged_data$normalised_intensity_log2, digits = 2), c(16.78, 16.94, 16.85, 16.81, 16.83, 16.82))
  expect_equal(nrow(protein_abundance), 288)
  expect_equal(ncol(protein_abundance), 4)

  arranged_data <- protein_abundance_all %>%
    dplyr::filter(protein == "protein_1" & peptide == "protein_intensity")
  expect_equal(round(arranged_data$normalised_intensity_log2, digits = 2), c(20.87, 20.97, 20.96, 20.81, 20.81, 20.86))
  expect_is(protein_abundance_all, "data.frame")
  expect_equal(nrow(protein_abundance_all), 4032)
  expect_equal(ncol(protein_abundance_all), 4)
})

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("deprecated plot_peptide_profiles works", {
    # not testing split_all = TRUE. Also not testing protein_abundance_plot = FALSE.
    expect_warning(p <- plot_peptide_profiles(
      data = protein_abundance_all,
      sample = sample,
      peptide = peptide,
      intensity_log2 = normalised_intensity_log2,
      grouping = protein,
      targets = c("protein_1", "protein_2"),
      protein_abundance_plot = TRUE
    ))
    expect_is(p, "list")
    expect_is(p[[1]], "ggplot")
    expect_error(print(p[[1]]), NA)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p_peptides <- plot_peptide_profiles(
        data = missing_data,
        sample = sample,
        peptide = peptide,
        intensity_log2 = normalised_intensity_log2,
        grouping = protein,
        targets = c("protein_12"),
        protein_abundance_plot = FALSE
      ))
    })

    expect_is(p_peptides, "list")
    expect_is(p_peptides[[1]], "ggplot")
    expect_error(print(p_peptides[[1]]), NA)
  })
}

test_that("peptide_profile_plot works", {
  # not testing split_all = TRUE. Also not testing protein_abundance_plot = FALSE.
  p <- peptide_profile_plot(
    data = protein_abundance_all,
    sample = sample,
    peptide = peptide,
    intensity_log2 = normalised_intensity_log2,
    grouping = protein,
    targets = c("protein_1", "protein_2"),
    protein_abundance_plot = TRUE
  )
  expect_is(p, "list")
  expect_is(p[[1]], "ggplot")
  expect_error(print(p[[1]]), NA)

  if (Sys.getenv("TEST_PROTTI") == "true") {
    p_peptides <- peptide_profile_plot(
      data = missing_data,
      sample = sample,
      peptide = peptide,
      intensity_log2 = normalised_intensity_log2,
      grouping = protein,
      targets = c("protein_12"),
      protein_abundance_plot = FALSE
    )

    expect_is(p_peptides, "list")
    expect_is(p_peptides[[1]], "ggplot")
    expect_error(print(p_peptides[[1]]), NA)
  }
})

diff <- calculate_diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity_log2 = normalised_intensity_log2, missingness = missingness, comparison = comparison, method = "t-test", retain_columns = c(protein))

test_that("calculate_diff_abundance works", {
  expect_is(diff, "data.frame")
  expect_equal(nrow(diff), 601)
  expect_equal(ncol(diff), 9)
  expect_equal(round(min(diff$adj_pval, na.rm = TRUE), digits = 9), 0.00758761)

  if (Sys.getenv("TEST_PROTTI") == "true") {
    data_mean_sd <- missing_data %>%
      tidyr::drop_na() %>%
      dplyr::group_by(condition, peptide, protein) %>%
      dplyr::summarise(mean = mean(normalised_intensity_log2, na.rm = TRUE), sd = sd(normalised_intensity_log2, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

    diff_mean_sd <- calculate_diff_abundance(data = data_mean_sd, condition = condition, grouping = peptide, mean = mean, sd = sd, n_samples = n, ref_condition = "condition_1", method = "t-test_mean_sd", retain_columns = c(protein))
    diff_moderated <- calculate_diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity_log2 = normalised_intensity_log2, missingness = missingness, comparison = comparison, method = "moderated_t-test", retain_columns = c(protein))
    diff_proDA <- calculate_diff_abundance(data = missing_data, sample = sample, condition = condition, grouping = peptide, intensity_log2 = normalised_intensity_log2, missingness = missingness, comparison = comparison, method = "proDA", retain_columns = c(protein))

    expect_is(diff_mean_sd, "data.frame")
    expect_is(diff_moderated, "data.frame")
    expect_is(diff_proDA, "data.frame")
    expect_equal(nrow(diff_mean_sd), 599)
    expect_equal(nrow(diff_moderated), 601)
    expect_equal(nrow(diff_proDA), 601)
    expect_equal(ncol(diff_mean_sd), 14)
    expect_equal(ncol(diff_moderated), 13)
    expect_equal(ncol(diff_proDA), 12)
    expect_equal(round(min(diff_mean_sd$adj_pval, na.rm = TRUE), digits = 9), 0.00758761)
    expect_equal(round(min(diff_moderated$adj_pval, na.rm = TRUE), digits = 9), 5.1129e-05)
    expect_equal(round(min(diff_proDA$adj_pval, na.rm = TRUE), digits = 5), 0.00125)
  }
})

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("deprecated diff_abundance works", {
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(diff_deprecated <- diff_abundance(
        data = missing_data,
        sample = sample,
        condition = condition,
        grouping = peptide,
        intensity_log2 = normalised_intensity_log2,
        missingness = missingness,
        comparison = comparison,
        method = "t-test",
        retain_columns = c(protein)
      ))
    })
    expect_is(diff_deprecated, "data.frame")
    expect_equal(nrow(diff_deprecated), 601)
    expect_equal(ncol(diff_deprecated), 9)
    expect_equal(round(min(diff_deprecated$adj_pval, na.rm = TRUE), digits = 9), 0.00758761)

    data_mean_sd <- missing_data %>%
      tidyr::drop_na() %>%
      dplyr::group_by(condition, peptide, protein) %>%
      dplyr::summarise(mean = mean(normalised_intensity_log2, na.rm = TRUE), sd = sd(normalised_intensity_log2, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(diff_mean_sd_deprecated <- diff_abundance(
        data = data_mean_sd,
        condition = condition,
        grouping = peptide,
        mean = mean,
        sd = sd,
        n_samples = n,
        ref_condition = "condition_1",
        method = "t-test_mean_sd",
        retain_columns = c(protein)
      ))
    })
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(diff_moderated_deprecated <- diff_abundance(
        data = missing_data,
        sample = sample,
        condition = condition,
        grouping = peptide,
        intensity_log2 = normalised_intensity_log2,
        missingness = missingness,
        comparison = comparison,
        method = "moderated_t-test",
        retain_columns = c(protein)
      ))
    })
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(diff_proDA_deprecated <- diff_abundance(
        data = missing_data,
        sample = sample,
        condition = condition,
        grouping = peptide,
        intensity_log2 = normalised_intensity_log2,
        missingness = missingness,
        comparison = comparison,
        method = "proDA",
        retain_columns = c(protein)
      ))
    })

    expect_is(diff_mean_sd_deprecated, "data.frame")
    expect_is(diff_moderated_deprecated, "data.frame")
    expect_is(diff_proDA_deprecated, "data.frame")
    expect_equal(nrow(diff_mean_sd_deprecated), 599)
    expect_equal(nrow(diff_moderated_deprecated), 601)
    expect_equal(nrow(diff_proDA_deprecated), 601)
    expect_equal(ncol(diff_mean_sd_deprecated), 14)
    expect_equal(ncol(diff_moderated_deprecated), 13)
    expect_equal(ncol(diff_proDA_deprecated), 12)
    expect_equal(round(min(diff_mean_sd_deprecated$adj_pval, na.rm = TRUE), digits = 9), 0.00758761)
    expect_equal(round(min(diff_moderated_deprecated$adj_pval, na.rm = TRUE), digits = 9), 5.1129e-05)
    expect_equal(round(min(diff_proDA_deprecated$adj_pval, na.rm = TRUE), digits = 5), 0.00125)
  })
}

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("depreciated plot_pval_distribution works", {
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- plot_pval_distribution(
        diff,
        peptide,
        pval
      ))
    })
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })
}

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("pval_distribution_plot works", {
    p <- pval_distribution_plot(
      diff,
      peptide,
      pval
    )
    expect_is(p, "ggplot")
    expect_error(print(p), NA)
  })
}

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("deprecated volcano_protti works", {
    sig_prots <- paste0("protein_", 1:25)
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- volcano_protti(
        data = diff,
        grouping = peptide,
        log2FC = diff,
        significance = adj_pval,
        method = "significant",
        target_column = protein,
        title = "Test tile",
        x_axis_label = "test x-Axis",
        y_axis_label = "test y-Axis",
        log2FC_cutoff = 1,
        significance_cutoff = 0.05,
        interactive = FALSE
      ))
    })
    expect_is(p, "ggplot")
    expect_error(print(p), NA)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p_interactive <- volcano_protti(
        data = diff,
        grouping = peptide,
        log2FC = diff,
        significance = adj_pval,
        method = "significant",
        target_column = protein,
        title = "Test tile",
        x_axis_label = "test x-Axis",
        y_axis_label = "test y-Axis",
        log2FC_cutoff = 1,
        significance_cutoff = 0.05,
        interactive = TRUE
      ))
    })
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p_target <- volcano_protti(
        data = diff,
        grouping = peptide,
        log2FC = diff,
        significance = adj_pval,
        method = "target",
        target = "protein_3",
        target_column = protein,
        title = "Test tile",
        x_axis_label = "test x-Axis",
        y_axis_label = "test y-Axis",
        log2FC_cutoff = 1,
        significance_cutoff = 0.05,
        interactive = FALSE
      ))
    })
    expect_is(p_target, "ggplot")
    expect_error(print(p_target), NA)
  })
}

test_that("volcano_plot works", {
  sig_prots <- paste0("protein_", 1:25)
  p <- volcano_plot(
    data = diff,
    grouping = peptide,
    log2FC = diff,
    significance = adj_pval,
    method = "significant",
    target_column = protein,
    title = "Test tile",
    x_axis_label = "test x-Axis",
    y_axis_label = "test y-Axis",
    log2FC_cutoff = 1,
    significance_cutoff = 0.05,
    interactive = FALSE
  )
  expect_is(p, "ggplot")
  expect_error(print(p), NA)

  if (Sys.getenv("TEST_PROTTI") == "true") {
    p_interactive <- volcano_plot(
      data = diff,
      grouping = peptide,
      log2FC = diff,
      significance = adj_pval,
      method = "significant",
      target_column = protein,
      title = "Test tile",
      x_axis_label = "test x-Axis",
      y_axis_label = "test y-Axis",
      log2FC_cutoff = 1,
      significance_cutoff = 0.05,
      interactive = TRUE
    )
    expect_is(p_interactive, "plotly")
    expect_error(print(p_interactive), NA)
  }

  p_target <- volcano_plot(
    data = diff,
    grouping = peptide,
    log2FC = diff,
    significance = adj_pval,
    method = "target",
    target = "protein_3",
    target_column = protein,
    title = "Test tile",
    x_axis_label = "test x-Axis",
    y_axis_label = "test y-Axis",
    log2FC_cutoff = 1,
    significance_cutoff = 0.05,
    interactive = FALSE
  )
  expect_is(p_target, "ggplot")
  expect_error(print(p_target), NA)
})

drc_fit <- fit_drc_4p(data = normalised_data_drc, sample = sample, grouping = peptide, response = normalised_intensity_log2, dose = concentration, log_logarithmic = TRUE, retain_columns = c(protein))

test_that("fit_drc_4p works", {
  # did not test the argument include_models = TRUE
  expect_is(drc_fit, "data.frame")
  expect_equal(nrow(drc_fit), 282)
  expect_equal(ncol(drc_fit), 18)
  expect_equal(round(max(drc_fit$correlation, na.rm = TRUE), digits = 3), 0.876)
  expect_equal(round(min(drc_fit$anova_pval, na.rm = TRUE), digits = 6), 0.007297)
})

if (Sys.getenv("TEST_PROTTI") == "true") {
  test_that("deprecated plot_drc_4p works", {
    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p <- plot_drc_4p(
        data = drc_fit,
        grouping = peptide,
        response = normalised_intensity_log2,
        dose = concentration,
        targets = c("peptide_1_2"),
        unit = "uM",
        y_axis_name = "test y-Axis"
      ))
    })

    expect_is(p, "ggplot")
    expect_warning(print(p), NA)

    rlang::with_options(lifecycle_verbosity = "warning", {
      expect_warning(p_facet <- plot_drc_4p(
        data = drc_fit,
        grouping = peptide,
        response = normalised_intensity_log2,
        dose = concentration,
        targets = c("peptide_1_2", "peptide_1_1"),
        unit = "uM"
      ))
    })

    expect_is(p_facet, "list")
    expect_warning(expect_error(print(p_facet), NA))
  })
}

test_that("drc_4p_plot works", {
  p <- drc_4p_plot(
    data = drc_fit,
    grouping = peptide,
    response = normalised_intensity_log2,
    dose = concentration,
    targets = c("peptide_1_2"),
    unit = "uM",
    y_axis_name = "test y-Axis"
  )

  expect_is(p, "ggplot")
  expect_warning(print(p), NA)

  p_facet <- drc_4p_plot(
    data = drc_fit,
    grouping = peptide,
    response = normalised_intensity_log2,
    dose = concentration,
    targets = c("peptide_1_2", "peptide_1_1"),
    unit = "uM"
  )

  expect_is(p_facet, "list")
  expect_warning(expect_error(print(p_facet), NA))
})

test_that("filter_cv works", {
  normalised_data_filtered <- normalised_data %>%
    filter_cv(peptide, condition, peptide_intensity_missing, cv_limit = 0.25, min_conditions = 2)

  normalised_data_filtered_cv_count <- normalised_data_filtered %>%
    dplyr::group_by(peptide, condition) %>%
    dplyr::summarise(cv_count = sum((sd(2^peptide_intensity_missing, na.rm = TRUE) / mean(2^peptide_intensity_missing, na.rm = TRUE)) > 0.25), .groups = "drop") %>%
    dplyr::filter(cv_count > 0)

  expect_is(normalised_data_filtered, "data.frame")
  expect_gt(nrow(normalised_data), nrow(normalised_data_filtered))
  expect_equal(nrow(normalised_data_filtered_cv_count), 0)

  normalised_data_drc_filtered <- normalised_data_drc %>%
    filter_cv(peptide, condition, peptide_intensity_missing, cv_limit = 0.25, min_conditions = 6)

  normalised_data_filtered_drc_cv_count <- normalised_data_drc_filtered %>%
    dplyr::group_by(peptide, condition) %>%
    dplyr::summarise(cv_count = sum((sd(2^peptide_intensity_missing, na.rm = TRUE) / mean(2^peptide_intensity_missing, na.rm = TRUE)) > 0.25), .groups = "drop") %>%
    dplyr::filter(cv_count > 0)

  expect_is(normalised_data_drc_filtered, "data.frame")
  expect_gt(nrow(normalised_data_drc), nrow(normalised_data_drc_filtered))
  expect_gt(nrow(normalised_data_filtered_drc_cv_count), 0)
})

test_that("calculate_aa_scores works", {
  peptide_data <- tibble::tibble(
    adj_pval = c(0.001, 0.0001, 0.4, 0.2, 0.001, 0.0001, 0.001, 0.0001, 0.4, 0.2, 0.001, 0.0001),
    diff = c(4.3, -5.8, 0.23, -0.5, 6.5, -7.3, 4.3, -5.8, 0.23, -0.5, 6.5, -7.3),
    pep_stripped_sequence = c(
      "MDVFMKGLS",
      "MDVFMK",
      "SKAKEGVV",
      "EGVVAAAE",
      "EGVVAAA",
      "EKTKQG",
      "MSLSNKL",
      "LSNKLTLD",
      "TLDKLDV",
      "DKLDV",
      "KGKRVV",
      "GKRVV"
    ),
    start = c(1, 1, 8, 12, 12, 18, 1, 3, 8, 10, 14, 15),
    end = c(10, 6, 16, 20, 19, 25, 7, 10, 15, 15, 20, 20),
    uniprot_id = c(
      "P37840",
      "P37840",
      "P37840",
      "P37840",
      "P37840",
      "P37840",
      "P00558",
      "P00558",
      "P00558",
      "P00558",
      "P00558",
      "P00558"
    )
  )

  aa_fingerprint <- calculate_aa_scores(peptide_data,
    protein = uniprot_id,
    grouping = pep_stripped_sequence,
    diff = diff,
    adj_pval = adj_pval,
    start = start,
    end = end
  )

  expect_is(aa_fingerprint, "data.frame")
  expect_equal(nrow(aa_fingerprint), 45)
  expect_equal(ncol(aa_fingerprint), 3)
})
