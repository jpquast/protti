#' Calculate differential abundance between conditions
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `calculate_diff_abundance()`.
#'
#' @keywords internal
#' @export
diff_abundance <-
  function(...) {
    # This function has been renamed and is therefore deprecated.
    lifecycle::deprecate_warn("0.2.0",
      "diff_abundance()",
      "calculate_diff_abundance()",
      details = "This function has been renamed."
    )

    calculate_diff_abundance(...)
  }
#' Calculate differential abundance between conditions
#'
#' Performs differential abundance calculations and statistical hypothesis tests on data frames
#' with protein, peptide or precursor data. Different methods for statistical testing are available.
#'
#' @param data a data frame containing at least the input variables that are required for the
#' selected method. Ideally the output of \code{assign_missingness} or \code{impute} is used.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' Is not required if \code{method = "t-test_mean_sd"}.
#' @param condition a character or numeric column in the \code{data} data frame that contains the
#' conditions.
#' @param grouping a character column in the \code{data} data frame that contains precursor or
#' peptide identifiers.
#' @param intensity_log2 a numeric column in the \code{data} data frame that contains intensity
#' values. The intensity values need to be log2 transformed. Is not required if
#' \code{method = "t-test_mean_sd"}.
#' @param missingness a character column in the \code{data} data frame that contains missingness
#' information. Can be obtained by calling \code{assign_missingness()}. Is not required if
#' \code{method = "t-test_mean_sd"}. The type of missingness assigned to a comparison does not have
#' any influence on the statistical test. However, if \code{filter_NA_missingness = TRUE} then
#' comparisons with missingness \code{NA} are filtered out prior to p-value adjustment.
#' @param comparison a character column in the \code{data} data frame that contains information of
#' treatment/reference condition pairs. Can be obtained by calling \code{assign_missingness}.
#' Comparisons need to be in the form condition1_vs_condition2, meaning two compared conditions are
#' separated by \code{"_vs_"}. This column determines for which condition pairs differential
#' abundances are calculated. Is not required if \code{method = "t-test_mean_sd"}, in that case
#' please provide a reference condition with the ref_condition argument.
#' @param mean a numeric column in the \code{data} data frame that contains mean values for two
#' conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param sd a numeric column in the \code{data} data frame that contains standard deviations for
#' two conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param n_samples a numeric column in the \code{data} data frame that contains the number of
#' samples per condition for two conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param ref_condition optional, character value providing the condition that is used as a
#' reference for differential abundance calculation. Only required for \code{method = "t-test_mean_sd"}.
#' Instead of providing one reference condition, "all" can be supplied, which will create all
#' pairwise condition pairs. By default \code{ref_condition = "all"}.
#' @param filter_NA_missingness a logical value, default is \code{TRUE}. For all methods except
#' \code{"t-test_mean_sd"} missingness information has to be provided. This information can be
#' for example obtained by calling \code{assign_missingness()}. If a reference/treatment pair has
#' too few samples to be considered robust based on user defined cutoffs, it is annotated with \code{NA}
#' as missingness by the \code{assign_missingness()} function. If this argument is \code{TRUE},
#' these \code{NA} reference/treatment pairs are filtered out after the testing and prior to p-value
#' adjustment.
#' @param method a character value, specifies the method used for statistical hypothesis testing.
#' Methods include Welch test (\code{"t-test"}), a Welch test on means, standard deviations and
#' number of replicates ("\code{"t-test_mean_sd"}) and a moderated t-test based on the \code{limma}
#' package (\code{"moderated_t-test"}). More information on the moderated t-test can be found in
#' the \code{limma} documentation. Furthermore, the \code{proDA} package specific method (\code{"proDA"})
#' can be used to infer means across samples based on a probabilistic dropout model. This
#' eliminates the need for data imputation since missing values are inferred from the model. More
#' information can be found in the \code{proDA} documentation. We do not recommend using the
#' \code{moderated_t-test} or \code{proDA} method if the data was filtered for low CVs or imputation
#' was performed. Default is \code{method = "moderated_t-test"}.
#' @param p_adj_method a character value, specifies the p-value correction method. Possible
#' methods are c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default
#' method is \code{"BH"}.
#' @param retain_columns a vector indicating if certain columns should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return A data frame that contains differential abundances (\code{diff}), p-values (\code{pval})
#' and adjusted p-values (\code{adj_pval}) for each protein, peptide or precursor (depending on
#' the \code{grouping} variable) and the associated treatment/reference pair. Depending on the
#' method the data frame contains additional columns:
#' \itemize{
#' \item{"t-test": }{The \code{std_error} column contains the standard error of the differential
#' abundances. \code{n_obs} contains the number of observations for the specific protein, peptide
#' or precursor (depending on the \code{grouping} variable) and the associated treatment/reference pair.}
#' \item{"t-test_mean_sd": }{Columns labeled as control refer to the second condition of the
#' comparison pairs. Treated refers to the first condition. \code{mean_control} and \code{mean_treated}
#' columns contain the means for the reference and treatment condition, respectively. \code{sd_control}
#' and \code{sd_treated} columns contain the standard deviations for the reference and treatment
#' condition, respectively. \code{n_control} and \code{n_treated} columns contain the numbers of
#' samples for the reference and treatment condition, respectively. The \code{std_error} column
#' contains the standard error of the differential abundances. \code{t_statistic} contains the
#' t_statistic for the t-test.}
#' \item{"moderated_t-test": }{\code{CI_2.5} and \code{CI_97.5} contain the 2.5% and 97.5%
#' confidence interval borders for differential abundances. \code{avg_abundance} contains average
#' abundances for treatment/reference pairs (mean of the two group means). \code{t_statistic}
#' contains the t_statistic for the t-test. \code{B} The B-statistic is the log-odds that the
#' protein, peptide or precursor (depending on \code{grouping}) has a differential abundance
#' between the two groups. Suppose B=1.5. The odds of differential abundance is exp(1.5)=4.48, i.e,
#' about four and a half to one. The probability that there is a differential abundance is
#' 4.48/(1+4.48)=0.82, i.e., the probability is about 82% that this group is differentially
#' abundant. A B-statistic of zero corresponds to a 50-50 chance that the group is differentially
#' abundant.\code{n_obs} contains the number of observations for the specific protein, peptide or
#' precursor (depending on the \code{grouping} variable) and the associated treatment/reference pair.}
#' \item{"proDA": }{The \code{std_error} column contains the standard error of the differential
#' abundances. \code{avg_abundance} contains average abundances for treatment/reference pairs
#' (mean of the two group means). \code{t_statistic} contains the t_statistic for the t-test.
#' \code{n_obs} contains the number of observations for the specific protein, peptide or precursor
#' (depending on the \code{grouping} variable) and the associated treatment/reference pair.}
#' }
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data enquo ensym as_name as_label expr := !!
#' @importFrom purrr map map2 map_df map_dbl map2_dbl reduce set_names
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames rownames_to_column tibble as_tibble
#' @importFrom stringr str_replace_all str_extract
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create synthetic data
#' data <- create_synthetic_data(
#'   n_proteins = 10,
#'   frac_change = 0.5,
#'   n_replicates = 4,
#'   n_conditions = 2,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' # Assign missingness information
#' data_missing <- assign_missingness(
#'   data,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   ref_condition = "all",
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' # Calculate differential abundances
#' # Using "moderated_t-test" and "proDA" improves
#' # true positive recovery progressively
#' diff <- calculate_diff_abundance(
#'   data = data_missing,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity_missing,
#'   missingness = missingness,
#'   comparison = comparison,
#'   method = "t-test",
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' head(diff, n = 10)
calculate_diff_abundance <-
  function(data,
           sample,
           condition,
           grouping,
           intensity_log2,
           missingness = missingness,
           comparison = comparison,
           mean = NULL,
           sd = NULL,
           n_samples = NULL,
           ref_condition = "all",
           filter_NA_missingness = TRUE,
           method = c("moderated_t-test", "t-test", "t-test_mean_sd", "proDA"),
           p_adj_method = "BH",
           retain_columns = NULL) {
    . <- NULL
    if (!(ref_condition %in% unique(pull(data, {{ condition }}))) & ref_condition != "all") {
      stop(strwrap("The name provided to ref_condition cannot be found in your conditions!
Please provide a valid reference condition.", prefix = "\n", initial = ""))
    }

    method <- match.arg(method)

    if (method != "t-test_mean_sd") {
      if (max(pull(data, {{ intensity_log2 }}), na.rm = TRUE) > 1000) {
        stop("Please log2 transform your intensities.")
      }
    }
    if (method == "t-test_mean_sd") {
      if (max(pull(data, {{ mean }}), na.rm = TRUE) > 1000) {
        stop("Please log2 transform your data.")
      }
    }

    if (method == "t-test") {
      message("[1/2] Create input for t-tests ... ", appendLF = FALSE)

      t_test_missingness_obs <- data %>%
        tidyr::drop_na({{ missingness }}, {{ intensity_log2 }}) %>%
        dplyr::group_by({{ comparison }}, {{ grouping }}) %>%
        dplyr::mutate(n_obs = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::distinct({{ grouping }}, {{ comparison }}, {{ missingness }}, .data$n_obs)

      t_test_input <- data %>%
        tidyr::drop_na({{ intensity_log2 }}) %>%
        dplyr::group_by({{ comparison }}, {{ grouping }}, {{ condition }}) %>%
        dplyr::summarize(intensity = list({{ intensity_log2 }}), .groups = "drop") %>%
        dplyr::mutate(type = ifelse({{ condition }} == stringr::str_extract({{ comparison }}, pattern = "(?<=_vs_).+"),
          "control",
          "treated"
        )) %>%
        dplyr::select(-{{ condition }}) %>%
        tidyr::pivot_wider(names_from = .data$type, values_from = .data$intensity, values_fill = list(NA))

      message("DONE", appendLF = TRUE)
      message("[2/2] Calculate t-tests ... ", appendLF = FALSE)

      t_test_result <- t_test_input %>%
        dplyr::mutate(pval = purrr::map2(
          .x = .data$treated,
          .y = .data$control,
          .f = function(.x, .y) {
            tryCatch(
              {
                suppressWarnings(stats::t.test(.x, .y))
              },
              error = function(error) {
                NA
              }
            )
          }
        )) %>%
        dplyr::mutate(std_error = purrr::map_dbl(
          .x = .data$pval,
          .f = ~ tryCatch(
            {
              .x$stderr
            },
            error = function(error) {
              NA
            }
          )
        )) %>%
        dplyr::mutate(pval = purrr::map_dbl(
          .x = .data$pval,
          .f = ~ tryCatch(
            {
              .x$p.value
            },
            error = function(error) {
              NA
            }
          )
        )) %>%
        dplyr::mutate(diff = map2_dbl(
          .x = .data$treated,
          .y = .data$control,
          .f = function(.x, .y) {
            suppressWarnings(mean(.x, na.rm = TRUE)) - suppressWarnings(mean(.y, na.rm = TRUE))
          }
        )) %>%
        dplyr::mutate(diff = ifelse(diff == "NaN", NA, diff)) %>%
        dplyr::group_by({{ comparison }}) %>%
        dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(.data$control, .data$treated)) %>%
        dplyr::left_join(t_test_missingness_obs, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison")) %>%
        dplyr::arrange(.data$adj_pval, .data$pval)

      message("DONE", appendLF = TRUE)

      if (!missing(retain_columns)) {
        t_test_result <- data %>%
          dplyr::ungroup() %>%
          dplyr::select(
            !!enquo(retain_columns),
            {{ intensity_log2 }},
            colnames(t_test_result)[!colnames(t_test_result) %in%
              c(
                "pval",
                "std_error",
                "diff",
                "adj_pval",
                "n_obs"
              )]
          ) %>%
          tidyr::drop_na({{ intensity_log2 }}) %>%
          dplyr::select(-{{ intensity_log2 }}) %>%
          dplyr::distinct() %>%
          dplyr::right_join(t_test_result, by = colnames(t_test_result)[!colnames(t_test_result) %in%
            c(
              "pval",
              "std_error",
              "diff",
              "adj_pval",
              "n_obs"
            )]) %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
      }

      if (filter_NA_missingness == TRUE) {
        t_test_result <- t_test_result %>%
          tidyr::drop_na({{ missingness }}) %>%
          dplyr::group_by({{ comparison }}) %>%
          dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method)) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
        return(t_test_result)
      }
      if (filter_NA_missingness == FALSE) {
        return(t_test_result)
      }
    }

    if (method == "t-test_mean_sd") {
      if (ref_condition == "all") {
        # creating all pairwise comparisons
        all_conditions <- unique(dplyr::pull(data, {{ condition }}))

        all_combinations <- tibble::as_tibble(t(combn(all_conditions, m = 2))) %>%
          dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2))

        message(
          strwrap("All pairwise comparisons are created from all conditions and their
missingness type is assigned.\n The created comparisons are: \n", prefix = "\n", initial = ""),
          paste(all_combinations$combinations, collapse = "\n")
        )
      }

      if (ref_condition != "all") {
        conditions_no_ref <- unique(dplyr::pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition]

        all_combinations <- tibble::tibble(V1 = conditions_no_ref, V2 = ref_condition) %>%
          dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2))
      }

      all_combinations <- all_combinations %>%
        tidyr::pivot_longer(cols = c(.data$V1, .data$V2), names_to = "name", values_to = rlang::as_name(rlang::enquo(condition))) %>%
        dplyr::select(-.data$name) %>%
        dplyr::group_by({{ condition }}) %>%
        dplyr::mutate(comparison = list(.data$combinations)) %>%
        dplyr::distinct(.data$comparison, {{ condition }})

      t_test_mean_sd_result <- data %>%
        dplyr::ungroup() %>%
        dplyr::distinct({{ condition }}, {{ grouping }}, {{ mean }}, {{ sd }}, {{ n_samples }}) %>%
        tidyr::drop_na() %>%
        dplyr::left_join(all_combinations, by = rlang::as_name(rlang::enquo(condition))) %>%
        tidyr::unnest(.data$comparison) %>%
        dplyr::rename(mean = {{ mean }}, sd = {{ sd }}, n = {{ n_samples }}) %>%
        dplyr::mutate({{ condition }} := ifelse({{ condition }} == stringr::str_extract(.data$comparison, pattern = "(?<=_vs_).+"),
          "control",
          "treated"
        )) %>%
        tidyr::pivot_wider(names_from = {{ condition }}, values_from = c(.data$mean, .data$sd, .data$n)) %>%
        dplyr::mutate(ttest_protti(
          mean1 = .data$mean_control,
          mean2 = .data$mean_treated,
          sd1 = .data$sd_control,
          sd2 = .data$sd_treated,
          n1 = .data$n_control,
          n2 = .data$n_treated
        )) %>%
        tidyr::drop_na(.data$pval) %>%
        dplyr::group_by(.data$comparison) %>%
        dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method)) %>%
        dplyr::arrange(.data$adj_pval, .data$pval)

      if (!missing(retain_columns)) {
        t_test_mean_sd_result <- data %>%
          dplyr::ungroup() %>%
          dplyr::select(!!enquo(retain_columns), colnames(t_test_mean_sd_result)[!colnames(t_test_mean_sd_result) %in%
            c(
              "mean_control",
              "mean_treated",
              "sd_control",
              "sd_treated",
              "n_control",
              "n_treated",
              "pval",
              "std_error",
              "diff",
              "adj_pval",
              "t_statistic",
              "comparison"
            )]) %>%
          dplyr::distinct() %>%
          dplyr::right_join(t_test_mean_sd_result, by = colnames(t_test_mean_sd_result)[!colnames(t_test_mean_sd_result) %in%
            c(
              "mean_control",
              "mean_treated",
              "sd_control",
              "sd_treated",
              "n_control",
              "n_treated",
              "pval",
              "std_error",
              "diff",
              "adj_pval",
              "t_statistic",
              "comparison"
            )]) %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
      }
      return(t_test_mean_sd_result)
    }

    if (method == "moderated_t-test") {
      if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package \"limma\" is needed for this function to work. Please install it.", call. = FALSE)
      }

      conditions_no_ref <- unique(pull(data, {{ condition }}))[!unique(pull(data, {{ condition }})) %in% ref_condition]

      message("[1/7] Creating moderated t-test input data ... ", appendLF = FALSE)
      moderated_t_test_input <- data %>%
        dplyr::distinct({{ grouping }}, {{ sample }}, {{ intensity_log2 }}) %>%
        tidyr::drop_na({{ intensity_log2 }}) %>%
        dplyr::arrange({{ sample }}) %>%
        tidyr::pivot_wider(names_from = {{ sample }}, values_from = {{ intensity_log2 }}) %>%
        tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
        as.matrix()

      message("DONE", appendLF = TRUE)
      message("[2/7] Defining moderated t-test design ... ", appendLF = FALSE)

      moderated_t_test_map <- data %>%
        dplyr::distinct({{ sample }}, {{ condition }}) %>%
        dplyr::mutate({{ condition }} := paste0("x", {{ condition }})) %>%
        dplyr::arrange({{ sample }})

      moderated_t_test_design <- stats::model.matrix(~ 0 + factor(
        stringr::str_replace_all(
          dplyr::pull(moderated_t_test_map, {{ condition }}),
          pattern = " ",
          replacement = "_"
        )
      ))

      colnames(moderated_t_test_design) <- levels(factor(
        stringr::str_replace_all(
          dplyr::pull(
            moderated_t_test_map,
            {{ condition }}
          ),
          pattern = " ",
          replacement = "_"
        )
      ))

      message("DONE", appendLF = TRUE)
      message("[3/7] Fitting lmFit model ... ", appendLF = FALSE)

      moderated_t_test_fit <- suppressWarnings(limma::lmFit(moderated_t_test_input, moderated_t_test_design))

      message("DONE", appendLF = TRUE)
      message("[4/7] Construct matrix of custom contrasts ... ", appendLF = FALSE)

      names <- paste0(
        "x", stringr::str_extract(unique(dplyr::pull(data, {{ comparison }})), pattern = ".+(?=_vs_)"), "_vs_x",
        stringr::str_extract(unique(dplyr::pull(data, {{ comparison }})), pattern = "(?<=_vs_).+")
      )

      comparisons <- paste0(
        "x", stringr::str_extract(
          stringr::str_replace_all(
            unique(dplyr::pull(
              data,
              {{ comparison }}
            )),
            pattern = " ",
            replacement = "_"
          ),
          pattern = ".+(?=_vs_)"
        ),
        "-x", stringr::str_extract(
          stringr::str_replace_all(
            unique(dplyr::pull(
              data,
              {{ comparison }}
            )),
            pattern = " ",
            replacement = "_"
          ),
          pattern = "(?<=_vs_).+"
        )
      )

      combinations <- purrr::map2(
        .x = names,
        .y = comparisons,
        .f = function(x, y) {
          rlang::exprs(!!rlang::as_name(x) := !!y)
        }
      )

      contrast_matrix <- eval(rlang::expr(limma::makeContrasts(!!!unlist(combinations), levels = moderated_t_test_design)))

      message("DONE", appendLF = TRUE)
      message("[5/7] Compute contrasts from linear model fit ... ", appendLF = FALSE)

      moderated_t_test_fit2 <- limma::contrasts.fit(moderated_t_test_fit, contrast_matrix)

      message("DONE", appendLF = TRUE)
      message("[6/7] Compute empirical Bayes statistics ... ", appendLF = FALSE)

      moderated_t_test_fit3 <- limma::eBayes(moderated_t_test_fit2)

      message("DONE", appendLF = TRUE)
      message("[7/7] Create result table ... ", appendLF = FALSE)

      moderated_t_test_missingness <-
        data %>%
        tidyr::drop_na({{ missingness }}, {{ intensity_log2 }}) %>%
        dplyr::group_by({{ comparison }}, {{ grouping }}) %>%
        dplyr::mutate(n_obs = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::distinct({{ grouping }}, {{ comparison }}, {{ missingness }}, .data$n_obs)

      moderated_t_test_result <- purrr::map_dfr(
        .x = names,
        .f = ~ limma::topTable(moderated_t_test_fit3,
          coef = .x,
          number = Inf,
          confint = TRUE,
          sort.by = "p",
          adjust.method = p_adj_method
        ) %>%
          tibble::rownames_to_column(rlang::as_name(rlang::enquo(grouping))) %>%
          dplyr::mutate(comparison = .x)
      ) %>%
        dplyr::mutate(comparison = stringr::str_replace_all({{ comparison }}, pattern = "^x|(?<=_vs_)x", replacement = "")) %>%
        dplyr::rename(
          diff = .data$logFC,
          CI_2.5 = .data$CI.L,
          CI_97.5 = .data$CI.R,
          t_statistic = .data$t,
          avg_abundance = .data$AveExpr,
          pval = .data$P.Value,
          adj_pval = .data$adj.P.Val
        ) %>%
        dplyr::left_join(moderated_t_test_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))

      message("DONE", appendLF = TRUE)

      if (!missing(retain_columns)) {
        moderated_t_test_result <- data %>%
          dplyr::ungroup() %>%
          dplyr::select(
            !!enquo(retain_columns),
            {{ intensity_log2 }},
            colnames(moderated_t_test_result)[!colnames(moderated_t_test_result) %in%
              c(
                "CI_2.5",
                "CI_97.5",
                "avg_abundance",
                "pval",
                "diff",
                "adj_pval",
                "t_statistic",
                "B",
                "n_obs"
              )]
          ) %>%
          tidyr::drop_na({{ intensity_log2 }}) %>%
          dplyr::select(-{{ intensity_log2 }}) %>%
          dplyr::distinct() %>%
          dplyr::right_join(moderated_t_test_result, by = colnames(moderated_t_test_result)[!colnames(moderated_t_test_result) %in%
            c(
              "CI_2.5",
              "CI_97.5",
              "avg_abundance",
              "pval",
              "diff",
              "adj_pval",
              "t_statistic",
              "B",
              "n_obs"
            )]) %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
      }

      if (filter_NA_missingness == TRUE) {
        moderated_t_test_result <- moderated_t_test_result %>%
          tidyr::drop_na({{ missingness }}) %>%
          dplyr::group_by(.data$comparison) %>%
          dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method)) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
        return(moderated_t_test_result)
      }
      if (filter_NA_missingness == FALSE) {
        return(moderated_t_test_result)
      }
    }

    if (method == "proDA") {
      if (!requireNamespace("proDA", quietly = TRUE)) {
        stop("Package \"proDA\" is needed for this function to work. Please install it.", call. = FALSE)
      }
      message("[1/5] Creating proDA input data ... ", appendLF = FALSE)

      proDA_input <- data %>%
        dplyr::distinct({{ grouping }}, {{ sample }}, {{ intensity_log2 }}) %>%
        dplyr::arrange({{ sample }}) %>%
        tidyr::pivot_wider(names_from = {{ sample }}, values_from = {{ intensity_log2 }}) %>%
        tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
        as.matrix()

      message("DONE", appendLF = TRUE)
      message("[2/5] Defining proDA design ... ", appendLF = FALSE)

      proDA_map <- data %>%
        dplyr::distinct({{ sample }}, {{ condition }}) %>%
        dplyr::arrange({{ sample }})

      proDA_design <- paste0("x", stringr::str_replace_all(dplyr::pull(proDA_map, {{ condition }}), pattern = " ", replacement = "_"))

      message("DONE", appendLF = TRUE)
      message("[3/5] Fitting proDA model (can take a few minutes) ... ", appendLF = FALSE)

      proDA_fit <-
        proDA::proDA(proDA_input,
          design = proDA_design
        )

      message("DONE", appendLF = TRUE)
      message("[4/5] Define missingness levels for filtering ... ", appendLF = FALSE)

      proDA_missingness <-
        data %>%
        tidyr::drop_na({{ missingness }}, {{ intensity_log2 }}) %>%
        dplyr::group_by({{ comparison }}, {{ grouping }}) %>%
        dplyr::mutate(n_obs = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::distinct({{ grouping }}, {{ comparison }}, {{ missingness }}, .data$n_obs)

      proDA_filter <- proDA_missingness %>%
        dplyr::distinct({{ grouping }}, {{ comparison }}) %>%
        split(.$comparison) %>%
        purrr::map(dplyr::select, -{{ comparison }})

      message("DONE", appendLF = TRUE)
      message("[5/5] Extracting differential abundance from model and apply filter ... ", appendLF = FALSE)

      names <- unique(dplyr::pull(data, {{ comparison }}))
      comparisons <- paste0(
        "x",
        stringr::str_extract(
          stringr::str_replace_all(
            unique(dplyr::pull(
              data,
              {{ comparison }}
            )),
            pattern = " ",
            replacement = "_"
          ),
          pattern = ".+(?=_vs_)"
        ), " - x",
        stringr::str_extract(
          stringr::str_replace_all(
            unique(dplyr::pull(
              data,
              {{ comparison }}
            )),
            pattern = " ",
            replacement = "_"
          ),
          pattern = "(?<=_vs_).+"
        )
      )

      proDA_result <- names %>%
        purrr::map(~proDA_fit) %>%
        purrr::set_names(nm = names) %>%
        purrr::map2(
          .y = comparisons,
          .f = ~ proDA::test_diff(.x, contrast = .y, sort_by = "adj_pval")
        )

      if (filter_NA_missingness == TRUE) {
        proDA_result <- proDA_result %>%
          purrr::map2(
            .y = proDA_filter,
            .f = ~ dplyr::inner_join(.x, .y, by = c("name" = as_label(enquo(grouping))))
          ) %>%
          purrr::map2(
            .y = names(.),
            .f = ~ dplyr::mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))
          ) %>%
          purrr::map_dfr(~ dplyr::mutate(.x, adj_pval = p.adjust(.data$pval, method = p_adj_method))) %>%
          dplyr::select(-.data$n_obs, -.data$n_approx) %>%
          dplyr::rename({{ grouping }} := .data$name, std_error = .data$se) %>%
          dplyr::left_join(proDA_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))

        message("DONE", appendLF = TRUE)
      }
      if (filter_NA_missingness == FALSE) {
        proDA_result <- proDA_result %>%
          purrr::map2(
            .y = names(.),
            .f = ~ dplyr::mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))
          ) %>%
          purrr::map_dfr(~ dplyr::mutate(.x, adj_pval = p.adjust(.data$pval, method = p_adj_method))) %>%
          dplyr::select(-.data$n_obs) %>%
          dplyr::select(-.data$n_obs, -.data$n_approx) %>%
          dplyr::rename({{ grouping }} := .data$name, std_error = .data$se) %>%
          dplyr::left_join(proDA_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))

        message("DONE", appendLF = TRUE)
      }

      if (!missing(retain_columns)) {
        proDA_result <- data %>%
          dplyr::ungroup() %>%
          dplyr::select(
            !!enquo(retain_columns),
            {{ intensity_log2 }},
            colnames(proDA_result)[!colnames(proDA_result) %in%
              c(
                "std_error",
                "avg_abundance",
                "pval",
                "diff",
                "adj_pval",
                "t_statistic",
                "df",
                "n_obs"
              )]
          ) %>%
          tidyr::drop_na({{ intensity_log2 }}) %>%
          dplyr::select(-{{ intensity_log2 }}) %>%
          dplyr::distinct() %>%
          dplyr::right_join(proDA_result, by = colnames(proDA_result)[!colnames(proDA_result) %in%
            c(
              "std_error",
              "avg_abundance",
              "pval",
              "diff",
              "adj_pval",
              "t_statistic",
              "df",
              "n_obs"
            )]) %>%
          dplyr::arrange(.data$adj_pval, .data$pval)
      }
      return(proDA_result)
    }
  }
