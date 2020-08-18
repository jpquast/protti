#' Calculate differential abundance between conditions
#'
#' Performs differential abundance calculations and statistical testing on data frames with protein, peptide or precursor data. Different methods for statistical testin are available.
#'
#' @param data A data frame containing at least the input variables that are required for the selected method. Ideally the output of \code{assign_missingness} or \code{impute} is used.
#' @param sample The column in the data frame containing the sample name. Is not required if \code{method = "t-test_mean_sd"}.
#' @param condition The column in the data frame containing the conditions.
#' @param grouping The column in the data frame containing precursor or peptide identifiers.
#' @param intensity The column in the data frame containing intensity values. The intensity values need to be log2 transformed. Is not required if \code{method = "t-test_mean_sd"}.
#' @param missingness The column in the data frame containing missingness information. Can be obtained by calling \code{assign_missingness}. 
#' Is not required if \code{method = "t-test_mean_sd"}.
#' @param comparison The column in the data frame containing comparison information of treatment/reference condition pairs. Can be obtained by 
#' calling \code{assign_missingness}. Is not required if \code{method = "t-test_mean_sd"}.
#' @param mean The column in the data frame containing mean values for two conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param sd The column in the data frame containing standard deviations for two conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param n_samples The column in the data frame containing the number of samples per condition for two conditions. Is only required if \code{method = "t-test_mean_sd"}.
#' @param ref_condition The condition that is used as a reference for differential abundance calculation.
#' @param filter_NA_missingness A logical, default is \code{TRUE}. For all methods except \code{"t-test_mean_sd"} missingness information provided. 
#' If a reference/treatment pair has too few samples to be considered robust, it is annotated with \code{NA} as missingness. If this argument 
#' is \code{TRUE}, these reference/treatment pairs are filtered out. 
#' @param method A character vector, specifies the method used for statistical testing. Methods include Welch test ("\code{t-test}"), a Welch test on means, 
#' standard deviations and number of replicates ("\code{t-test_mean_sd}") and a moderated t-test based on the \code{limma} package ("\code{moderated_t-test}"). 
#' More information on the moderated t-test can be found in the \code{limma} documentation. Furthermore, the \code{proDA} package specific method ("\code{proDA}") can 
#' be used to infer means across samples based on a probabilistic dropout model. This eliminates the need for data imputation since missing values are infered from the 
#' model. More information can be found in the \code{proDA} documentation.
#' @param p_adj_method A character vector, specifies the p-value correction method. Possible methods are c("holm", "hochberg", "hommel", "bonferroni", "BH", 
#' "BY", "fdr", "none"). Default method is \code{"BH"}.
#'
#' @return A data frame that contains differential abundances (\code{diff}), p-values (\code{pval}) and adjusted p-values (\code{adj_pval}) for each protein, 
#' peptide or precursor (depending on the \code{grouping} variable) and the associated treatment/reference pair.
#' Depending on the method the data frame contains additional columns:
#' \itemize{
#' \item{"t-test": }{The \code{std_error} column contains the standard error of the differential abundances. \code{n_obs} contains the number of 
#' observations for the specific protein, peptide or precursor (depending on the \code{grouping} variable) and the associated treatment/reference pair.}
#' \item{"t-test_mean_sd": }{\code{mean_control} and \code{mean_treated} columns contain the means for the reference and treatment condition, respectively. 
#' \code{sd_control} and \code{sd_treated} columns contain the standard deviations for the reference and treatment condition, respectively. 
#' \code{n_control} and \code{n_treated} columns contain the numbers of samples for the reference and treatment condition, respectively. The \code{std_error} 
#' column contains the standard error of the differential abundances. \code{t_statistic} contains the t_statistic for the t-test.}
#' \item{"moderated_t-test": }{\code{CI_2.5} and \code{CI_97.5} give the 2.5% and 97.5% confidence interval borders for the differential abundance. \code{avg_abundance} 
#' contains average abundances for treatment/reference pairs (mean of the two group means). \code{t_statistic} contains the t_statistic for the t-test. \code{B} The 
#' B-statistic is the log-odds that the protein, peptide or precursor (depending on \code{grouping}) has a differential abundance between the two groups. Suppose B=1.5. 
#' The odds of differential abundance is exp(1.5)=4.48, i.e, about four and a half to one. The probability that there is a differentail abundance is 4.48/(1+4.48)=0.82, 
#' i.e., the probability is about 82% that this group is differentially abundant. A B-statistic of zero corresponds to a 50-50 chance that the group is differentially 
#' abundant.\code{n_obs} contains the number of observations for the specific protein, peptide or precursor (depending on the \code{grouping} variable) and the 
#' associated treatment/reference pair.}
#' \item{"proDA": }{The \code{std_error} column contains the standard error of the differential abundances. \code{avg_abundance} contains average abundances for 
#' treatment/reference pairs (mean of the two group means). \code{t_statistic} contains the t_statistic for the t-test. \code{n_obs} contains the number of 
#' observations for the specific protein, peptide or precursor (depending on the \code{grouping} variable) and the associated treatment/reference pair.}
#' } 
#' @import dplyr
#' @import tidyr
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom proDA result_names proDA test_diff
#' @importFrom rlang .data enquo ensym as_name as_label expr := !!
#' @importFrom purrr map map2 map_df map_dbl map_chr map2_dbl reduce set_names
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_replace_all
#' @export
#'
#' @examples
#' \dontrun{
#' diff_abundance(
#' data,
#' sample = r_file_name,
#' condition = r_condition,
#' grouping = eg_precursor_id,
#' intensity = normalised_intensity_log2,
#' missingness = missingness,
#' comparison = comparison,
#' ref_condition = "control",
#' method = "t-test"
#' )
#' }
diff_abundance <-
  function(data,
           sample,
           condition,
           grouping,
           intensity,
           missingness,
           comparison,
           mean = NULL,
           sd = NULL,
           n_samples = NULL,
           ref_condition,
           filter_NA_missingness = TRUE,
           method = c("t-test", "t-test_mean_sd", "moderated_t-test", "proDA"),
           p_adj_method = "BH") {
  . = NULL  
    
  method <- match.arg(method)
  
  if(method != "t-test_mean_sd") {
  if(max(pull(data, {{intensity}}), na.rm = TRUE) > 1000) {stop("Please log2 transform your intensities.")}
  }
  if(method == "t-test_mean_sd") {
    if(max(pull(data, {{mean}}), na.rm = TRUE) > 1000) {stop("Please log2 transform your data.")}
  }
  
  if(method == "t-test"){

    message("[1/2] Create input for t-tests ... ", appendLF = FALSE)
    
    t_test_missingness_obs <- data %>%
      tidyr::drop_na({{missingness}}, {{intensity}}) %>%
      dplyr::group_by({{comparison}}, {{grouping}}) %>%
      dplyr::mutate(n_obs = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::distinct({{grouping}}, {{comparison}}, {{missingness}}, .data$n_obs)
    
    t_test_input <- data %>%
      dplyr::group_by({{comparison}}, {{grouping}}, {{condition}}) %>%
      dplyr::summarize(intensity = list({{intensity}}), .groups = "drop") %>%
      dplyr::mutate(type = ifelse({{condition}} == "control_AE", "control", "treated")) %>%
      dplyr::select(-{{condition}}) %>%
      tidyr::pivot_wider(names_from = .data$type, values_from = .data$intensity)
    
    message("DONE", appendLF = TRUE)
    message("[2/2] Calculate t-tests ... ", appendLF = FALSE)
    
    t_test_result <- t_test_input %>%
      dplyr::mutate(pval = purrr::map2(
        .x = .data$treated,
        .y = .data$control,
        .f = function(.x, .y) {
          tryCatch({
            stats::t.test(.x, .y)
          },
          error = function(error) {
            NA
          })
        }
      )) %>%
      dplyr::mutate(std_error = purrr::map_dbl(
        .x = .data$pval,
        .f = ~ tryCatch({.x$stderr
        },
        error = function(error) {
          NA
        })
      )) %>%
      dplyr::mutate(pval = purrr::map_dbl(
        .x = .data$pval,
        .f = ~ tryCatch({.x$p.value
        },
        error = function(error) {
          NA
        })
      )) %>%
      dplyr::mutate(diff = map2_dbl(
        .x = .data$treated,
        .y = .data$control,
        .f = function(.x, .y) {
            mean(.x, na.rm = TRUE) - mean(.y, na.rm = TRUE)
        }
      )) %>%
      dplyr::mutate(diff = ifelse(diff == "NaN", NA, diff)) %>%
      dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method)) %>%
      dplyr::select(-c(.data$control, .data$treated)) %>%
      dplyr::left_join(t_test_missingness_obs, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
    
    message("DONE", appendLF = TRUE)
    
    if(filter_NA_missingness == TRUE){
      t_test_result <- t_test_result %>%
        tidyr::drop_na({{missingness}})
      return(t_test_result)
    }
    if(filter_NA_missingness == FALSE) {
      return(t_test_result)
    }
  }
  
  if(method == "t-test_mean_sd") {
    conditions_no_ref <- unique(dplyr::pull(data, {{condition}}))[!unique(dplyr::pull(data, {{condition}})) %in% ref_condition]
     t_test_mean_sd_result <- data %>%
       dplyr::mutate(comparison = ifelse(
         {{condition}} %in% conditions_no_ref,
         paste0({{condition}}, "_vs_", ref_condition),
         list(paste0(conditions_no_ref, "_vs_",  ref_condition))
       )) %>%
       tidyr::unnest(.data$comparison) %>%
       dplyr::rename(mean = {{mean}}, sd = {{sd}}, n = {{n_samples}}) %>%
       split(.$comparison) %>%
       purrr::map_dfr( ~ .x %>%
                dplyr::mutate({{condition}} := ifelse({{condition}} == ref_condition, "control", "treated")) %>%
                tidyr::pivot_wider(names_from = {{condition}}, values_from = c(.data$mean, .data$sd, .data$n)) %>%
                  dplyr::mutate(ttest_protti(.data$mean_treated, .data$mean_control, .data$sd_control, .data$sd_treated, .data$n_control, .data$n_treated))
               ) %>%
       tidyr::drop_na(.data$pval) %>%
       dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = p_adj_method))

    return(t_test_mean_sd_result)
  }
  
  if(method == "moderated_t-test"){
  conditions_no_ref <- unique(pull(data, {{condition}}))[!unique(pull(data, {{condition}})) %in% ref_condition]
  
  message("[1/7] Creating moderated t-test input data ... ", appendLF = FALSE)
  moderated_t_test_input <- data %>%
    dplyr::distinct({{grouping}}, {{sample}}, {{intensity}}) %>%
    dplyr::arrange({{sample}}) %>% 
    tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
    as.matrix()   
    
  message("DONE", appendLF = TRUE)
  message("[2/7] Defining moderated t-test design ... ", appendLF = FALSE)
  
  moderated_t_test_map <- data %>%
    dplyr::distinct({{sample}}, {{condition}}) %>% 
    dplyr::mutate({{condition}} := paste0("x", {{condition}})) %>%
    dplyr::arrange({{sample}})
  
  moderated_t_test_design <- stats::model.matrix(~ 0 + factor(dplyr::pull(moderated_t_test_map, {{condition}})))
  
  colnames(moderated_t_test_design) <- levels(factor(dplyr::pull(moderated_t_test_map, {{condition}})))
  
  message("DONE", appendLF = TRUE)
  message("[3/7] Fitting lmFit model ... ", appendLF = FALSE)
  
  moderated_t_test_fit <- suppressWarnings(limma::lmFit(moderated_t_test_input, moderated_t_test_design))
  
  message("DONE", appendLF = TRUE)
  message("[4/7] Construc matrix of custom contrasts ... ", appendLF = FALSE)
  
  names <- purrr::map_chr(conditions_no_ref, ~paste0("x", .x, "_vs_x", ref_condition))
  comparisons <- purrr::map_chr(conditions_no_ref, ~paste0("x", .x, "-x", ref_condition))
  combinations <- purrr::map2(.x = names,
                      .y = comparisons,
                      .f = function(x, y){
                         rlang::exprs(!!rlang::as_name(x) := !!y)
                       })
  
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
    tidyr::drop_na({{missingness}}, {{intensity}}) %>%
    dplyr::group_by({{comparison}}, {{grouping}}) %>%
    dplyr::mutate(n_obs = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct({{grouping}}, {{comparison}}, {{missingness}}, .data$n_obs)
  
  moderated_t_test_result <- purrr::map_dfr(.x = names,
                .f = ~ limma::topTable(moderated_t_test_fit3, coef = .x, number = Inf, confint = TRUE, sort.by = "p", adjust.method = p_adj_method) %>%
                  tibble::rownames_to_column(rlang::as_name(rlang::enquo(grouping))) %>%
                  dplyr::mutate(comparison = .x)
  ) %>%
    dplyr::mutate(comparison = stringr::str_replace_all({{comparison}}, pattern = "^x|(?<=_vs_)x", replacement = "")) %>%
    dplyr::rename(diff = .data$logFC, CI_2.5 = .data$CI.L, CI_97.5 = .data$CI.R, t_statistic = .data$t, avg_abundance = .data$AveExpr, pval = .data$P.Value, adj_pval = .data$adj.P.Val) %>%
    dplyr::left_join(moderated_t_test_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
  
  message("DONE", appendLF = TRUE)
  
  if(filter_NA_missingness == TRUE){
    moderated_t_test_result <- moderated_t_test_result %>%
      tidyr::drop_na({{missingness}})
    return(moderated_t_test_result)
  }
  if(filter_NA_missingness == FALSE) {
    return(moderated_t_test_result)
  }
  }
  
  if(method == "proDA"){
  message("[1/5] Creating proDA input data ... ", appendLF = FALSE)
    
  proDA_input <- data %>%
    dplyr::distinct({{grouping}}, {{sample}}, {{intensity}}) %>%
    dplyr::arrange({{sample}}) %>% 
    tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(rlang::enquo(grouping))) %>%
    as.matrix()
  
  message("DONE", appendLF = TRUE)
  message("[2/5] Defining proDA design ... ", appendLF = FALSE)
  
  proDA_map <- data %>%
    dplyr::distinct({{sample}}, {{condition}}) %>% 
    dplyr::arrange({{sample}})
  
  proDA_design <- dplyr::pull(proDA_map, {{condition}})
  
  message("DONE", appendLF = TRUE)
  message("[3/5] Fitting proDA model (can take a few minutes) ... ", appendLF = FALSE)
  
  proDA_fit <- 
    proDA::proDA(proDA_input,
          design = proDA_design,
          reference_level = ref_condition)
  
  proDA_result_names <-
    proDA::result_names(proDA_fit)[!proDA::result_names(proDA_fit) %in% c("Intercept")]
  
  message("DONE", appendLF = TRUE)
  message("[4/5] Define missingness levels for filtering ... ", appendLF = FALSE)
  
  proDA_missingness <- 
    data %>% 
    tidyr::drop_na({{missingness}}, {{intensity}}) %>%
    dplyr::group_by({{comparison}}, {{grouping}}) %>%
    dplyr::mutate(n_obs = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct({{grouping}}, {{comparison}}, {{missingness}}, .data$n_obs)
  
  proDA_filter <- proDA_missingness %>%
    dplyr::distinct({{grouping}}, {{comparison}}) %>% 
    split(.$comparison) %>%
    purrr::map(dplyr::select, -{{comparison}})
  
  message("DONE", appendLF = TRUE)
  message("[5/5] Extracting differential abundance from model and apply filter ... ", appendLF = FALSE)
  
  proDA_result <- proDA_result_names %>%
    purrr::map(~ proDA_fit) %>%
    purrr::set_names(nm = proDA_result_names) %>%
    purrr::map2(.y = names(.),
         .f = ~ proDA::test_diff(.x, contrast = .y, sort_by = "adj_pval")) 
  
  if(filter_NA_missingness == TRUE) {
  proDA_result <- proDA_result %>% 
    purrr::map2(.y = proDA_filter, 
         .f = ~ dplyr::inner_join(.x, .y, by = c("name" = as_label(enquo(grouping))))) %>%
    purrr::map2(.y = names(.), 
         .f = ~ dplyr::mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))) %>% 
    purrr::map_dfr(~ dplyr::mutate(.x, adj_pval = p.adjust(.data$pval, method = p_adj_method))) %>%
    dplyr::select(-.data$n_obs, -.data$n_approx) %>%
    dplyr::rename({{grouping}} := .data$name, std_error = .data$se) %>%
    dplyr::left_join(proDA_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
  
  message("DONE", appendLF = TRUE)
  
  return(proDA_result)
  }
  if(filter_NA_missingness == FALSE) {
    proDA_result <- proDA_result %>% 
      purrr::map2(.y = names(.), 
           .f = ~ dplyr::mutate(.x, comparison = str_replace_all(.y, pattern = "`", replacement = ""))) %>% 
      purrr::map_dfr(~ dplyr::mutate(.x, adj_pval = p.adjust(.data$pval, method = p_adj_method))) %>%
      dplyr::select(-.data$n_obs) %>%
      dplyr::select(-.data$n_obs, -.data$n_approx) %>%
      dplyr::rename({{grouping}} := .data$name, std_error = .data$se) %>%
      dplyr::left_join(proDA_missingness, by = c(rlang::as_name(rlang::enquo(grouping)), "comparison"))
    
    message("DONE", appendLF = TRUE)
    
    return(proDA_result)
  }
  }
}