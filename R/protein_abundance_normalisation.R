#' Perform normalization of LiP peptide data by protein abundance from tryptic controls
#'
#' LiP peptide intensities are normalised by protein abundances from tryptic controls to correct for protein fold changes due to treatment. 
#'
#' @param lip_peptides A data frame containing LiP peptide level data. The \code{sample, protein_id, grouping, condition, peptide_intensity} variables 
#' are required. Ideally the output of \code{assign_missingness} or \code{impute} is used for the method \code{"ref_vs_rest"}, since here also the variable 
#' \code{missingness} is required in contrast to method \code{"all_vs_all"}.
#' @param tc_proteins A data frame containing tryptic control protein level data. The \code{sample, protein_id, condition, protein_intensity} variables
#' are required. Ideally the output of \code{assign_missingness} or \code{impute} is used for the method \code{"ref_vs_rest"}, since here also the variable 
#' \code{missingness} is required in contrast to method \code{"all_vs_all"}.
#' @param sample The column in the data frame containing the sample name. Is not required if \code{method = "t-test_mean_sd"}.
#' @param protein_id The column in the data frame containing the protein accession numbers. 
#' @param grouping The column in the data frame containing precursor or peptide identifiers.
#' @param condition The column in the data frame containing the conditions.
#' @param peptide_intensity The column in the data frame containing peptide intensity values. The intensity values need to be log2 transformed.
#' @param protein_intensity The column in the data frame containing protein intensity values. The intensity values need to be log2 transformed.
#' @param method A character vector specifying which method should be used. If \code{method = "all_vs_all"} an ANOVA based method will be performed that 
#' compares all conditions with each other followed by Tukey's honestly significant difference test to identify pairs of conditions that are significantly different.
#' If \code{mehtod = "ref_vs_rest"} a reference condition will be compared against all other conditions using a Welch t-test in order to identify significantly
#' changing peptides or precursors.
#' @param paired A logical indicating if LiP and tryptic control data are paired. In other words it indicates if replicates of LiP and tryptic control belong together.
#' @param replicate_index The column in the data frame containing numbering for condition replicates. This column is only required if LiP and tryptic control samples are paired.
#' In that case the corresponding samples in both files should have the same numbering.
#' @param comparison The column in the data frame containing comparison information of treatment/reference condition pairs. Can be obtained by 
#' calling \code{assign_missingness}. 
#' @param completeness The minimal degree of data completeness to be considered for evaluation. Value has to be between 0 and 1, default is 0.7. It is multiplied 
#' with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each condition to be considered for
#' evaluation.
#'
#' @return A data frame that contains significance test statistics and differential abundances (\code{diff}) for each peptide or precursor and the 
#' corresponding comparison of conditions. If \code{method = "all_vs_all"} was used, ANOVA statistics are reported (\code{ms_group, ms_error, f, pval, adj_pval, n_conditions}). 
#' \code{ms_group} is the error within a group and \code{ms_error} is the error between groups. Also Tukey's honestly significant difference test statistics are 
#' reported (\code{df, q_value, pval_tukey}). If \code{method = "ref_vs_rest"} Welch t-test statistics are reported (\code{std_error, t_statistic, pval, adj_pval}).
#' @import dplyr
#' @import tidyr
#' @importFrom stats sd
#' @importFrom naniar replace_with_na
#' @importFrom rlang .data enquo ensym as_name expr exprs := !!
#' @importFrom purrr map_df 
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' \dontrun{
#' protein_abundance_normalisation(
#' lip_peptides = data_lip, 
#' tc_proteins = data_tc, 
#' sample = r_file_name,
#' protein_id = pg_protein_accessions, 
#' grouping = eg_precursor_id, 
#' condition = r_condition,
#' peptide_intensity = normalised_intensity_log2,
#' protein_intensity = protein_intensity_log2,
#' method = "all_vs_all")
#' }
protein_abundance_normalisation <- function(lip_peptides, tc_proteins, sample, grouping, condition, protein_id, peptide_intensity, protein_intensity, method, paired = FALSE, replicate_index = NULL, comparison = NULL, completeness = 0.7){
  .=NULL
  # keep copy of input data
  input_lip_peptides <- lip_peptides
  
  message("General data preparation ... ", appendLF = FALSE)
  
  # filter data to only include proteins identified in both datasets
  
  proteins_lip <- unique(dplyr::pull(lip_peptides, {{protein_id}}))
  proteins_tc <- unique(dplyr::pull(tc_proteins, {{protein_id}}))
  
  proteins_overlap <- dplyr::intersect(proteins_lip, proteins_tc)
  
  # subset both datasets based on overlapping proteins
  
  lip_peptides <- lip_peptides %>% 
    dplyr::filter({{protein_id}} %in% proteins_overlap) %>% 
    dplyr::distinct({{grouping}}, {{condition}}, {{protein_id}}, {{peptide_intensity}}, {{comparison}}, {{replicate_index}})

  tc_proteins <- tc_proteins %>% 
    dplyr::filter({{protein_id}} %in% proteins_overlap) %>% 
    dplyr::distinct({{protein_id}}, {{condition}}, {{protein_intensity}}, {{comparison}}, {{replicate_index}})
  
  if(paired == FALSE){
  # calculate mean, sd and n_replicates
  
  lip_peptides <- lip_peptides %>% 
    naniar::replace_with_na(replace = eval(rlang::expr(list(!!!rlang::exprs(!!rlang::as_name(rlang::enquo(peptide_intensity)) := 0))))) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by({{grouping}}, {{condition}}, {{comparison}}) %>% 
    dplyr::mutate(n_replicates_peptide = dplyr::n(),
                  mean_peptide = mean({{peptide_intensity}}),
                  sd_peptide = stats::sd({{peptide_intensity}})) %>% 
    dplyr::select(-c({{peptide_intensity}})) %>% 
    dplyr::distinct()
  
  tc_proteins <- tc_proteins %>% 
    naniar::replace_with_na(replace = eval(rlang::expr(list(!!!rlang::exprs(!!rlang::as_name(rlang::enquo(protein_intensity)) := 0))))) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by({{protein_id}}, {{condition}}, {{comparison}}) %>%
    dplyr::mutate(n_replicates_protein = dplyr::n(),
                  mean_protein = mean({{protein_intensity}}),
                  sd_protein = stats::sd({{protein_intensity}})) %>% 
    dplyr::select(-c({{protein_intensity}})) %>% 
    dplyr::distinct()
  
  # combine data frames and calculate peptide/protein ratios and propagated standard deviation
  
  combined <- lip_peptides %>% 
    dplyr::left_join(tc_proteins, by = c(rlang::as_name(rlang::enquo(protein_id)), rlang::as_name(rlang::enquo(condition)), rlang::as_name(rlang::enquo(comparison)))) %>% 
    dplyr::mutate(peptide_protein_ratio = .data$mean_peptide - .data$mean_protein,
                  peptide_protein_sd = sqrt(.data$sd_peptide^2 + .data$sd_protein^2))
  }
  
  if(paired == TRUE){
    #combine lip_peptides and tc_proteins dataframes based on replicate index, condition and peptide. Calculate mean, sd and n_replicates
    
    lip_peptides <- lip_peptides %>% 
      naniar::replace_with_na(replace = eval(rlang::expr(list(!!!rlang::exprs(!!rlang::as_name(rlang::enquo(peptide_intensity)) := 0))))) %>% 
      tidyr::drop_na()
    
    tc_proteins <- tc_proteins %>% 
      naniar::replace_with_na(replace = eval(rlang::expr(list(!!!rlang::exprs(!!rlang::as_name(rlang::enquo(protein_intensity)) := 0))))) %>% 
      tidyr::drop_na()
    
    combined <- lip_peptides %>% 
      dplyr::left_join(tc_proteins, by = c(rlang::as_name(rlang::enquo(protein_id)), rlang::as_name(rlang::enquo(condition)), rlang::as_name(rlang::enquo(replicate_index)))) %>% #if only peptides for which a protein intensity was measured should be included change to inner_join
      dplyr::group_by({{grouping}}, {{condition}}, {{comparison}}) %>% 
      dplyr::mutate(n_replicates_peptide = dplyr::n(),
                    mean_peptide = mean({{peptide_intensity}}),
                    sd_peptide = stats::sd({{peptide_intensity}})) %>% 
      tidyr::drop_na({{protein_intensity}}) %>% # to remove peptide intensities that do not have associated protein intensities. This is important for proper n_replicates_protein calculation.
      dplyr::mutate(n_replicates_protein = dplyr::n(),
                    mean_protein = mean({{protein_intensity}}),
                    sd_protein = stats::sd({{protein_intensity}})) %>% 
      dplyr::select(-c({{protein_intensity}}, {{peptide_intensity}}, {{replicate_index}})) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(peptide_protein_ratio = .data$mean_peptide - .data$mean_protein,
                    peptide_protein_sd = sqrt(.data$sd_peptide^2 + .data$sd_protein^2))
  }
  
  message("DONE", appendLF = TRUE)
  
  if(method == "all_vs_all"){
    
    message("[1/4] Prepare data for ANOVA ... ", appendLF = FALSE)
    
    #calculate the minimal number of replicates that should be present for a sample to be kept
    n_all_conditions <- length(unique(dplyr::pull(input_lip_peptides, {{condition}})))
    n_replicates <- nrow(dplyr::distinct(input_lip_peptides, {{sample}}, {{condition}})) / n_all_conditions
    min_replicates <- floor(n_replicates * completeness)

    #prepare input data for anova test
    anova_input <- combined %>%
      dplyr::filter(.data$n_replicates_peptide >= min_replicates) %>%
      tidyr::drop_na() %>% 
      dplyr::group_by({{grouping}}) %>% 
      dplyr::mutate(n_conditions = dplyr::n_distinct(!!rlang::ensym(condition))) %>% 
      dplyr::filter(.data$n_conditions >= 2)
    
    message("DONE", appendLF = TRUE)
    message("[2/4] Perform ANOVA ... ", appendLF = FALSE)
    
    #perform anova test
    anova_output <- anova_protti(data = anova_input,
                                 grouping = {{grouping}},
                                 condition = {{condition}},
                                 mean_ratio = .data$peptide_protein_ratio,
                                 sd = .data$peptide_protein_sd,
                                 n = .data$n_replicates_peptide) %>% 
      dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>% 
      dplyr::left_join(anova_input %>% 
                         dplyr::distinct({{condition}}, {{grouping}}, {{protein_id}}, .data$peptide_protein_ratio, .data$n_replicates_peptide, .data$n_conditions), 
                       by = rlang::as_name(rlang::enquo(grouping)))
    
    message("DONE", appendLF = TRUE)
    message("[3/4] Prepare data for Tukey's honestly significant difference test ... ", appendLF = FALSE)
    
    #prepare data for tukey's honestly significant difference test
    
    conditions <- unique(dplyr::pull(input_lip_peptides, {{condition}}))
    combinations <- suppressWarnings(tibble::as_tibble(t(utils::combn(conditions, 2)))) %>% 
      dplyr::mutate(comparisons = paste0(.data$V1, "_vs_", .data$V2)) %>% 
      tidyr::pivot_longer(cols = c(.data$V1, .data$V2), values_to = rlang::as_name(rlang::enquo(condition))) %>% 
      dplyr::select(-.data$name)
    
    input_tukey <- anova_output %>% 
      dplyr::group_by({{grouping}}) %>% 
      dplyr::mutate(df = sum(.data$n_replicates_peptide) - .data$n_conditions) %>% 
      dplyr::select(-.data$n_replicates_peptide) %>% 
      dplyr::left_join(combinations, by = rlang::as_name(rlang::enquo(condition))) %>% 
      split(.$comparisons) 
    
    message("DONE", appendLF = TRUE)
    
    #perform tukey's honestly significant difference test
    
    pb <- progress::progress_bar$new(total = length(input_tukey), format = "[4/4] Perform Tukey's honestly significant difference test [:bar] :current/:total (:percent) :eta")
    
    output_tukey <- input_tukey %>% 
      purrr::map_df(.f = ~{pb$tick(); .x %>%
                      dplyr::mutate({{condition}} := ifelse({{condition}} == unique(stringr::str_extract(.data$comparisons, pattern = ".+(?=_vs_)")), "condition_1", "condition_2")) %>%
                      tidyr::pivot_wider(names_from = {{condition}}, values_from = c(.data$peptide_protein_ratio)) %>% 
                      dplyr::mutate(q_value = (.data$condition_1 - .data$condition_2) / sqrt(.data$ms_error / .data$n_conditions)) %>% 
                      dplyr::mutate(diff = .data$condition_1 - .data$condition_2) %>% 
                      dplyr::mutate(pval_tukey = ptukey(abs(.data$q_value), nmeans = n_replicates, df = .data$df, lower.tail = FALSE))
                    }) %>% 
      dplyr::arrange(.data$pval_tukey)
    
    return(output_tukey)
  }
  
  if(method == "ref_vs_rest"){
    #extract reference condition
    ref_condition = stringr::str_extract(dplyr::pull(input_lip_peptides, {{comparison}})[1], pattern = "(?<=_vs_).+")
    
    message("[1/2] Prepare data for t-test ... ", appendLF = FALSE)
    
    ttest_input <- combined %>% 
      tidyr::drop_na() %>% 
      dplyr::group_by({{grouping}}, {{comparison}}) %>% 
      dplyr::mutate(n_conditions = dplyr::n_distinct({{condition}})) %>% 
      dplyr::filter(.data$n_conditions == 2) %>% 
      dplyr::distinct({{condition}}, {{comparison}}, {{grouping}}, {{protein_id}}, .data$n_replicates_peptide, .data$peptide_protein_ratio, .data$peptide_protein_sd) %>% 
      split(dplyr::pull(., {{comparison}})) 
    
    message("DONE", appendLF = TRUE)
    message("[2/2] Perform t-test ... ", appendLF = FALSE)
    
    ttest_output <- ttest_input %>% 
      purrr::map_df(.f = ~ .x %>%
                        dplyr::mutate({{condition}} := ifelse({{condition}} == ref_condition, "control", "treated")) %>%
                        tidyr::pivot_wider(names_from = {{condition}}, values_from = c(.data$peptide_protein_ratio, .data$peptide_protein_sd, .data$n_replicates_peptide)) %>% 
                        dplyr::mutate(ttest_protti(.data$peptide_protein_ratio_treated, 
                                                            .data$peptide_protein_ratio_control, 
                                                            .data$peptide_protein_sd_control, 
                                                            .data$peptide_protein_sd_treated, 
                                                            .data$n_replicates_peptide_control, 
                                                            .data$n_replicates_peptide_treated))
                    ) %>% 
      tidyr::drop_na(.data$pval) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(adj_pval = stats::p.adjust(.data$pval, method = "BH")) %>% 
      dplyr::arrange(.data$pval)
    
    message("DONE", appendLF = TRUE)
    
    return(ttest_output)
  }
}
# 
# 
# check2 <- protein_abundance_normalisation(lip_peptides = missingness_lip, 
#                               tc_proteins = tc_prot, 
#                               sample = r_file_name,
#                               protein_id = pg_protein_accessions, 
#                               grouping = eg_precursor_id, 
#                               condition = r_condition,
#                               peptide_intensity = normalised_intensity_log2,
#                               protein_intensity = protein_intensity,
#                               method = "ref_vs_rest",
#                               comparison = comparison,
#                               replicate_index = replicate_index,
#                               paired = FALSE)
