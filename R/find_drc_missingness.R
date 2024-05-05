#' Find missing not at random cases in dose-response type data
#'
#' Dose-response data can suffer from loss of detection due to the detection limit
#' of the device used to record the data. On mass spectrometers there will always
#' be precursor ions detected at the lower end of the intensity normal distribution.
#' If one of these precursors responds to a dose treatment with a decrease in intensity,
#' it might fall below the limit of detection and therefore no more observations are 
#' made in the upper or lower end of the treatment concentrations. While it is not possible
#' to fit a sigmoidal dose-response curve and to extract useful parameters in these cases,
#' it is at least possible to classify them as responders. For that we are looking for 
#' dose-related consecutive not at random missingness (MNAR) of observations. 
#' 
#' We first cluster the data in two groups with the most unequal number of observations possible.
#' For that we optimize for the minimal summed squared error of two clusters for each possible dataset 
#' split along the dose axis. On this split we perform a fisher's exact test, for the enrichment of
#' detected and undetected samples. The obtained p-value leaves us with a good ranking of possible
#' missing not at random cases.
#' 
#' @param data a data frame containing at least the input columns.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' @param grouping a character column in the \code{data} data frame that contains e.g. protein, precursor or
#' peptide identifiers. Can also be any other identifiers as long as there is one potential dose-response 
#' curve per entity in this column.
#' @param response a numeric column in the \code{data} data frame that contains e.g. intensity values that
#' make up the response of the dose-response curve.
#' @param dose a numeric column in the \code{data} data frame that contains dose information. This does not 
#' need to be a concentration but can also simply be the dose order.
#' @param min_cluster_size a numeric value that specifies the minimal size of either of the two clusters
#' the dose-response curve was split into. Default is 3.
#' @param fisher_pval_cutoff a numeric value that specifies the cutoff for the p-value obtained from
#' the fisher's exact test. Everything with a p-value below or equal to this threshold will be 
#' considered to pass.
#' @param complete_doses optional but recommended, a numeric vector that contains all theoretically 
#' observable doses for the given data. The doses need to be in ascending order. This argument is 
#' important for datasets that are small and might therefore not contain all the actual doses.
#' @param max_dose_replicates optional but recommended (required if `complete_doses` was provided), 
#' a numeric vector that contains the maximum number of observations for each dose point. Needs 
#' to be in the same order as `complete_doses`.
#' @param for_plot a logical value that determines the output of the function. Default is FALSE.
#' 
#' @return A data frame that contains the fisher's exact test p-value as well as information on 
#' the dataset split. The classification of MNAR cases is done based on the specified p-value
#' cutoff. If `for_plot` is `TRUE` a data frame with information required for a plot of each case 
#' is returned.
#' 
#' @import dplyr
#' @import tidyr
#' @importFrom tibble tibble 
#' @importFrom rlang .data enquo !! as_name
#' @importFrom purrr map_dfr
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' 
#' # Load libraries
#' library(dplyr)
#' library(ggplot2)
#' 
#' # Set seed
#' set.seed(123) # Makes example reproducible
#' 
#' # Create example data
#' data_test <- create_synthetic_data(
#'   n_proteins = 2,
#'   frac_change = 1,
#'   n_replicates = 3,
#'   n_conditions = 8,
#'   method = "dose_response",
#'   concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'   additional_metadata = FALSE
#' ) %>% 
#'   # remove random sample to simulate uneven sample numbers
#'   filter(!sample %in% c("sample_4", "sample_19")) 
#' 
#' # Run MNAR classification
#' MNAR_classification <- data_test %>% 
#'   find_drc_missingness(
#'     sample = sample,
#'     grouping = peptide,
#'     response = peptide_intensity_missing,
#'     dose = concentration,
#'     # set p-value cutoff low for highest confidence cases
#'     fisher_pval_cutoff = 0.0001, 
#'     complete_doses = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'     # some concentrations have less replicates due to removed samples
#'     max_dose_replicates = c(3, 2, 3, 3, 3, 3, 2, 3) 
#'   ) 
#' 
#' MNAR_cases <- MNAR_classification %>% 
#'   filter(MNAR) %>% 
#'   pull(peptide)
#' 
#' # Export cases for plotting
#' MNAR_for_plot <- data_test %>% 
#'   find_drc_missingness(
#'     sample = sample,
#'     grouping = peptide,
#'     response = peptide_intensity_missing,
#'     dose = concentration,
#'     fisher_pval_cutoff = 0.0001,
#'     complete_doses = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'     # some concentrations have less replicates due to removed samples
#'     max_dose_replicates = c(3, 2, 3, 3, 3, 3, 2, 3),
#'     for_plot = TRUE
#'   ) %>% 
#'   filter(peptide %in% MNAR_cases)
#' 
#' # Plot MNAR cases
#' ggplot(MNAR_for_plot %>% arrange(desc(type))) +
#'   geom_col(aes(dose_order, count, fill = type), position = "identity") +
#'   geom_line(MNAR_for_plot %>% distinct(peptide, cluster_split_id, score_total_error), mapping = aes(cluster_split_id, score_total_error / 10), linewidth = 1.5, alpha = 0.8, col = protti::protti_colours[1]) +
#'   facet_wrap(~ peptide) + 
#'   scale_fill_manual(values = c(protti::protti_colours[2], "grey"))+
#'   scale_y_continuous(sec.axis = sec_axis(~ .*10, name = "Score")) +
#'   labs(x = "Concentration Order", y = "Count", fill = "Count Type") + 
#'   theme_bw() +
#'   theme(plot.title = ggplot2::element_text(size = 18),
#'         axis.text.x = ggplot2::element_text(size = 10),
#'         axis.text.y = ggplot2::element_text(size = 10, color = protti::protti_colours[2]),
#'         axis.text.y.right = element_text(color = protti::protti_colours[1]),
#'         axis.title.y.right = element_text(color= protti::protti_colours[1]),
#'         axis.title.y = ggplot2::element_text(size = 12, color = protti::protti_colours[2]),
#'         axis.title.x = ggplot2::element_text(size = 12),
#'         legend.title = ggplot2::element_text(size = 12),
#'         legend.text = ggplot2::element_text(size = 12),
#'         strip.text.x = ggplot2::element_text(size = 10),
#'         strip.background = element_blank())
#'         
find_drc_missingness <- function(data,
                                 sample,
                                 grouping,
                                 response,
                                 dose,
                                 min_cluster_size = 3,
                                 fisher_pval_cutoff,
                                 complete_doses = NULL,
                                 max_dose_replicates = NULL,
                                 for_plot = FALSE){
  
  # Check how many groups there are in the dataset
  groups_n <- data %>% 
    dplyr::distinct({{grouping}}) %>% 
    nrow()
  
  # If the number of entities in cases is low and complete_doses or
  # complte_dose_replicates were not provided warn the user.
  if(groups_n <= 10 & (missing(complete_doses) | missing(max_dose_replicates))){
    warning(strwrap('You have less than 10 entities for the grouping variable and did not supply 
                    "complete_doses" or "max_dose_replicates". There is a chance that you do not represent 
                    all doses and replicates in the provided entities so it is strongly recommended 
                    to provide these two arguments to the function for small datasets.', prefix = "\n", initial = ""))
  }
  
  # Extract dose replicate information from data
  # Is only used if not provided by the user
  dose_replicate_data <- data %>% 
    tidyr::drop_na({{ response }}) %>% 
    dplyr::distinct({{ dose }}, {{ sample }}) %>% 
    dplyr::group_by({{ dose }}) %>% 
    dplyr::summarise(max_replicates = dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange({{ dose }})
  
  # Extract dose information from data if not provided
  if(missing(complete_doses)){
    if (!missing(max_dose_replicates)) {
      # It is not possible to only provide max_dose_replicates
      # without complete doses
      stop('Please also provide "complete_doses" if you provide "max_dose_replicates"!')
    }
    
    complete_doses <- dose_replicate_data %>% 
      dplyr::pull({{ dose }})
    
  } else {
    
    dose_replicate_data <- dose_replicate_data %>% 
      tidyr::complete({{ dose }} := complete_doses) %>% 
      dplyr::mutate(max_replicates = ifelse(is.na(.data$max_replicates), 0, .data$max_replicates))
    
  }
  
  # Extract replicate information from data if not provided
  if(!missing(max_dose_replicates)){
    dose_replicate_data <- tibble::tibble(max_replicates = max_dose_replicates,
                                      {{ dose }} := complete_doses)
  }

  # Add doses not detected for each grouping entity to the data
  # Then add the max replicates to the data
  data_complete <- data %>% 
    dplyr::distinct({{grouping}}, {{dose}}, {{sample}}, {{response}}) %>% 
    tidyr::drop_na({{ response }}) %>% 
    tidyr::complete({{ grouping }}, {{ dose }} := complete_doses) %>% 
    dplyr::left_join(dose_replicate_data, by = c(rlang::as_name(rlang::enquo(dose))))
    
  # Create all possible splits of the data and calculate the score
  # The best clustering is based on the minimal summed cluster variance
  cluster_splits <- 1:(nrow(dose_replicate_data)-1)
  
  cluster_prep <- cluster_splits %>% 
    purrr::map_dfr(.f = ~{
      data_complete %>% 
        dplyr::mutate(clusters = {{ dose }} %in% complete_doses[1:.x]) %>% 
        dplyr::mutate(cluster_split_id = .x)
    }) %>% 
    dplyr::group_by(.data$cluster_split_id, .data$clusters, {{dose}}, {{grouping}}) %>% 
    dplyr::mutate(observed_replicates = ifelse(is.na({{response}}), 0, dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct(.data$observed_replicates, .data$cluster_split_id, .data$clusters, {{dose}}, {{grouping}}, .data$max_replicates) %>% 
    dplyr::group_by(.data$cluster_split_id, .data$clusters, {{grouping}}) %>% 
    dplyr::mutate(
      mean = mean(.data$observed_replicates),
      max_mean = mean(.data$max_replicates),
      n_0 = sum(.data$observed_replicates == 0),
      detected = sum(.data$observed_replicates),
      undetected = sum(.data$max_replicates) - .data$detected,
      cluster_total_error = sum((.data$observed_replicates - mean(.data$observed_replicates))^2)) 
  
  # cluster_prep is also used later when data is used for plotting, because it still contains observed_replicates
  
  clusters <- cluster_prep %>% 
    dplyr::distinct( 
      {{grouping}},
      .data$cluster_split_id, 
      .data$clusters, 
      .data$mean,
      .data$max_mean,
      .data$cluster_total_error,
      .data$detected,
      .data$undetected,
      .data$n_0
    ) %>% 
    dplyr::ungroup(.data$clusters) %>% 
    dplyr::mutate(score_total_error = sum(.data$cluster_total_error)) %>% 
    dplyr::arrange(.data$score_total_error) %>% 
    dplyr::group_by({{ grouping }}) %>% 
    dplyr::mutate(min_score = min(.data$score_total_error)) %>% 
    dplyr::select(-"cluster_total_error") %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(clusters = ifelse(.data$clusters, "cluster1", "cluster2")) %>% 
    tidyr::pivot_wider(values_from = c(.data$n_0, .data$mean, .data$max_mean, .data$detected, .data$undetected), 
                       names_from = .data$clusters) 
  
  # Perform a fisher's exact test on the best cluster split for each 
  # entity in the grouping column
  test_data_prep <- clusters %>% 
    dplyr::filter(.data$cluster_split_id >= min_cluster_size & 
             .data$cluster_split_id <= (nrow(dose_replicate_data) - min_cluster_size)) %>%
    dplyr::filter(.data$min_score == .data$score_total_error) %>% 
    dplyr::select(-"min_score") 
  
  # Find all precursors missing from the data
  missing <- clusters %>% 
    dplyr::filter(!{{ grouping }} %in% (test_data_prep %>% dplyr::pull({{grouping}}))) %>% 
    dplyr::pull({{ grouping}}) %>% 
    unique()
  
  if (nrow(test_data_prep) != 0){
    test_data <- test_data_prep %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(pval = fisher.test(matrix(c(.data$detected_cluster1, .data$detected_cluster2,
                                                .data$undetected_cluster1, .data$undetected_cluster2), nrow = 2))$p.value) %>%
      dplyr::ungroup() %>% 
      bind_rows(tibble::tibble({{grouping}} := missing))
    
  } else {
    test_data <- tibble::tibble({{grouping}} := missing,
                                cluster_split_id = NA,
                                score_total_error = NA,
                                n_0_cluster1 = NA,
                                n_0_cluster2 = NA,
                                pval = NA,
                                mean_cluster1 = NA, 
                                mean_cluster2 = NA,
                                max_mean_cluster1 = NA,
                                max_mean_cluster2 = NA,
                                detected_cluster1 = NA,
                                detected_cluster2 = NA,
                                undetected_cluster1 = NA,
                                undetected_cluster2 = NA)
  }
  
  # Apply cutoffs to data
  result <- test_data %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(incomplete_mean = min(.data$mean_cluster1, .data$mean_cluster2),
           complete_mean = max(.data$mean_cluster1, .data$mean_cluster2),
           incomplete_completeness = .data$incomplete_mean / ifelse(.data$incomplete_mean == .data$mean_cluster1, 
                                                                    .data$max_mean_cluster1, 
                                                                    .data$max_mean_cluster2),
           complete_completeness = .data$complete_mean / ifelse(.data$complete_mean == .data$mean_cluster1, 
                                                           .data$max_mean_cluster1, 
                                                           .data$max_mean_cluster2),
           delta_completeness = .data$complete_completeness - .data$incomplete_completeness) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(incomplete_cluster_size = ifelse(.data$incomplete_mean == .data$mean_cluster1, 
                                          .data$cluster_split_id, 
                                          nrow(dose_replicate_data) - .data$cluster_split_id),
           incomplete_n_0 = ifelse(.data$incomplete_mean == .data$mean_cluster1, 
                                   .data$n_0_cluster1, 
                                   .data$n_0_cluster2)) %>% 
    # dplyr::mutate(min_n_cutoff = (.data$incomplete_cluster_size - 2 ) / .data$incomplete_cluster_size) %>%
    # dplyr::mutate(test = .data$incomplete_mean <= .data$min_n_cutoff) %>% 
    dplyr::mutate(MNAR = .data$complete_completeness > 0.5 & 
                    .data$incomplete_completeness < 0.25 &
                    .data$pval <= fisher_pval_cutoff) %>% 
    dplyr::mutate(MNAR = ifelse(is.na(MNAR), FALSE, MNAR)) %>% 
    dplyr::arrange(.data$pval)
  
  if(for_plot == FALSE){
    return(result)
  }

  # create data for plotting
  cluster_prep_subset <- cluster_prep %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(tibble::tibble({{ dose }} := complete_doses, dose_order = 1:length(complete_doses)), by = rlang::as_name(rlang::enquo(dose))) %>%
    dplyr::distinct({{ grouping }}, .data$cluster_split_id, .data$observed_replicates, .data$max_replicates, .data$dose_order)
  
  cluster_subset <- clusters %>% 
    dplyr::distinct({{ grouping }}, .data$cluster_split_id, .data$score_total_error)
  
  result_for_plot <- cluster_prep_subset %>% 
    dplyr::left_join(cluster_subset, by = c(rlang::as_name(rlang::enquo(grouping)), 
                                            "cluster_split_id")) %>% 
    tidyr::pivot_longer(cols = c(.data$observed_replicates, .data$max_replicates), values_to = "count", names_to = "type") %>% 
    dplyr::mutate(type = ifelse(type == "observed_replicates", "Observed", "Possible")) %>% 
    dplyr::mutate(cluster_split_id = .data$cluster_split_id + 0.5)
  
  result_for_plot
}