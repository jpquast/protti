#' Assignment of missingness types
#'
#' The type of missingness (missing at random, missing not at random) is assigned based on the comparison of a reference condition and every other condition.
#'
#' @param data A data frame containing at least the input variables.
#' @param sample The column in the data data frame containing the sample name.
#' @param condition The column in the data data frame containing the conditions.
#' @param grouping The column in the data data frame containing precursor or peptide identifiers.
#' @param intensity The column in the data data frame containing intensity values.
#' @param noise optional, the column in the data data frame containing noise values for imputation. Only needed if noise imputation is performed.
#' @param ref_condition The condition that is used as a reference for missingness determination. By default \code{ref_condition = "control"}.
#' @param completeness_MAR The minimal degree of data completeness to be considered as MAR. Value has to be between 0 and 1, default is 0.7. 
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each 
#' condition to be considered as MAR. This number is always at least 1. 
#' @param completeness_MNAR The maximal degree of data completeness to be considered as MNAR. Value has to be between 0 and 1, default is 0.25. 
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number is the maximal number of observations for one 
#' condition to be considered as MNAR when the other condition is complete.
#'
#' @return A data frame that contains the reference condition paired with each treatment condition. The \code{comparison} column contains the comparison 
#' name for the specific treatment/reference pair. The \code{missingness} column reports the type of missingness. 
#' \itemize{
#' \item{"complete": }{No missing values for every replicate of this reference/treatment pair for the specific grouping variable.}
#' \item{"MNAR": }{Missing not at random. All replicates of either the reference or treatment condition have missing values for the specific grouping variable.}
#' \item{"MAR": }{Missing at random. At least n-1 replicates have missing values for the reference/treatment pair for the specific grouping varible.}
#' } 
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
#' @importFrom purrr map_df
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' assign_missingness(
#' data,
#' sample = r_file_name,
#' condition = r_condition,
#' grouping = eg_precursor_id,
#' intensity = normalised_intensity_log2
#' )
#' }
assign_missingness <- function(data, sample, condition, grouping, intensity, noise = NULL, ref_condition = "control", completeness_MAR = 0.7, completeness_MNAR = 0.25){
  . = NULL
  conditions_no_ref <- unique(pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition]
  
  data <-  data %>%
    dplyr::distinct({{sample}}, {{condition}}, {{grouping}}, {{intensity}}, {{noise}}) %>%
    tidyr::complete(nesting(!!ensym(sample), !!ensym(condition)), !!ensym(grouping)) %>%
    dplyr::group_by({{grouping}}, {{condition}}) %>%
    dplyr::mutate(n_detect = sum(!is.na({{intensity}}))) %>%
    dplyr::mutate(n_replicates = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comparison = ifelse(
      {{condition}} %in% conditions_no_ref,
      paste0({{condition}}, "_vs_", ref_condition),
      list(paste0(conditions_no_ref, "_vs_",  ref_condition))
    ))%>%
    tidyr::unnest(.data$comparison) 
  
  data %>%
    dplyr::mutate(type = ifelse({{condition}} == ref_condition, "control", "treated")) %>%
    split(.$comparison) %>%
    purrr::map_df(~ .x %>%
                    tidyr::pivot_wider(names_from = .data$type, values_from = .data$n_detect) %>%
                    dplyr::group_by({{grouping}}) %>%
                    tidyr::fill(treated, control, .direction = "updown") %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(missingness = dplyr::case_when(.data$control == .data$n_replicates & .data$treated == .data$n_replicates ~ "complete",
                                                                 .data$control <= floor(n_replicates * completeness_MNAR) & .data$treated == .data$n_replicates ~ "MNAR",
                                                                 .data$control == .data$n_replicates & .data$treated <= floor(n_replicates * completeness_MNAR) ~ "MNAR",
                                                                 .data$control >= max(floor(.data$n_replicates * completeness_MAR), 1) & .data$treated >= max(floor(.data$n_replicates * completeness_MAR), 1) ~ "MAR"))
    ) %>%
    dplyr::select(-c(.data$control, .data$n_replicates, .data$treated)) %>%
    dplyr::ungroup()
}