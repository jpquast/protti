#' Assignment of missingness types
#'
#' The type of missingness (missing at random, missing not at random) is assigned based on the comparison of a reference condition and every other condition.
#'
#' @param data a data frame containing at least the input variables.
#' @param sample the name of the column containing the sample names.
#' @param condition the name of the column containing the conditions.
#' @param grouping the name of the column containing precursor or peptide identifiers.
#' @param intensity the name of the column containing intensity values.
#' @param noise optional, the name of the column containing noise values for imputation. Only needed if noise imputation is performed.
#' @param ref_condition the condition that is used as a reference for missingness determination. By default \code{ref_condition = "control"}.
#' @param completeness_MAR the minimal degree of data completeness to be considered as MAR. Value has to be between 0 and 1, default is 0.7. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each condition to be considered as MAR. This number is always at least 1. 
#' @param completeness_MNAR the maximal degree of data completeness to be considered as MNAR. Value has to be between 0 and 1, default is 0.25. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the maximal number of observations for one condition to be considered as MNAR when the other condition is complete.
#'
#' @return A data frame that contains the reference condition paired with each treatment condition. The \code{comparison} column contains the comparison name for the specific reference/treatment pair. The \code{missingness} column reports the type of missingness. 
#' \itemize{
#' \item{"complete": }{No missing values for every replicate of this reference/treatment pair for the specific grouping variable.}
#' \item{"MNAR": }{Missing not at random. All replicates of either the reference or treatment condition have missing values for the specific grouping variable.}
#' \item{"MAR": }{Missing at random. At least n-1 replicates have missing values for the reference/treatment pair for the specific grouping varible.}
#' } 
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
#' @importFrom rlang as_name
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
  data <-  data %>%
    dplyr::distinct({{sample}}, {{condition}}, {{grouping}}, {{intensity}}, {{noise}}) %>%
    tidyr::complete(nesting(!!ensym(sample), !!ensym(condition)), !!ensym(grouping)) %>%
    dplyr::mutate(detect = !is.na({{intensity}})) %>%
    dplyr::group_by({{grouping}}, {{condition}}) %>%
    dplyr::mutate(n_detect = sum(.data$detect)) %>%
    dplyr::mutate(n_replicates = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comparison = ifelse(
      {{condition}} %in% unique(pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition],
      paste0(ref_condition, "_vs_", {{condition}}),
      list(paste0(ref_condition, "_vs_", unique(pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition]))
    ))%>%
    tidyr::unnest(.data$comparison) 
  
  data %>%
    dplyr::distinct({{condition}}, {{grouping}}, .data$n_detect, .data$n_replicates, .data$comparison) %>%
    dplyr::group_by({{grouping}}, .data$comparison) %>%
    dplyr::mutate(missingness = dplyr::case_when(.data$n_detect == .data$n_replicates & lag(.data$n_detect) == .data$n_replicates ~ "complete",
                                                 .data$n_detect <= floor(n_replicates * completeness_MNAR) & lag(.data$n_detect) == .data$n_replicates ~ "MNAR",
                                                 .data$n_detect == .data$n_replicates & lag(.data$n_detect) <= floor(n_replicates * completeness_MNAR) ~ "MNAR",
                                                 .data$n_detect >= max(floor(.data$n_replicates * completeness_MAR), 1) & lag(.data$n_detect) >= max(floor(.data$n_replicates * completeness_MAR), 1) ~ "MAR")) %>%
    tidyr::fill(.data$missingness, .direction = "updown") %>%
    dplyr::select(-.data$n_detect, -.data$n_replicates) %>%
    dplyr::left_join(data, by = c("comparison", rlang::as_name(enquo(condition)), rlang::as_name(enquo(grouping)))) %>%
    dplyr::select(-.data$n_detect, -.data$n_replicates, -.data$detect) %>%
    dplyr::ungroup()
}