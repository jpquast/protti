#' Assignment of missingness types
#'
#' The type of missingness (missing at random, missing not at random) is assigned based on the comparison of a reference condition and every other condition.
#'
#' @param data a data frame containing at least the input variables.
#' @param sample The column in the data frame containing the sample name.
#' @param condition The column in the data frame containing the conditions.
#' @param grouping The column in the data frame containing precursor or peptide identifiers.
#' @param intensity The column in the data frame containing intensity values.
#' @param ref_condition a character vector providing the condition that is used as a reference for missingness determination.
#' Instead of providing one reference condition, "all" can be supplied, which will create all pairwise condition pairs. By default \code{ref_condition = "all"}.
#' @param completeness_MAR The minimal degree of data completeness to be considered as MAR. Value has to be between 0 and 1, default is 0.7.
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each
#' condition to be considered as MAR. This number is always at least 1.
#' @param completeness_MNAR The maximal degree of data completeness to be considered as MNAR. Value has to be between 0 and 1, default is 0.20.
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number is the maximal number of observations for one
#' condition to be considered as MNAR when the other condition is complete.
#' @param retain_columns A vector indicating if certain columns should be retained from the input data frame. Default is not retaining
#' additional columns \code{retain_columns = NULL}. Specific columns can be retained by providing their names (not in quotations marks,
#' just like other column names, but in a vector).
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
#' @importFrom tibble tibble as_tibble
#' @importFrom stringr str_extract
#' @importFrom rlang .data enquo !! as_name
#' @importFrom purrr map_df
#' @importFrom magrittr %>%
#' @importFrom utils combn
#' @export
#'
#' @examples
#' \dontrun{
#' assign_missingness(
#'   data,
#'   sample = r_file_name,
#'   condition = r_condition,
#'   grouping = eg_precursor_id,
#'   intensity = normalised_intensity_log2,
#'   retain_columns = c(pg_protein_accessions)
#' )
#' }
assign_missingness <- function(data, sample, condition, grouping, intensity, ref_condition = "all", completeness_MAR = 0.7, completeness_MNAR = 0.20, retain_columns = NULL) {
  . <- NULL
  if (!(ref_condition %in% unique(dplyr::pull(data, {{ condition }}))) & ref_condition != "all") {
    stop("The name provided to ref_condition cannot be found in your conditions! Please provide a valid reference condition.")
  }
  
  if(ref_condition == "all"){
    # creating all pairwise comparisons
    all_conditions <- unique(dplyr::pull(data, {{ condition }}))
    
    all_combinations <- tibble::as_tibble(t(utils::combn(all_conditions, m = 2))) %>% 
      dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2)) 
    
    message('"all" was provided as reference condition. All pairwise comparisons are created from the conditions and assigned their missingness.\n The created comparisons are: \n', paste(all_combinations$combinations, collapse =  "\n"))
  }

  if(ref_condition != "all"){
  conditions_no_ref <- unique(pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition]
  
  all_combinations <- tibble::tibble(V1 = conditions_no_ref, V2 = ref_condition) %>% 
    dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2)) 
  }
  
  # create dataframe that contains all combinations to be tested
  all_combinations <- all_combinations %>% 
    tidyr::pivot_longer(cols = c(.data$V1, .data$V2), names_to = "name", values_to = rlang::as_name(rlang::enquo(condition))) %>% 
    dplyr::select(-.data$name) %>% 
    dplyr::group_by({{ condition }}) %>% 
    dplyr::mutate(comparison = list(.data$combinations)) %>% 
    dplyr::distinct(.data$comparison, {{ condition }}) 

  data_prep <- data %>%
    dplyr::ungroup() %>% 
    dplyr::distinct({{ sample }}, {{ condition }}, {{ grouping }}, {{ intensity }}) %>%
    tidyr::complete(nesting(!!ensym(sample), !!ensym(condition)), !!ensym(grouping)) %>%
    dplyr::group_by({{ grouping }}, {{ condition }}) %>%
    dplyr::mutate(n_detect = sum(!is.na({{ intensity }}))) %>%
    dplyr::mutate(n_replicates = dplyr::n()) %>%
    dplyr::ungroup() %>% 
    dplyr::left_join(all_combinations, by = rlang::as_name(rlang::enquo(condition))) %>% 
    tidyr::unnest(.data$comparison)

  result <- data_prep %>%
    dplyr::mutate(type = ifelse({{ condition }} == stringr::str_extract(.data$comparison, pattern = "(?<=_vs_).+"), "control", "treated")) %>%
    split(.$comparison) %>%
    purrr::map_df(~ .x %>%
      tidyr::pivot_wider(names_from = .data$type, values_from = .data$n_detect) %>%
      dplyr::group_by({{ grouping }}) %>%
      tidyr::fill(treated, control, .direction = "updown") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(missingness = dplyr::case_when(
        .data$control == .data$n_replicates & .data$treated == .data$n_replicates ~ "complete",
        .data$control <= floor(n_replicates * completeness_MNAR) & .data$treated == .data$n_replicates ~ "MNAR",
        .data$control == .data$n_replicates & .data$treated <= floor(n_replicates * completeness_MNAR) ~ "MNAR",
        .data$control >= max(floor(.data$n_replicates * completeness_MAR), 1) & .data$treated >= max(floor(.data$n_replicates * completeness_MAR), 1) ~ "MAR"
      ))) %>%
    dplyr::select(-c(.data$control, .data$n_replicates, .data$treated)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange({{ grouping }})

  if (missing(retain_columns)) {
    return(result)
  } else {
    join_result <- data %>%
      dplyr::select(!!enquo(retain_columns), colnames(result)[!colnames(result) %in% c("comparison", "missingness")]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c("comparison", "missingness")]) %>%
      dplyr::arrange({{ grouping }})

    return(join_result)
  }
}
