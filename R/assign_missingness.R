#' Assignment of missingness types
#'
#' The type of missingness (missing at random, missing not at random) is assigned based on the
#' comparison of a reference condition and every other condition.
#'
#' @param data a data frame containing at least the input variables.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' @param condition a character or numeric column in the \code{data} data frame that contains the
#' conditions.
#' @param grouping a character column in the \code{data} data frame that contains protein, precursor or
#' peptide identifiers.
#' @param intensity a numeric column in the \code{data} data frame that contains intensity values that
#' relate to the \code{grouping} variable.
#' @param ref_condition a character vector providing the condition that is used as a reference for
#' missingness determination. Instead of providing one reference condition, "all" can be supplied,
#' which will create all pairwise condition pairs. By default \code{ref_condition = "all"}.
#' @param completeness_MAR a numeric value that specifies the minimal degree of data completeness to
#' be considered as MAR. Value has to be between 0 and 1, default is 0.7. It is multiplied with
#' the number of replicates and then adjusted downward. The resulting number is the minimal number
#' of observations for each condition to be considered as MAR. This number is always at least 1.
#' @param completeness_MNAR a numeric value that specifies the maximal degree of data completeness to
#' be considered as MNAR. Value has to be between 0 and 1, default is 0.20. It is multiplied with
#' the number of replicates and then adjusted downward. The resulting number is the maximal number
#' of observations for one condition to be considered as MNAR when the other condition is complete.
#' @param retain_columns a vector that indicates columns that should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return A data frame that contains the reference condition paired with each treatment condition.
#' The \code{comparison} column contains the comparison name for the specific treatment/reference
#' pair. The \code{missingness} column reports the type of missingness.
#' * "complete": No missing values for every replicate of this reference/treatment pair for
#' the specific grouping variable.
#' * "MNAR": Missing not at random. All replicates of either the reference or treatment
#' condition have missing values for the specific grouping variable.
#' * "MAR": Missing at random. At least n-1 replicates have missing values for the
#' reference/treatment pair for the specific grouping varible.
#' * NA: The comparison is not complete enough to fall into any other category. It will not
#' be imputed if imputation is performed. For statistical significance testing these comparisons
#' are filtered out after the test and prior to p-value adjustment. This can be prevented by setting
#' `filter_NA_missingness = FALSE` in the `calculate_diff_abundance()` function.
#'
#' The type of missingness has an influence on the way values are imputeted if imputation is
#' performed subsequently using the `impute()` function. How each type of missingness is
#' specifically imputed can be found in the function description. The type of missingness
#' assigned to a comparison does not have any influence on the statistical test in the
#' `calculate_diff_abundance()` function.
#' @import dplyr
#' @import tidyr
#' @importFrom tibble tibble as_tibble
#' @importFrom stringr str_extract str_sort
#' @importFrom rlang .data enquo !! as_name
#' @importFrom purrr map_df
#' @importFrom magrittr %>%
#' @importFrom utils combn capture.output
#' @export
#'
#' @examples
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 10,
#'   frac_change = 0.5,
#'   n_replicates = 4,
#'   n_conditions = 2,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' head(data, n = 24)
#'
#' # Assign missingness information
#' data_missing <- assign_missingness(
#'   data,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   ref_condition = "all",
#'   retain_columns = c(protein)
#' )
#'
#' head(data_missing, n = 24)
assign_missingness <- function(data,
                               sample,
                               condition,
                               grouping,
                               intensity,
                               ref_condition = "all",
                               completeness_MAR = 0.7,
                               completeness_MNAR = 0.20,
                               retain_columns = NULL) {
  . <- NULL
  if (!(ref_condition %in% unique(dplyr::pull(data, {{ condition }}))) & ref_condition != "all") {
    stop(strwrap("The name provided to ref_condition cannot be found in your conditions!
Please provide a valid reference condition.", prefix = "\n", initial = ""))
  }

  if (ref_condition == "all") {
    # creating all pairwise comparisons
    all_conditions <- unique(dplyr::pull(data, {{ condition }}))

    all_combinations <- tibble::as_tibble(t(utils::combn(all_conditions, m = 2)),
      .name_repair = ~ make.names(., unique = TRUE)
    ) %>%
      dplyr::rename(
        V1 = .data$X,
        V2 = .data$X.1
      ) %>%
      dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2))

    message(
      strwrap('"all" was provided as reference condition. All pairwise comparisons are created
from the conditions and assigned their missingness. The created comparisons are:', prefix = "\n", initial = ""), "\n",
      paste(all_combinations$combinations, collapse = "\n")
    )
  }

  if (ref_condition != "all") {
    conditions_no_ref <- unique(pull(data, !!ensym(condition)))[!unique(pull(data, !!ensym(condition))) %in% ref_condition]

    all_combinations <- tibble::tibble(V1 = conditions_no_ref, V2 = ref_condition) %>%
      dplyr::mutate(combinations = paste0(.data$V1, "_vs_", .data$V2))
  }

  # create dataframe that contains all combinations to be tested
  all_combinations <- all_combinations %>%
    tidyr::pivot_longer(cols = c("V1", "V2"), names_to = "name", values_to = rlang::as_name(rlang::enquo(condition))) %>%
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
    tidyr::unnest("comparison")

  # check if there are any unequal replicate comparisons
  unequal_replicates <- data_prep %>%
    dplyr::arrange({{ condition }}) %>%
    dplyr::distinct(.data$n_replicates, .data$comparison) %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::filter(.data$n > 1) %>%
    dplyr::mutate(n_replicates = paste0(.data$n_replicates, collapse = "/"))

  if (any(unequal_replicates$n > 2)) {
    stop(
      "\n",
      strwrap('Some created comparisons seem to have more than two unequal number of replicates.
              This usually only happens if the wrong grouping variable was selected. Please check this!
              The grouping variable should split the dataset so that each sample of each condition only
              appears once for each element of the grouping. E.g. grouping peptide: Each peptide should
              only have sample_1 associated once with condition_1 and not twice or more often. If in this
              case grouping "protein" was inadvertently selected a protein might have multiple peptides, each
              containing sample_1 of condition_1, which means it appears more than once (appears as many times
              as there are peptides per protein). This means each condition can have an unequal number of replicates
              that is as high as the max number of proteins, which is not the correct calculation for replicates.',
        prefix = "\n", initial = ""
      ), "\n"
    )
  }

  unequal_replicates <- unequal_replicates %>%
    dplyr::distinct(.data$n_replicates, .data$comparison)

  if (nrow(unequal_replicates) != 0) {
    message("\n")
    message(
      strwrap("The following comparisons have been detected to have unequal replicate numbers.
              If this is intended please ignore this message. This function can appropriately deal
              with unequal replicate numbers.", prefix = "\n", initial = ""), "\n"
    )
    message(paste0(utils::capture.output(unequal_replicates), collapse = "\n"))
  }

  result <- data_prep %>%
    dplyr::mutate(type = ifelse({{ condition }} == stringr::str_extract(.data$comparison, pattern = "(?<=_vs_).+"),
      "control",
      "treated"
    )) %>%
    split(.$comparison) %>%
    purrr::map_df(.f = ~ .x %>%
      tidyr::pivot_wider(names_from = "type", values_from = c("n_detect", "n_replicates")) %>%
      dplyr::group_by({{ grouping }}) %>%
      tidyr::fill("n_detect_treated", "n_detect_control", "n_replicates_treated", "n_replicates_control", .direction = "updown") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(missingness = dplyr::case_when(
        .data$n_detect_control == .data$n_replicates_control &
          .data$n_detect_treated == .data$n_replicates_treated ~ "complete",
        .data$n_detect_control <= floor(n_replicates_control * completeness_MNAR) &
          .data$n_detect_treated == .data$n_replicates_treated ~ "MNAR",
        .data$n_detect_control == .data$n_replicates_control &
          .data$n_detect_treated <= floor(n_replicates_treated * completeness_MNAR) ~ "MNAR",
        .data$n_detect_control >= max(floor(.data$n_replicates_control * completeness_MAR), 1) &
          .data$n_detect_treated >= max(floor(.data$n_replicates_control * completeness_MAR), 1) ~ "MAR"
      ))) %>%
    dplyr::select(-c("n_detect_control", "n_detect_treated", "n_replicates_control", "n_replicates_treated")) %>%
    # Arrange by grouping but in a numeric order of the character vector.
    dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))

  if (missing(retain_columns)) {
    return(result)
  } else {
    join_result <- data %>%
      dplyr::ungroup() %>%
      dplyr::select(!!enquo(retain_columns), colnames(result)[!colnames(result) %in% c("comparison", "missingness")]) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c("comparison", "missingness")]) %>%
      # Arrange by grouping but in a numeric order of the character vector.
      dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE)))) %>%
      # propagation of consistent values to NA places
      dplyr::group_by({{ grouping }}) %>%
      dplyr::mutate(dplyr::across(!!enquo(retain_columns), ~ {
        # Check if all non-NA values are the same
        if (any(is.na(.x)) & dplyr::n_distinct(na.omit(.x)) == 1 & !any(is.na(.x) & !is.na({{ intensity }}))) {
          # Replace NA with the consistent value
          tidyr::replace_na(.x, unique(na.omit(.x)))
        } else {
          # Leave as is
          .x
        }
      })) %>%
      dplyr::ungroup()

    # Annotate sample related retained columns. These have a unique value for every sample.
    # Above we annotated any columns that had a consistent for every group, here the inconsistent ones are annotated

    sample_annotations <- join_result %>%
      dplyr::select(!!enquo(retain_columns), {{ intensity }}, {{ sample }}) %>%
      dplyr::select(
        dplyr::where(~ !any(is.na(.x) & !is.na(dplyr::pull(join_result, {{ intensity }}))) & any(is.na(.x))),
        {{ sample }},
        -{{ intensity }}
      ) %>%
      tidyr::drop_na() %>%
      dplyr::distinct() %>%
      dplyr::group_by({{ sample }}) %>%
      # drop the columns that contain multiple values per group
      # grouping doesn't work with selection so first we need to find the columns with the non-distinct values with the summary bellow
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ if (dplyr::n_distinct(.x) == 1) dplyr::first(.x) else NA), .groups = "drop") %>%
      dplyr::select(-dplyr::where(~ any(is.na(.x)))) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    join_result <- join_result %>%
      dplyr::rows_update(sample_annotations, by = rlang::as_name(rlang::enquo(sample)))

    return(join_result)
  }
}
