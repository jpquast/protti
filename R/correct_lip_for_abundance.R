#' Protein abundance correction
#'
#' Performs the correction of LiP-peptides for changes in protein abundance and
#' calculates their significance using a t-test
#'
#' @param lip_data a data frame containing at least the input variables. Ideally,
#' the result from the \code{calculate_diff_abundance} function is used.
#' @param trp_data a data frame containing at least the input variables minus the grouping column. Ideally,
#' the result from the \code{calculate_diff_abundance} function is used.
#' @param protein_id a character column in the \code{lip_data} and \code{trp_data} data frames
#' that contains protein identifiers.
#' @param grouping a character column in the \code{lip_data} data frame that contains precursor or
#' peptide identifiers.
#' @param comparison a character column in the \code{lip_data} and \code{trp_data} data frames
#' that contains the comparisons between conditions.
#' @param diff a numeric column in the \code{lip_data} and \code{trp_data} data frames
#' that contains log2-fold changes for peptide or protein quantities.
#' @param n_obs a numeric column in the \code{lip_data} and \code{trp_data} data frames
#' containing the number of observations used to calculate fold changes.
#' @param std_error a numeric column in the \code{lip_data} and \code{trp_data} data frames
#' containing the standard error of fold changes.
#' @param p_adj_method a character value, specifies the p-value correction method. Possible
#' methods are c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default
#' method is \code{"BH"}.
#' @param retain_columns a vector indicating if certain columns should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector). Please note that if you retain columns that have multiple
#' rows per grouped variable there will be duplicated rows in the output.
#' @param method a character value, specifies the method used to estimate the degrees of freedom.
#' Possible methods are c("satterthwaite", "no_df_approximation"). \code{satterthwaite} uses the Welch-Satterthwaite
#' equation to estimate the pooled degrees of freedom, as described in https://doi.org/10.1016/j.mcpro.2022.100477 and
#' implemented in the MSstatsLiP package. This approach respects the number of protein measurements for the degrees of freedom.
#' \code{no_df_approximation} just takes the number of peptides into account when calculating the degrees of freedom.
#'
#' @return a data frame containing corrected differential abundances (\code{adj_diff}, adjusted
#' standard errors (\code{adj_std_error}), degrees of freedom (\code{df}), pvalues (\code{pval}) and
#' adjusted p-values (\code{adj_pval})
#'
#' @author Aaron Fehr
#' @import dplyr
#' @importFrom rlang .data enquo sym as_name expr := !!
#' @export
#'
#' @examples
#'
#' # Load libraries
#'
#' library(dplyr)
#'
#' # Load example data and simulate tryptic data by summing up precursors
#'
#' data <- rapamycin_10uM
#'
#' data_trp <- data %>%
#'   dplyr::group_by(pg_protein_accessions, r_file_name) %>%
#'   dplyr::mutate(pg_quantity = sum(fg_quantity)) %>%
#'   dplyr::distinct(
#'     r_condition,
#'     r_file_name,
#'     pg_protein_accessions,
#'     pg_quantity
#'   )
#'
#'
#' # Calculate differential abundances for LiP and Trp data
#'
#' diff_lip <- data %>%
#'   dplyr::mutate(fg_intensity_log2 = log2(fg_quantity)) %>%
#'   assign_missingness(
#'     sample = r_file_name,
#'     condition = r_condition,
#'     intensity = fg_intensity_log2,
#'     grouping = eg_precursor_id,
#'     ref_condition = "control",
#'     retain_columns = "pg_protein_accessions"
#'   ) %>%
#'   calculate_diff_abundance(
#'     sample = r_file_name,
#'     condition = r_condition,
#'     grouping = eg_precursor_id,
#'     intensity_log2 = fg_intensity_log2,
#'     comparison = comparison,
#'     method = "t-test",
#'     retain_columns = "pg_protein_accessions"
#'   )
#'
#'
#' diff_trp <- data_trp %>%
#'   dplyr::mutate(pg_intensity_log2 = log2(pg_quantity)) %>%
#'   assign_missingness(
#'     sample = r_file_name,
#'     condition = r_condition,
#'     intensity = pg_intensity_log2,
#'     grouping = pg_protein_accessions,
#'     ref_condition = "control"
#'   ) %>%
#'   calculate_diff_abundance(
#'     sample = r_file_name,
#'     condition = r_condition,
#'     grouping = pg_protein_accessions,
#'     intensity_log2 = pg_intensity_log2,
#'     comparison = comparison,
#'     method = "t-test"
#'   )
#'
#' # Correct for abundance changes
#'
#' corrected <- correct_lip_for_abundance(
#'   lip_data = diff_lip,
#'   trp_data = diff_trp,
#'   protein_id = pg_protein_accessions,
#'   grouping = eg_precursor_id,
#'   retain_columns = c("missingness"),
#'   method = "satterthwaite"
#' )
#'
#' head(corrected, n = 10)
correct_lip_for_abundance <- function(
    lip_data,
    trp_data,
    protein_id,
    grouping,
    comparison = comparison,
    diff = diff,
    n_obs = n_obs,
    std_error = std_error,
    p_adj_method = "BH",
    retain_columns = NULL,
    method = c("satterthwaite", "no_df_approximation")) {
  method <- match.arg(method)

  se_pep <- rlang::sym(paste0(rlang::as_name(rlang::enquo(std_error)), "_pep"))
  se_prot <- rlang::sym(paste0(rlang::as_name(rlang::enquo(std_error)), "_prot"))
  diff_pep <- rlang::sym(paste0(rlang::as_name(rlang::enquo(diff)), "_pep"))
  diff_prot <- rlang::sym(paste0(rlang::as_name(rlang::enquo(diff)), "_prot"))
  n_pep <- rlang::sym(paste0(rlang::as_name(rlang::enquo(n_obs)), "_pep"))
  n_prot <- rlang::sym(paste0(rlang::as_name(rlang::enquo(n_obs)), "_prot"))

  temp_lip_data <- lip_data %>%
    dplyr::select(!!enquo(retain_columns), {{ comparison }}, {{ protein_id }}, {{ grouping }}, {{ diff }}, {{ n_obs }}, {{ std_error }}) %>%
    dplyr::distinct()

  temp_trp_data <- trp_data %>%
    dplyr::distinct({{ comparison }}, {{ protein_id }}, {{ diff }}, {{ n_obs }}, {{ std_error }})


  test <- temp_lip_data %>%
    dplyr::distinct({{ comparison }}, {{ protein_id }}, {{ grouping }})

  if (nrow(test) != nrow(temp_lip_data)) {
    message("Warning: Your data frame contains dublicated values due to retained columns. This will affect the multiple testing correction.")
  }

  combined_data <- dplyr::left_join(
    x = temp_lip_data,
    y = temp_trp_data,
    by = c(rlang::as_name(rlang::enquo(comparison)), rlang::as_name(rlang::enquo(protein_id))),
    suffix = c("_pep", "_prot")
  )

  n_unmatched <- combined_data %>%
    dplyr::filter(is.na(!!diff_prot) == T) %>%
    nrow()

  percent_unmatched <- round(n_unmatched / nrow(combined_data) * 100, 2)

  if (n_unmatched != 0) {
    message(paste0("No protein data was available for ", n_unmatched, " peptides (", percent_unmatched, "% of dataset)."))
  }


  if (method == "satterthwaite") {
    corrected_data <- combined_data %>%
      dplyr::mutate(
        adj_diff = !!diff_pep - !!diff_prot,
        adj_std_error = sqrt((!!se_pep)**2 + (!!se_prot)**2),
        numer = ((!!se_pep)**2 + (!!se_prot)**2)**2,
        denom = ((!!se_pep)**4 / (!!n_pep - 2) + (!!se_prot)**4 / (!!n_prot - 2)),
        df = numer / denom,
        tval = adj_diff / adj_std_error,
        pval = 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
      )

    adjusted_data <- dplyr::left_join(
      x = corrected_data,
      y = corrected_data %>%
        dplyr::filter(is.na(pval) == FALSE) %>%
        dplyr::group_by({{ comparison }}) %>%
        dplyr::mutate(adj_pval = p.adjust(pval, method = {{ p_adj_method }})),
      by = colnames(corrected_data)
    ) %>%
      dplyr::select(-numer, -denom, -tval)

    return(adjusted_data)
  }

  if (method == "no_df_approximation") {
    corrected_data <- combined_data %>%
      dplyr::mutate(
        adj_diff = !!diff_pep - !!diff_prot,
        adj_std_error = sqrt((!!se_pep)**2 + (!!se_prot)**2),
        df = !!n_pep - 2,
        tval = adj_diff / adj_std_error,
        pval = 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
      )

    adjusted_data <- dplyr::left_join(
      x = corrected_data,
      y = corrected_data %>%
        dplyr::filter(is.na(pval) == FALSE) %>%
        dplyr::mutate(adj_pval = p.adjust(pval, method = {{ p_adj_method }})),
      by = colnames(corrected_data)
    ) %>%
      dplyr::select(-tval)


    return(adjusted_data)
  }
}
