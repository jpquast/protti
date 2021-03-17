#' Label-free protein quantification
#'
#' Determines relative protein abundances from ion quantification. Only proteins with at least 3 peptides are considered for quantification.
#'
#' @param data A data frame that contains at least the input variables.
#' @param sample The name of the column containing the sample name.
#' @param protein_id The name of the column containing the protein accession numbers.
#' @param precursor The name of the column containing precursors. 
#' @param peptide The name of the column containing peptide sequences. This column is needed to filter for proteins with at least 3 
#' unique peptides. This can equate to more than three precursors. The quantification is done on the precursor level. 
#' @param intensity_log2 The name of the column containing log2 transformed precursor intensities. 
#' @param method A character vector specifying with which method protein quantities should be calculated. Possible options include \code{"sum"},
#' which takes the sum of all precursor intensities as the protein abundance. Another option is \code{"iq"}, which performs protein 
#' quantification based on a maximal peptide ratio extraction algorithm that is adapted from the MaxLFQ algorithm of the MaxQuant software. 
#' Functions from the \href{https://academic.oup.com/bioinformatics/article/36/8/2611/5697917}{\code{iq}} package are used. Default is \code{"iq"}. 
#' @param for_plot A logical indicating whether the result should be only protein intensities or protein intensities together with precursor 
#' intensities that can be used for plotting using \code{qc_protein_abundance}. Default is \code{FALSE}.
#' @param retain_columns A vector indicating if certain columns should be retained from the input data frame. Default is not retaining 
#' additional columns \code{retain_columns = NULL}. Specific columns can be retained by providing their names (not in quotations marks, 
#' just like other column names, but in a vector).
#' 
#' @return If \code{for_plot = FALSE}, protein abundances are returned, if \code{for_plot = TRUE} also precursor intensities are returned. The 
#' later output is ideal for plotting with \code{qc_protein_abundance} and can be filtered to only include protein abundances.
#'
#' @import dplyr
#' @import progress
#' @importFrom tidyr complete pivot_wider drop_na
#' @importFrom rlang .data := !! ensym as_name enquo
#' @importFrom tibble column_to_rownames as_tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom purrr map map2_df discard pluck
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_protein_abundance(
#' data,
#' sample = r_file_name,
#' protein_id = pg_protein_accessions,
#' precursor = eg_precursor_id,
#' peptide = pep_stripped_sequence,
#' intensity_log2 = normalised_intensity_log2,
#' method = "iq",
#' retain_columns = c(pg_protein_accessions)
#' )
#' }
calculate_protein_abundance <- function(data, sample, protein_id, precursor, peptide, intensity_log2, method = "iq", for_plot = FALSE, retain_columns = NULL){
  . = NULL
  
  # Filter out any proteins with less than 3 peptides
  input <- data %>%
    dplyr::distinct({{sample}}, {{protein_id}}, {{precursor}}, {{peptide}}, {{intensity_log2}}) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by({{protein_id}}) %>%
    dplyr::mutate(n_peptides = dplyr::n_distinct(!!rlang::ensym(peptide))) %>%
    dplyr::filter(.data$n_peptides > 2) %>% 
    dplyr::ungroup()
  
  if(method == "sum"){
    result <- input %>%
      dplyr::group_by({{sample}}, {{protein_id}}) %>% 
      dplyr::summarise({{intensity_log2}} := log2(sum(2^{{intensity_log2}})), .groups = "drop") 
    
    if(missing(retain_columns) & for_plot == FALSE) return(result)

    combined <- result %>% 
      dplyr::mutate({{precursor}} := "protein_intensity") %>% 
      dplyr::bind_rows(input)
  }
  if(method == "iq"){
    if (!requireNamespace("iq", quietly = TRUE)) {
      stop("Package \"iq\" is needed for this function to work. Please install it.", call. = FALSE)
    }
  pb <- progress::progress_bar$new(total = length(unique(dplyr::pull(input, {{protein_id}}))), format = "  Preparing data [:bar] :current/:total (:percent) :eta")
  
  input <- input %>%
    dplyr::distinct({{sample}}, {{protein_id}}, {{precursor}}, {{intensity_log2}}) %>%
    tidyr::complete(!!rlang::ensym(sample), nesting(!!rlang::ensym(precursor), !!rlang::ensym(protein_id))) %>%
    split(dplyr::pull(., {{protein_id}})) %>%
    purrr::map(.f = ~{pb$tick(); .x %>%
        dplyr::select(-{{protein_id}}) %>%
        tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity_log2}}) %>%
        tibble::column_to_rownames(rlang::as_name(rlang::enquo(precursor))) %>%
        as.matrix()}
    ) 

  pb <- progress::progress_bar$new(total = length(input), format = "  Applying maximal peptide ratio extraction algorithm [:bar] :current/:total (:percent) :eta")
  
  combined <- input %>% 
    purrr::map2_df(.y = names(.),
            .f = ~ {pb$tick(); iq::maxLFQ(.x) %>% 
                purrr::pluck("estimate") %>%
                matrix(ncol = ncol(.x), 
                       nrow = 1,
                       dimnames = list("protein_intensity", colnames(.x))) %>% 
                rbind(.x) %>%
                tibble::as_tibble(rownames = NA) %>% 
                tibble::rownames_to_column(var = rlang::as_name(rlang::enquo(precursor))) %>% 
                tidyr::pivot_longer(-{{precursor}}, names_to = rlang::as_name(rlang::enquo(sample)), values_to = rlang::as_name(rlang::enquo(intensity_log2))) %>% 
                dplyr::mutate({{protein_id}} := .y)}
    ) %>% 
    tidyr::drop_na()
  
  if(missing(retain_columns) & for_plot == TRUE) return(combined)
  
  result <- combined %>% 
    dplyr::filter({{precursor}} == "protein_intensity") %>% 
    dplyr::select(-{{precursor}})
  }
  
  if (!missing(retain_columns) & for_plot == FALSE) {
    result <- data %>% 
      dplyr::select(!!enquo(retain_columns), colnames(result)[!colnames(result) %in% c(rlang::as_name(rlang::enquo(intensity_log2)))]) %>% 
      dplyr::distinct() %>% 
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c(rlang::as_name(rlang::enquo(intensity_log2)))])
    
    return(result)
  } 
  if (!missing(retain_columns) & for_plot == TRUE) {
    combined <- data %>% 
      dplyr::select(!!enquo(retain_columns), colnames(combined)[!colnames(combined) %in% c(rlang::as_name(rlang::enquo(intensity_log2)), rlang::as_name(rlang::enquo(precursor)))]) %>% 
      dplyr::distinct() %>% 
      dplyr::right_join(combined, by = colnames(combined)[!colnames(combined) %in% c(rlang::as_name(rlang::enquo(intensity_log2)), rlang::as_name(rlang::enquo(precursor)))])
    
    return(combined)
  } 
  
  return(result)
}