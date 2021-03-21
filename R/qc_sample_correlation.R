#' Correlation based hirachical clustering of samples
#'
#' A correlation heatmap is created that uses hirachical clustering to determine sample similarity.
#'
#' @param data a dataframe contains at least the input variables.
#' @param sample the name of the column containing the sample names.
#' @param grouping the name of the column containing precursor or peptide identifiers.
#' @param intensity_log2 the name of the column containing log2 intensity values.
#' @param condition the name of the column containing the conditions.
#' @param digestion optional, the name of the column containing information about the digestion method used. Eg. "LiP" or "tryptic control".
#' @param run_order optional, the name of the column containing the order in which samples were measured. Useful to investigate batch effects due to run order.
#' @param method the method to be used for correlation. \code{"spearman"} is the default but can be changed to \code{"pearson"} or \code{"kendall"}.
#' @param interactive logical, default is \code{FALSE}. Determines if an interactive or static heatmap should be created using \code{heatmaply} or \code{pheatmap}, respectively.
#'
#' @return A correlation heatmap that compares each sample. The dendrogram is sorted by optimal leaf ordering.
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data as_name
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_sample_correlation(
#'   data,
#'   sample = r_file_name,
#'   grouping = eg_precursor_id,
#'   intensity_log2 = intensity_log2,
#'   condition = r_condition
#' )
#' }
qc_sample_correlation <- function(data, sample, grouping, intensity_log2, condition, digestion = NULL, run_order = NULL, method = "spearman", interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data

  correlation <- data %>%
    dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity_log2 }}) %>%
    tidyr::pivot_wider(names_from = {{ sample }}, values_from = {{ intensity_log2 }}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(grouping))) %>%
    stats::cor(method = {{ method }}, use = "complete.obs")

  annotation <- data %>%
    dplyr::mutate({{ condition }} := as.character({{ condition }})) %>%
    dplyr::distinct({{ sample }}, {{ condition }}, {{ digestion }}, {{ run_order }}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(sample)))

  # Create list for colouring of annotations. Looks rather complicated but is simple.
  # Just takes into account that arguments might not be supplied by the user.
  n_conditions <- 0
  n_digest <- 0
  n_run_ord <- 0
  conditions_colour <- c()
  digest_colours <- c()
  run_ord_colours <- c()
  if (!missing(condition)) {
    conditions <- unique(dplyr::pull(annotation, {{ condition }}))
    n_conditions <- length(conditions)
    conditions_colours <- protti_colours[1:n_conditions]
    names(conditions_colours) <- conditions
  }
  if (!missing(digestion)) {
    digest <- unique(dplyr::pull(annotation, {{ digestion }}))
    n_digest <- length(digest)
    digest_colours <- protti_colours[(n_conditions + 1):(n_digest + n_conditions)]
    names(digest_colours) <- digest
  }
  if (!missing(run_order)) {
    run_ord <- unique(dplyr::pull(annotation, {{ run_order }}))
    n_run_ord <- length(run_ord)
    run_ord_colours <- protti_colours[(n_digest + n_conditions + 1):(n_run_ord + n_digest + n_conditions)]
    names(run_ord_colours) <- run_ord
  }
  annotation_colours <- list(conditions_colours, digest_colours, run_ord_colours)
  names(annotation_colours) <- c(if (!missing(condition)) {
    rlang::as_name(enquo(condition))
  } else {
    "condition"
  }, if (!missing(digestion)) {
    rlang::as_name(enquo(digestion))
  } else {
    "digestion"
  }, if (!missing(run_order)) {
    rlang::as_name(enquo(run_order))
  } else {
    "run_order"
  })

  if (interactive == TRUE) {
    if (!requireNamespace("heatmaply", quietly = TRUE)) {
      stop("Package \"heatmaply\" is needed for this function to work. Please install it.", call. = FALSE)
    }
    heatmap_interactive <-
      heatmaply::heatmaply(
        correlation,
        main = "Correlation based hirachical clustering of samples",
        col_side_colors = annotation,
        col_side_palette = c(annotation_colours[[1]], annotation_colours[[2]], annotation_colours[[3]]),
        k_col = NA,
        k_row = NA,
        plot_method = "plotly"
      )
    return(heatmap_interactive)
  }
  if (interactive == FALSE) {
    dependency_test <- c(
      dendextend = !requireNamespace("dendextend", quietly = TRUE),
      pheatmap = !requireNamespace("pheatmap", quietly = TRUE),
      viridis = !requireNamespace("viridis", quietly = TRUE),
      seriation = !requireNamespace("seriation", quietly = TRUE)
    )
    if (any(dependency_test)) {
      dependency_name <- names(dependency_test[dependency_test == TRUE])
      if (length(dependency_name) == 1) {
        stop("Package \"", paste(dependency_name), "\" is needed for this function to work. Please install it.", call. = FALSE)
      } else {
        stop("Packages \"", paste(dependency_name, collapse = "\" and \""), "\" are needed for this function to work. Please install them.", call. = FALSE)
      }
    }

    # Create heatmaply dendrogram for pheatmap
    distance <- stats::dist(correlation)
    hierachical_clustering <- stats::hclust(distance)
    dendrogram <- stats::as.dendrogram(hierachical_clustering)
    dendrogram_row <- dendextend::seriate_dendrogram(dendrogram, distance, method = "OLO")
    dendrogram_column <- dendextend::rotate(dendrogram_row, order = rev(labels(distance)[seriation::get_order(stats::as.hclust(dendrogram_row))]))

    heatmap_static <-
      pheatmap::pheatmap(
        correlation,
        cluster_rows = stats::as.hclust(dendrogram_row),
        cluster_cols = stats::as.hclust(dendrogram_column),
        annotation = annotation,
        annotation_colors = annotation_colours,
        main = "Correlation based hirachical clustering of samples",
        color = viridis::viridis(100)
      )
    return(heatmap_static)
  }
}
