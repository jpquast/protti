#' Correlation based hirachical clustering of samples
#'
#' A correlation heatmap is created that uses hirachical clustering to determine sample similarity.
#'
#' @param data a dataframe contains at least the input variables.
#' @param sample the name of the column containing the sample names.
#' @param grouping the name of the column containing precursor or peptide identifiers.
#' @param intensity the name of the column containing intensity values. Note: The input intensities should be log2 transformed.
#' @param condition the name of the column containing the conditions.
#' @param digestion optional, the name of the column containing information about the digestion method used. Eg. "LiP" or "tryptic control".
#' @param run_order optional, the name of the column containing the order in which samples were measured. Useful to investigate batch effects due to run order.
#' @param method the method to be used for correlation. \code{"spearman"} is the default but can be changed to \code{"pearson"} or \code{"kendall"}.
#' @param interactive logical, default is \code{TRUE}. Determines if an interactive or static heatmap should be created using \code{heatmaply} or \code{pheatmap}, respectively.
#'
#' @return A correlation heatmap that compares each sample. The dendrogram is sorted by optimal leaf ordering.
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap
#' @importFrom heatmaply heatmaply
#' @importFrom viridis viridis
#' @importFrom seriation get_order
#' @importFrom dendextend seriate_dendrogram
#' @importFrom dendextend rotate
#' @importFrom rlang .data
#' @importFrom rlang as_name
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' qc_sample_correlation(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' intensity = intensity_log2,
#' condition = r_condition
#' )
#' }
qc_sample_correlation <- function(data, sample, grouping, intensity, condition, digestion = NULL, run_order = NULL, method = "spearman", interactive = TRUE){
  correlation <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{intensity}}) %>%
    tidyr::pivot_wider(names_from = {{sample}}, values_from = {{intensity}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(grouping))) %>%
    stats::cor(method = {{method}}, use = "complete.obs")

  annotation <- data %>%
    dplyr::distinct({{sample}}, {{condition}}, {{digestion}}, {{run_order}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(sample)))

  # Create heatmaply dendrogram for pheatmap
  distance <- stats::dist(correlation)
  hierachical_clustering <- stats::hclust(distance)
  dendrogram <- stats::as.dendrogram(hierachical_clustering)
  dendrogram_row <- dendextend::seriate_dendrogram(dendrogram, distance, method="OLO")
  dendrogram_column <- dendextend::rotate(dendrogram_row, order = rev(labels(distance)[seriation::get_order(stats::as.hclust(dendrogram_row))]))

  if (interactive == TRUE) {
    heatmap_interactive <-
      heatmaply::heatmaply(
        correlation,
        main = "Correlation based hirachical clustering of samples",
        col_side_colors = annotation,
        k_col = NA,
        k_row = NA,
        plot_method = "plotly"
      )
    return(heatmap_interactive)
  }
  if (interactive == FALSE) {
    heatmap_static <-
      pheatmap::pheatmap(
        correlation,
        cluster_rows = stats::as.hclust(dendrogram_row),
        cluster_cols = stats::as.hclust(dendrogram_column),
        annotation = annotation,
        main = "Correlation based hirachical clustering of samples",
        color = viridis::viridis(100)
      )
    return(heatmap_static)
  }
}
