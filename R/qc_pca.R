#' Plot principal component analysis
#'
#' Plots a principal component analysis based on peptide or precursor intensities.
#'
#' @param data A dataframe containing  sample names, peptide or precursor identifiers, corresponding intensities and a condition column indicating e.g. the treatment.
#' @param sample The column in the data dataframe containing the sample name.
#' @param grouping The column in the data dataframe containing either precursor or peptide identifiers.
#' @param intensity Column containing the corresponding intensity values for each peptide or precursor.
#' @param condition Column indicating the treatment or condition for each sample.
#' @param digestion Optional column indicating the mode of digestion (limited proteolysis or tryptic digest).
#'
#' @return A plotted principal component analysis showing PC1 and PC2
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#' @importFrom stats prcomp
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' \dontrun{
#' qc_PCA(
#' data,
#' sample = r_file_name,
#' grouping = eg_precursor_id,
#' intensity = normalised_intensity_log2,
#' condition = r_condition
#' )
#' }
#'
qc_pca <-
  function(data, sample, grouping, intensity, condition, digestion = NULL)
  {
    . = NULL

    pca_input <- data %>%
      dplyr::distinct({{sample}}, {{grouping}}, {{intensity}}) %>%
      dplyr::group_by({{sample}}, {{grouping}}) %>%
      dplyr::summarise(intensity = sum({{intensity}})) %>%
      tidyr::pivot_wider(names_from = {{sample}}, values_from = intensity) %>%
      tidyr::drop_na() %>%
      dplyr::select(-{{grouping}}) %>%
      t(.)

    annotation <- data %>%
      dplyr::distinct({{sample}}, {{condition}}, {{digestion}})

    pca <- prcomp(pca_input, center = TRUE)

    pca_df <- as.data.frame(pca$x) %>%
      dplyr::mutate({{sample}} := factor(row.names(.))) %>%
      dplyr::left_join(annotation, by = rlang::as_name(enquo(sample)))

    pca_sdev_df <- as.data.frame(pca$sdev)

    pca_sdev_df <- pca_sdev_df %>%
      dplyr::mutate(percent_variance = (pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100),
                    dimension = row.names(.))

    plot <- pca_df %>%
      ggplot2::ggplot(aes(x = .data$PC1, y = .data$PC2, col = {{condition}}, shape = {{digestion}})) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
        title = "Principal component analysis",
        x = paste("PC1", "(", round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == 1], 1), "%)"),
        y = paste("PC2", "(", round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == 2], 1), "%)")
      ) +
      ggrepel::geom_text_repel(aes(label = paste(
        stringr::str_replace_all(as.character({{sample}}), fixed("_"), " ")
      )),
      size = 4,
      show.legend = FALSE) +
      ggplot2::theme(
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title = element_text(face = "italic"),
        title = element_text(face = "bold")
      )
    return(plot)
  }
