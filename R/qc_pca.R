#' Plot principal component analysis
#'
#' Plots a principal component analysis based on peptide or precursor intensities.
#'
#' @param data a data frame containing  sample names, peptide or precursor identifiers, corresponding intensities and a condition column indicating e.g. the treatment.
#' @param sample the column in the data data frame containing the sample name.
#' @param grouping the column in the data data frame containing either precursor or peptide identifiers.
#' @param intensity the column in the data data frame containing containing the corresponding intensity values for each peptide or precursor.
#' @param condition the column in the data data frame indicating the treatment or condition for each sample.
#' @param components character vector indicating the two components that should be displayed in the plot. By default these are PC1 and
#' PC2. You can provide these using a character vector of the form c("PC1", "PC2").
#' @param digestion optional column indicating the mode of digestion (limited proteolysis or tryptic digest).
#' @param plot_style character vector specifying what plot should be returned. If `plot_style = "pca"` is selected the two PCA
#' components supplied with the `components` argument are plottet against each other. This is the default. `plot_style = "scree"` returns
#' a scree plot that displays the variance explained by each principal component in percent. The scree is useful for checking if any other
#' than the default first two components should be plotted.
#'
#' @return A plotted principal component analysis showing PC1 and PC2
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace str_extract str_sort
#' @importFrom stats prcomp
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' qc_pca(
#'   data,
#'   sample = r_file_name,
#'   grouping = eg_precursor_id,
#'   intensity = normalised_intensity_log2,
#'   condition = r_condition,
#'   components = c("PC2", "PC3"),
#'   plot_style = "scree"
#' )
#' }
#'
qc_pca <-
  function(data, sample, grouping, intensity, condition, components = c("PC1", "PC2"), digestion = NULL, plot_style = "pca") {
    protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colours", envir = environment()) # then overwrite it with real data
    . <- NULL

    pca_input <- data %>%
      dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
      dplyr::group_by({{ sample }}, {{ grouping }}) %>%
      dplyr::summarise(intensity = sum({{ intensity }})) %>%
      tidyr::pivot_wider(names_from = {{ sample }}, values_from = intensity) %>%
      tidyr::drop_na() %>%
      dplyr::select(-{{ grouping }}) %>%
      t(.)

    annotation <- data %>%
      dplyr::distinct({{ sample }}, {{ condition }}, {{ digestion }})

    pca <- prcomp(pca_input, center = TRUE)

    pca_df <- as.data.frame(pca$x) %>%
      dplyr::mutate({{ sample }} := factor(row.names(.))) %>%
      dplyr::left_join(annotation, by = rlang::as_name(enquo(sample)))

    pca_sdev_df <- as.data.frame(pca$sdev)

    pca_sdev_df <- pca_sdev_df %>%
      dplyr::mutate(
        percent_variance = (pca$sdev^2 / sum(pca$sdev^2) * 100),
        dimension = row.names(.)
      ) %>%
      dplyr::mutate(dimension = factor(.data$dimension, levels = unique(stringr::str_sort(.data$dimension, numeric = TRUE))))
    if (plot_style == "pca") {
      plot <- pca_df %>%
        ggplot2::ggplot(aes(x = !!rlang::sym(components[1]), y = !!rlang::sym(components[2]), col = as.character({{ condition }}), shape = {{ digestion }})) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(
          title = "Principal component analysis",
          x = paste(components[1], "(", round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == stringr::str_extract(components[1], "\\d")], 1), "%)"),
          y = paste(components[2], "(", round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == stringr::str_extract(components[2], "\\d")], 1), "%)"),
          color = "Condition"
        ) +
        ggrepel::geom_text_repel(aes(label = paste(
          stringr::str_replace_all(as.character({{ sample }}), fixed("_"), " ")
        )),
        size = 4,
        show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = protti_colours) +
        ggplot2::theme(
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(face = "italic"),
          title = element_text(face = "bold"),
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 15),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15)
        )
    }
    if (plot_style == "scree") {
      plot <- pca_sdev_df %>%
        ggplot2::ggplot(aes(x = .data$dimension, y = .data$percent_variance)) +
        ggplot2::geom_col(col = "black", size = 1, fill = protti_colours[1]) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_line(size = 1, group = 1) +
        ggplot2::labs(
          title = "Principal component scree plot",
          x = "Dimension",
          y = "Explained variance [%]"
        ) +
        ggplot2::geom_text(aes(label = paste0(
          as.character(round(.data$percent_variance, digits = 1)), "%"
        )),
        size = 4,
        vjust = -0.6,
        hjust = -0.1
        ) +
        ggplot2::scale_y_continuous(limits = NULL, expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 15),
          axis.title.y = ggplot2::element_text(size = 15)
        )
    }


    return(plot)
  }
