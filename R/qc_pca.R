#' Plot principal component analysis
#'
#' Plots a principal component analysis based on peptide or precursor intensities.
#'
#' @param data a data frame that contains sample names, peptide or precursor identifiers,
#' corresponding intensities and a condition column indicating e.g. the treatment.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' @param grouping a character column in the \code{data} data frame that contains either precursor
#' or peptide identifiers.
#' @param intensity a numeric column in the \code{data} data frame that contains the corresponding
#' intensity values for each peptide or precursor.
#' @param condition a numeric or character column in the \code{data} data frame that contains condition information
#' (e.g. "treated" and "control").
#' @param components a character vector indicating the two components that should be displayed in
#' the plot. By default these are PC1 and PC2. You can provide these using a character vector of
#' the form c("PC1", "PC2").
#' @param digestion optional, a character column in the \code{data} data frame that indicates the
#' mode of digestion (limited proteolysis or tryptic digest). Alternatively, any other variable
#' by which the data should be split can be provided.
#' @param plot_style a character value that specifies what plot should be returned. If
#' `plot_style = "pca"` is selected the two PCA components supplied with the `components` argument
#' are plottet against each other. This is the default. `plot_style = "scree"` returns a scree
#' plot that displays the variance explained by each principal component in percent. The scree is
#' useful for checking if any other than the default first two components should be plotted.
#'
#' @return A principal component analysis plot showing PC1 and PC2. If `plot_style = "scree"`, a
#' scree plot for all dimensions is returned.
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
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 100,
#'   frac_change = 0.05,
#'   n_replicates = 3,
#'   n_conditions = 2,
#' )
#'
#' # Plot scree plot
#' qc_pca(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   condition = condition,
#'   plot_style = "scree"
#' )
#'
#' # Plot principal components
#' qc_pca(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   condition = condition
#' )
qc_pca <-
  function(data,
           sample,
           grouping,
           intensity,
           condition,
           components = c("PC1", "PC2"),
           digestion = NULL,
           plot_style = "pca") {
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
      dplyr::mutate(dimension = factor(.data$dimension,
        levels = unique(stringr::str_sort(.data$dimension, numeric = TRUE))
      ))

    if (plot_style == "pca") {
      plot <- pca_df %>%
        ggplot2::ggplot(aes(
          x = !!rlang::sym(components[1]),
          y = !!rlang::sym(components[2]),
          col = {{ condition }},
          shape = {{ digestion }}
        )) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(
          title = "Principal component analysis",
          x = paste(
            components[1],
            "(",
            round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == stringr::str_extract(components[1], "\\d")], 1),
            "%)"
          ),
          y = paste(
            components[2],
            "(",
            round(pca_sdev_df$percent_variance[pca_sdev_df$dimension == stringr::str_extract(components[2], "\\d")], 1),
            "%)"
          ),
          color = "Condition"
        ) +
        ggrepel::geom_text_repel(
          aes(label = paste(
            stringr::str_replace_all(as.character({{ sample }}), fixed("_"), " ")
          )),
          size = 4,
          show.legend = FALSE
        ) +
        {
          if (is.numeric(unique(dplyr::pull(pca_df, {{ condition }})))) {
            ggplot2::scale_color_gradientn(colours = c(
              "#0D0887", "#2E0595", "#46039F", "#5C01A6", "#7201A8", "#8707A6", "#9A169F",
              "#AC2694", "#BC3587", "#CA457A", "#D6556D", "#E26561", "#EB7655", "#F48849",
              "#FA9B3D", "#FDAF31", "#FDC527", "#F9DC24", "#F0F921"
            ))
          } else {
            ggplot2::scale_color_manual(values = protti_colours)
          }
        } +
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
        ggplot2::geom_text(
          aes(label = paste0(
            as.character(round(.data$percent_variance, digits = 1)), "%"
          )),
          size = 4,
          vjust = -0.6,
          hjust = -0.1
        ) +
        ggplot2::scale_y_continuous(
          limits = NULL,
          expand = ggplot2::expansion(mult = c(0.05, 0.08))
        ) +
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
