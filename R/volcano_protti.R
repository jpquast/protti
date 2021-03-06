#' Volcano plot
#'
#' Plots a volcano plot for the given input.
#'
#' @param data a data frame containing at least the input variables.
#' @param grouping the column in the data data frame containing either precursor or peptide identifiers.
#' @param log2FC the column in the data frame containing the log2 transfromed fold changes between two conditions.
#' @param significance the column containing the p-value or adjusted p-value for the corresponding fold changes. P-value is ideally adjusted using e.g. Benjamini-Hochberg correction.
#' @param method character verctor with the method used for the plot. \code{method = "target"} highlights your protein, proteins or any other entities of
#' interest (specified in the `target` argument) in the volcano plot. \code{method = "significant"} highlights all significantly changing entities.
#' @param target_column optional column required for \code{method = "target"}, can contain for example protein identifiers or a logical that marks
#' certain proteins such as proteins that are known to interact with the treatment. Can also be provided if \code{method = "significant"}
#' to label data points in an interactive plot.
#' @param target optional character vector argument required for \code{method = "target"}. It can contain one or more specific entities of the
#' column provided in \code{target_column}. This can be for example a protein ID if \code{target_column} contains protein IDs or TRUE or FALSE for a logical column.
#' @param facet_by optional argument specifying a column that contains information by which the data should be faceted into multiple plots.
#' @param title optional argument specifying the title of the volcano plot. Default is "Volcano plot".
#' @param x_axis_label optional argument specifying the x-axis label. Default is "log2(fold change)".
#' @param y_axis_label optional argument specifying the y-axis label. Default is "-log10(q-value)".
#' @param legend_label optional argument specifying the legend label. Default is "Target".
#' @param log2FC_cutoff optional argument specifying the log2 transformed fold change cutoff used for assessing whether changes are significant. Default value is 1.
#' @param significance_cutoff optional argument specifying the p-value cutoff used for assessing significance of changes. Default is 0.01.
#' @param interactive logical, indicating whether the plot should be interactive or not. Default is \code{interactive = FALSE}.
#'
#' @return Depending on the method used a volcano plot with either highlighted targets (\code{method = "target"}) or highlighted significant proteins (\code{method = "significant"}) is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang .data new_formula enquo
#' @importFrom magrittr %>%
#' @importFrom tidyr drop_na
#' @importFrom forcats fct_inorder
#' @importFrom plotly ggplotly
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_protti(
#'   data,
#'   grouping = pep_stripped_sequence,
#'   log2FC = log2FC,
#'   significance = p_value,
#'   method = "target",
#'   target_column = uniprot_id,
#'   target = "Q9Y6K9",
#'   facet_by = comparison,
#'   title = "Finding Nemo",
#'   x_axis_label = "log2(fold change) treated vs untreated",
#'   y_axis_label = "-log10(p-value)",
#'   legend_label = "Target Protein",
#'   log2FC_cutoff = 2,
#'   significance_cutoff = 0.05,
#'   interactive = TRUE
#' )
#' }
volcano_protti <- function(data, grouping, log2FC, significance, method, target_column = NULL, target = NULL, facet_by = NULL, title = "Volcano plot", x_axis_label = "log2(fold change)", y_axis_label = "-log10(q-value)", legend_label = "Target", log2FC_cutoff = 1, significance_cutoff = 0.01, interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data

  data <- data %>%
    tidyr::drop_na({{ log2FC }}, {{ significance }})

  if (method == "target") {
    data <- data %>%
      dplyr::mutate(target = ifelse({{ target_column }} %in% target, TRUE, FALSE))

    plot <- data %>%
      dplyr::filter(.data$target == FALSE) %>%
      ggplot2::ggplot(aes(
        label1 = {{ target_column }},
        label2 = {{ grouping }}
      )) +
      geom_point(aes(
        x = {{ log2FC }},
        y = -1 * log10({{ significance }})
      ),
      colour = "grey60"
      ) +
      geom_point(
        data = dplyr::filter(data, .data$target == TRUE),
        aes(
          x = {{ log2FC }},
          y = -1 * log10({{ significance }}),
          color = {{ target_column }}
        ),
        size = 3
      ) +
      labs(
        title = title,
        x = x_axis_label,
        y = y_axis_label,
        color = legend_label
      ) +
      geom_hline(yintercept = -log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = log2FC_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -log2FC_cutoff, linetype = "dashed") +
      {
        if (!missing(facet_by)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet_by)), scales = "free")
      } +
      theme_bw() +
      theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        strip.text = ggplot2::element_text(size = 15),
        legend.title = ggplot2::element_text(size = 15),
        legend.text = ggplot2::element_text(size = 15),
        strip.background = element_blank()
      ) +
      ggplot2::scale_color_manual(values = protti_colours) +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0)))

    if (interactive == FALSE) {
      return(plot)
    }
    return(plotly::ggplotly(plot, tooltip = c("x", "y", "label1", "label2")))
  }
  if (method == "significant") {
    plot <- data %>%
      dplyr::filter(!((({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < significance_cutoff)) | (({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < significance_cutoff)))) %>%
      ggplot2::ggplot(aes(
        label1 = {{ target_column }},
        label2 = {{ grouping }}
      )) +
      geom_point(aes(
        x = {{ log2FC }},
        y = -log10({{ significance }})
      ),
      colour = "grey60"
      ) +
      geom_point(
        data = dplyr::filter(data, (({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < significance_cutoff)) | (({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < significance_cutoff))),
        aes(
          x = {{ log2FC }},
          y = -log10({{ significance }})
        ),
        size = 3,
        colour = "#5680C1"
      ) +
      labs(
        title = title,
        x = x_axis_label,
        y = y_axis_label
      ) +
      geom_hline(yintercept = -1 * log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = log2FC_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -1 * log2FC_cutoff, linetype = "dashed") +
      {
        if (!missing(facet_by)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet_by)), scales = "free")
      } +
      theme_bw() +
      theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        strip.text = ggplot2::element_text(size = 15),
        strip.background = element_blank(),
        legend.position = "none"
      ) +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0)))

    if (interactive == FALSE) {
      return(plot)
    }
    return(plotly::ggplotly(plot))
  }
}
