#' Volcano plot
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `volcano_plot()`.
#'
#' @keywords internal
#' @export
volcano_protti <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "volcano_protti()",
    "volcano_plot()",
    details = "This function has been renamed."
  )
  volcano_plot(...)
}
#' Volcano plot
#'
#' Plots a volcano plot for the given input.
#'
#' @param data a data frame that contains at least the input variables.
#' @param grouping a character column in the \code{data} data frame that contains either precursor
#' or peptide identifiers.
#' @param log2FC a character column in the \code{data} data frame that contains the log2
#' transfromed fold changes between two conditions.
#' @param significance a character column in the \code{data} data frame that contains the p-value
#' or adjusted p-value for the corresponding fold changes. P-value is ideally adjusted using e.g.
#' Benjamini-Hochberg correction.
#' @param method a character value that specifies the method used for the plot.
#' \code{method = "target"} highlights your protein, proteins or any other entities of interest
#' (specified in the `target` argument) in the volcano plot. \code{method = "significant"}
#' highlights all significantly changing entities.
#' @param target_column optional, a column required for \code{method = "target"}, can contain for
#' example protein identifiers or a logical that marks certain proteins such as proteins that are
#' known to interact with the treatment. Can also be provided if \code{method = "significant"}
#' to label data points in an interactive plot.
#' @param target optional, a vector required for \code{method = "target"}. It
#' can contain one or more specific entities of the column provided in \code{target_column}. This
#' can be for example a protein ID if \code{target_column} contains protein IDs or TRUE or FALSE
#' for a logical column.
#' @param facet_by optional, a character column that contains information by which the data should
#' be faceted into multiple plots.
#' @param title optional, a character value that specifies the title of the volcano plot. Default
#' is "Volcano plot".
#' @param x_axis_label optional, a character value that specifies the x-axis label. Default is
#' "log2(fold change)".
#' @param y_axis_label optional, a character value that specifies the y-axis label. Default is
#' "-log10(q-value)".
#' @param legend_label optional, a character value that specifies the legend label. Default is
#' "Target".
#' @param log2FC_cutoff optional, a numeric value that specifies the log2 transformed fold change
#' cutoff used for assessing whether changes are significant. Default value is 1.
#' @param significance_cutoff optional, a numeric value that specifies the p-value cutoff used for
#' assessing significance of changes. Default is 0.01.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return Depending on the method used a volcano plot with either highlighted targets
#' (\code{method = "target"}) or highlighted significant proteins (\code{method = "significant"})
#' is returned.
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
#' set.seed(123) # Makes example reproducible
#'
#' # Create synthetic data
#' data <- create_synthetic_data(
#'   n_proteins = 10,
#'   frac_change = 0.5,
#'   n_replicates = 4,
#'   n_conditions = 2,
#'   method = "effect_random",
#'   additional_metadata = FALSE
#' )
#'
#' # Assign missingness information
#' data_missing <- assign_missingness(
#'   data,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity = peptide_intensity_missing,
#'   ref_condition = "all",
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' # Calculate differential abundances
#' diff <- calculate_diff_abundance(
#'   data = data_missing,
#'   sample = sample,
#'   condition = condition,
#'   grouping = peptide,
#'   intensity_log2 = peptide_intensity_missing,
#'   missingness = missingness,
#'   comparison = comparison,
#'   method = "moderated_t-test",
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' volcano_plot(
#'   data = diff,
#'   grouping = peptide,
#'   log2FC = diff,
#'   significance = pval,
#'   method = "target",
#'   target_column = change_peptide,
#'   target = TRUE,
#'   facet_by = comparison,
#'   significance_cutoff = 0.05
#' )
volcano_plot <- function(data,
                         grouping,
                         log2FC,
                         significance,
                         method,
                         target_column = NULL,
                         target = NULL,
                         facet_by = NULL,
                         title = "Volcano plot",
                         x_axis_label = "log2(fold change)",
                         y_axis_label = "-log10(q-value)",
                         legend_label = "Target",
                         log2FC_cutoff = 1,
                         significance_cutoff = 0.01,
                         interactive = FALSE) {
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
        if (!missing(facet_by)) {
          ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet_by)),
            scales = "free"
          )
        }
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
      scale_x_continuous(breaks = seq(
        round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0), 1
      )) +
      coord_cartesian(xlim = c(
        round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0)
      ))

    if (interactive == FALSE) {
      return(plot)
    }
    return(plotly::ggplotly(plot, tooltip = c("x", "y", "label1", "label2")))
  }
  if (method == "significant") {
    plot <- data %>%
      dplyr::filter(!((({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < significance_cutoff)) |
        (({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < significance_cutoff)))) %>%
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
        data = dplyr::filter(data, (({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < significance_cutoff)) |
          (({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < significance_cutoff))),
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
        if (!missing(facet_by)) {
          ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet_by)),
            scales = "free"
          )
        }
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
      scale_x_continuous(breaks = seq(
        round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0), 1
      )) +
      coord_cartesian(xlim = c(
        round(-1 * max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data, {{ log2FC }})), na.rm = TRUE) + 0.5, 0)
      ))

    if (interactive == FALSE) {
      return(plot)
    }
    return(plotly::ggplotly(plot))
  }
}
