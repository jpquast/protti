#' Volcano plot
#'
#' Plots a volcano plot for the given input.
#'
#' @param data Dataframe containing at least the input variables.
#' @param grouping Column in the data dataframe containing either precursor or peptide identifiers.
#' @param log2FC Column in the dataframe containing the log2 transfromed fold changes between two conditions.
#' @param significance Column containing the p-value or adjusted p-value for the corresponding fold changes. P-value is ideally adjusted using e.g. Benjamini-Hochberg correction.
#' @param method Method used for the plot. \code{method = "target"} highlights your protein of interest in the volcano plot, \code{method = "significant"} highlights all significantly changing entities.
#' @param protein_identifier Optional column required for \code{method = "target"}, contains protein identifiers (e.g. UniProt IDs).
#' @param target Optional argument required for \code{method = "target"}, protein identifier for your protein of interest.
#' @param title Optional argument specifying the title of the volcano plot. Default is "Volcano plot".
#' @param x_axis_label Optional argument specifying the x-axis label. Default is "log2(fold change)".
#' @param y_axis_label Optional argument specifying the y-axis label. Default is -log10(q-value)".
#' @param log2FC_cutoff Optional argument specifying the log2 transformed fold change cutoff used for assessing whether changes are significant. Default value is 1.
#' @param significance_cutoff Optional argument specifying the p-value cutoff used for assessing significance of changes. Default is 0.01.
#' @param interactive Logical, indicating whether the plot should be interactive or not. Default is \code{interactive = FALSE}.
#'
#' @return Depending on the method used a volcano plot with either highlighted target protein (\code{method = "target"}) or highlighted significant proteins (\code{method = "significant"}) is returned.
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom forcats fct_inorder
#' @importFrom rlang as_name
#' @importFrom rlang :=
#' @importFrom rlang enquo
#' @importFrom rlang ensym
#' @importFrom plotly ggplotly
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_protti(
#' data,
#' grouping = pep_stripped_sequence,
#' log2FC = log2FC,
#' significance = p_value,
#' method = "target",
#' protein_identifier = uniprot_id,
#' target = "Q9Y6K9",
#' title = "Finding Nemo",
#' x_axis_label = "log2(fold change) treated vs untreated",
#' y_axis_label = "-log10(p-value)",
#' log2FC_cutoff = 2,
#' significance_cutoff = 0.05,
#' interactive = TRUE
#' )
#' }
volcano_protti <- function(data, grouping, log2FC, significance, method, protein_identifier = NULL, target = NULL, title = "Volcano plot", x_axis_label = "log2(fold change)", y_axis_label = "-log10(q-value)", log2FC_cutoff = 1, significance_cutoff = 0.01, interactive = FALSE)
{
  if (method == "target")
  {
    data <- data %>%
      mutate(!!rlang::ensym(target) := ifelse({{protein_identifier}} == rlang::as_name(rlang::enquo(target)), TRUE, FALSE))

    plot <- data %>%
      dplyr::arrange(!!rlang::ensym(target)) %>%
      dplyr::filter(!!rlang::ensym(target) == FALSE) %>%
      ggplot2::ggplot(aes(
        x = {{log2FC}},
        y = -1 * log10({{significance}}),
        colour = !!rlang::ensym(target),
        label1 = {{protein_identifier}},
        label2 = {{grouping}}
      )) +
      geom_point() +
      geom_point(
        data = dplyr::filter(data, !!rlang::ensym(target) == TRUE),
        aes(x = {{log2FC}},
            y = -1 * log10({{significance}})),
        size = 3
      ) +
      scale_color_manual(values = c("grey60", "#5680C1")) +
      labs(
        title = title,
        x = x_axis_label,
        y = y_axis_label
      ) +
      geom_hline(yintercept = -1 * log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = log2FC_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -1 * log2FC_cutoff, linetype = "dashed") +
      theme(legend.position = "none")  +
      theme_bw() +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) + 0.5, 0)))

    if(interactive == FALSE) return(plot)
    return(plotly::ggplotly(plot))
  }
  if (method == "significant")
  {
    plot <- data %>%
      ggplot2::ggplot(aes(
        x = {{log2FC}},
        y = - 1 * log10({{significance}}),
        label1 = {{protein_identifier}},
        label2 = {{grouping}}
      )) +
      geom_point(col = "grey60") +
      geom_point(
        data = dplyr::filter(data, (( {{log2FC}} > log2FC_cutoff) & ({{significance}} < significance_cutoff) ) | ( ({{log2FC}} < -1 * log2FC_cutoff) & ({{significance}} < significance_cutoff) )),
        aes(x = {{log2FC}},
            y = -1 * log10({{significance}})),
        size = 3,
        col = "#5680C1"
      ) +
      labs(
        title = title,
        x = x_axis_label,
        y = y_axis_label
      ) +
      geom_hline(yintercept = -1 * log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = log2FC_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -1 *log2FC_cutoff, linetype = "dashed") +
      theme(legend.position = "none")  +
      theme_bw() +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{log2FC}})), na.rm = TRUE) + 0.5, 0)))

    if(interactive == FALSE) return(plot)
    return(plotly::ggplotly(plot))
  }
}
