#' Volcano plot
#'
#' Plots a volcano plot for the given input.
#'
#' @param data Dataframe containing at least the input variables.
#' @param grouping Column in the data dataframe containing either precursor or peptide identifiers.
#' @param foldchange Column in the dataframe containing the fold changes (ideally log2 transformed) between two conditions
#' @param significance Column containing the p-value or adjusted p-value for the corresponding fold changes. P-value is ideally adjusted using e.g. Benjamini-Hochberg correction.
#' @param method Method used for the plot. \code{method = "target"} highlights your protein of interest in the volcano plot, \code{method = "significant"} highlights all significantly changing entities.
#' @param protein_identifier Optional column required for \code{method = "target"}, contains protein identifiers (e.g. UniProt IDs).
#' @param target Optional argument required for \code{method = "target"}, protein identifier for your protein of interest.
#' @param title Optional argument specifying the title of the volcano plot. Default is "Volcano plot".
#' @param x_axis_label Optional argument specifying the x-axis label. Default is "log2(fold change)".
#' @param foldchange_cutoff Optional argument specifying the fold change cutoff used for assessing whether changes are significant. As a log2 scale should be used, the default value is 1.
#' @param significance_cutoff Optional argument specifying the p-value cutoff used for assessing significance of changes. Default is 0.01.
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
#' foldchange = Log2FC,
#' significance = p_value,
#' method = "target",
#' protein_identifier = uniprot_id,
#' target = "Q9Y6K9",
#' title = "Finding Nemo",
#' x_axis_label = "log2(fold change) treated vs untreated",
#' foldchange_cutoff = 0.5,
#' significance_cutoff = 0.05
#' )
#' }
volcano_protti <- function(data, grouping, foldchange, significance, method, protein_identifier = NULL, target = NULL, title = "Volcano plot", x_axis_label = "log2(fold change)", foldchange_cutoff = 1, significance_cutoff = 0.01)
{
  if (method == "target")
  {
    data <- data %>%
      mutate(!!rlang::ensym(target) := ifelse({{protein_identifier}} == rlang::as_name(rlang::enquo(target)), TRUE, FALSE))

    plot <- data %>%
      dplyr::arrange(!!rlang::ensym(target)) %>%
      dplyr::filter(!!rlang::ensym(target) == FALSE) %>%
      ggplot2::ggplot(aes(
        x = {{foldchange}},
        y = -1 * log10({{significance}}),
        colour = !!rlang::ensym(target),
        label1 = {{protein_identifier}},
        label2 = {{grouping}}
      )) +
      geom_point() +
      geom_point(
        data = dplyr::filter(data, !!rlang::ensym(target) == TRUE),
        aes(x = {{foldchange}},
            y = -1 * log10({{significance}})),
        size = 3
      ) +
      scale_color_manual(values = c("grey60", "#5680C1")) +
      labs(
        title = title,
        x = x_axis_label,
        y = "-log10(p-value)"
      ) +
      geom_hline(yintercept = -1 * log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = foldchange_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -1 * foldchange_cutoff, linetype = "dashed") +
      theme(legend.position = "none")  +
      theme_bw() +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) + 0.5, 0)))
    return(plotly::ggplotly(plot))
  }
  if (method == "significant")
  {
    plot <- data %>%
      ggplot2::ggplot(aes(
        x = {{foldchange}},
        y = - 1 * log10({{significance}}),
        label1 = {{protein_identifier}},
        label2 = {{grouping}}
      )) +
      geom_point(col = "grey60") +
      geom_point(
        data = dplyr::filter(data, (( {{foldchange}} > foldchange_cutoff) & ({{significance}} < significance_cutoff) ) | ( ({{foldchange}} < -1 * foldchange_cutoff) & ({{significance}} < significance_cutoff) )),
        aes(x = {{foldchange}},
            y = -1 * log10({{significance}})),
        size = 3,
        col = "#5680C1"
      ) +
      labs(
        title = title,
        x = x_axis_label,
        y = "-log10(p-value)"
      ) +
      geom_hline(yintercept = -1 * log10(significance_cutoff), linetype = "dashed") +
      geom_vline(xintercept = foldchange_cutoff, linetype = "dashed") +
      geom_vline(xintercept = -1 *foldchange_cutoff, linetype = "dashed") +
      theme(legend.position = "none")  +
      theme_bw() +
      scale_x_continuous(breaks = seq(round(-1 * max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) + 0.5, 0), 1)) +
      coord_cartesian(xlim = c(round(-1 * max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) - 0.5, 0), round(max(abs(dplyr::pull(data, {{foldchange}})), na.rm = TRUE) + 0.5, 0)))
    return(plotly::ggplotly(plot))
  }
}
