#' Check CV distribution
#'
#' Calculates and plots the coefficients of variation for the selected grouping.
#'
#' @param data A dataframe containing at least peptide, precursor or protein identifiers, information on conditions and intensity values for each peptide, precursor or protein.
#' @param grouping The column in the input dataframe containing the grouping variables (e.g. peptides, precursors or proteins).
#' @param condition The column in the data dataframe containing condition information (e.g. "treated" and "control").
#' @param intensity Column containing the corresponding raw intensity values (not log2 transformed) for each peptide or precursor.
#' @param plot A logical indicating whether the result should be plotted. Default is TRUE.
#' @param plot_style Variable indicating the plotting style. \code{plot_style = "boxplot"} plots a boxplot, whereas \code{plot_style = "density"} plots the CV density distribution. \code{plot_style = "violin"} returns a violin plot. Default is \code{plot_style = "density"}.
#'
#' @return Either the median CVs in % or a plot showing the distribution of the CVs.
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats sd
#' @importFrom forcats fct_relevel
#' @importFrom utils data
#' @export
#'
#' @examples
#' \dontrun{
#' qc_cvs(
#' data,
#' grouping = pep_stripped_sequence,
#' condition = condition,
#' intensity = fg_quantity,
#' plot = TRUE,
#' plot_style = "density"
#' )
#' }
qc_cvs <-
  function(data, grouping, condition, intensity, plot = TRUE, plot_style = "density") {
    protti_colors <- "placeholder" # assign a placeholder to prevent a missing global variable warning
    utils::data("protti_colors", envir=environment()) # then overwrite it with real data
    if(plot == FALSE)
    {
      if(max(dplyr::pull(data, {{intensity}}), na.rm = TRUE) < 1000)
      {
      stop("Please backtransform your data or use raw values. The function does not handle log2 transformed data.")
      }

      result <- data %>%
        dplyr::distinct({{grouping}}, {{condition}}, {{intensity}}) %>%
        tidyr::drop_na({{intensity}}) %>%
        dplyr::group_by({{grouping}}) %>%
        dplyr::mutate(cv_combined = (stats::sd({{intensity}}) / mean({{intensity}})) * 100) %>%
        dplyr::group_by({{condition}}, {{grouping}}) %>%
        dplyr::mutate(cv = (stats::sd({{intensity}}) / mean({{intensity}})) * 100) %>%
        dplyr::distinct({{condition}}, {{grouping}}, .data$cv_combined, .data$cv) %>%
        tidyr::drop_na() %>%
        dplyr::group_by({{condition}}) %>%
        dplyr::mutate(median_cv = stats::median(.data$cv)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(median_cv_combined = stats::median(.data$cv_combined)) %>%
        dplyr::select(- {{grouping}}, - .data$cv_combined, - .data$cv) %>%
        dplyr::distinct()

      return(result)
    }

    if(plot == TRUE)
    {
      if(max(dplyr::pull(data, {{intensity}}), na.rm = TRUE) < 1000)
      {
        stop("Please backtransform your data or use raw values. The function does not handle log2 transformed data.")
      }
      result <- data %>%
        dplyr::distinct({{grouping}}, {{condition}}, {{intensity}}) %>%
        tidyr::drop_na({{intensity}}) %>%
        dplyr::group_by({{grouping}}) %>%
        dplyr::mutate(cv_combined = (stats::sd({{intensity}}) / mean({{intensity}})) * 100) %>%
        dplyr::group_by({{condition}}, {{grouping}}) %>%
        dplyr::mutate(cv = (stats::sd({{intensity}}) / mean({{intensity}})) * 100) %>%
        dplyr::ungroup() %>%
        dplyr::distinct({{condition}}, {{grouping}}, .data$cv_combined, .data$cv) %>%
        tidyr::drop_na() %>%
        tidyr::pivot_longer(cols = starts_with("cv"), names_to = "type", values_to = "values") %>%
        dplyr::mutate(type = ifelse(.data$type == "cv", {{condition}}, "combined")) %>%
        dplyr::mutate(type = forcats::fct_relevel(as.factor(.data$type), "combined")) %>%
        dplyr::select(-{{condition}}) %>%
        dplyr::distinct()

      if (max(result$values) > 200) {
        cv_too_high <- result %>% dplyr::filter(.data$values > 200) %>% nrow()
        warning(paste(cv_too_high), " values were exluded from the plot (CV > 200 %)")
      }

    if(plot_style == "boxplot")
    {
      plot <- ggplot2::ggplot(result) +
        ggplot2::geom_boxplot(aes(x = .data$type, y = .data$values, fill = .data$type), na.rm = TRUE) +
        ggplot2::theme_bw() +
        ggplot2::ylab("Coefficient of variation [%]") +
        ggplot2::labs(title = "Coefficients of variation") +
        ggplot2::theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
        ggplot2::scale_y_continuous(limits = c(0, 200)) +
        ggplot2::scale_fill_manual(values = protti_colors)

      return(plot)
    }
    if(plot_style == "density")
    {
      plot <- ggplot2::ggplot(result) +
        ggplot2::geom_density(aes(x = .data$values, col = .data$type), size = 1, na.rm = TRUE) +
        ggplot2::theme_bw() +
        ggplot2::ylab("Density") +
        ggplot2::labs(title = "Coefficients of variation") +
        ggplot2::xlab("Coefficient of variation [%]") +
        ggplot2::scale_x_continuous(limits = c(0, 200)) +
        ggplot2::scale_color_manual(values = protti_colors)

      return(plot)
    }
      if(plot_style == "violin")
      {
        plot <- ggplot2::ggplot(result, aes(x = .data$type, y = .data$values, fill = .data$type)) +
          ggplot2::geom_violin(na.rm = TRUE) +
          ggplot2::geom_boxplot(width = 0.15, fill = "white", na.rm = TRUE, alpha = 0.6) +
          ggplot2::theme_bw() +
          ggplot2::ylab("Coefficient of variation [%]") +
          ggplot2::labs(title = "Coefficients of variation") +
          ggplot2::theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
          ggplot2::scale_y_continuous(limits = c(0, 200)) +
          ggplot2::scale_fill_manual(values = protti_colors)
        return(plot)
      }
    }
  }
