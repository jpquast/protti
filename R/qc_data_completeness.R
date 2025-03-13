#' Data completeness with condition-based coloring
#'
#' Calculates the percentage of data completeness with optional condition-based coloring.
#' This is a modified version of the original qc_data_completeness function.
#'
#' @param data a data frame containing at least the input variables.
#' @param sample a character or factor column in the \code{data} data frame that contains the sample names.
#' @param grouping a character column in the \code{data} data frame that contains either precursor
#' or peptide identifiers.
#' @param intensity a numeric column in the \code{data} data frame that contains any intensity
#' intensity values that missingness should be determined for.
#' @param condition optional, a character column in the \code{data} data frame that contains condition information
#' (e.g. "treated" and "control"). If provided, it will be used to color the bars.
#' @param digestion optional, a character column in the \code{data} data frame that indicates the
#' mode of digestion (limited proteolysis or tryptic digest). Alternatively, any other variable
#' by which the data should be split can be provided.
#' @param plot a logical value that indicates whether the result should be plotted.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A bar plot that displays the percentage of data completeness over all samples.
#' If \code{plot = FALSE} a data frame is returned. If \code{interactive = TRUE}, the plot is
#' interactive.
#' @import dplyr
#' @import ggplot2
#' @importFrom purrr map2_df
#' @importFrom tidyr complete
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom rlang .data := new_formula enquo
#' @importFrom stringr str_sort
#' @export
qc_data_completeness <- function(data,
                                        sample,
                                        grouping,
                                        intensity,
                                        condition = NULL,
                                        digestion = NULL,
                                        plot = TRUE,
                                        interactive = FALSE) {
  . <- NULL
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data  
  # Create a copy of the data with condition information if provided
  if (!missing(condition)) {
    data_with_condition <- data %>%
      dplyr::select({{ sample }}, {{ grouping }}, {{ intensity }}, {{ condition }}) %>%
      dplyr::distinct()
  } else {
    data_with_condition <- data %>%
      dplyr::select({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
      dplyr::distinct()
  }
  
  if (!missing(digestion)) {
    result <- data %>%
      split(dplyr::pull(data, {{ digestion }})) %>%
      purrr::map2_df(
        .y = names(.),
        .f = ~ .x %>%
          dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
          tidyr::complete({{ sample }}, {{ grouping }}) %>%
          dplyr::group_by({{ sample }}) %>%
          dplyr::summarise(completeness = sum(!is.na({{ intensity }})) / dplyr::n() * 100, .groups = "drop") %>%
          dplyr::mutate({{ digestion }} := .y)
      ) %>%
      filter(.data$completeness > 0)
  } else {
    result <- data %>%
      dplyr::distinct({{ sample }}, {{ grouping }}, {{ intensity }}) %>%
      tidyr::complete({{ sample }}, {{ grouping }}) %>%
      dplyr::group_by({{ sample }}) %>%
      dplyr::summarise(completeness = sum(!is.na({{ intensity }})) / dplyr::n() * 100, .groups = "drop")
  }
  
  # Add condition information to result if provided
  if (!missing(condition)) {
    sample_condition_map <- data_with_condition %>%
      dplyr::select({{ sample }}, {{ condition }}) %>%
      dplyr::distinct()
    
    result <- result %>%
      dplyr::left_join(sample_condition_map, by = rlang::as_name(enquo(sample)))
  }
  
  if (plot == FALSE) {
    return(result)
  }
  
  if (is(dplyr::pull(result, {{ sample }}), "character")) {
    result <- result %>%
      dplyr::mutate({{ sample }} := factor({{ sample }},
                                           levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
      ))
  }
  
  # Create the plot with condition-based coloring if condition is provided
  if (!missing(condition)) {
    completeness_plot <- result %>%
      ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$completeness, fill = {{ condition }})) +
      ggplot2::geom_col(col = "black", size = 1) +
      ggplot2::scale_fill_manual(values = protti_colours) +
      {
        if (interactive == FALSE) {
          geom_text(
            data = result,
            aes(label = round(.data$completeness, digits = 1)),
            position = position_stack(vjust = 0.5)
          )
        }
      }
  } else {
    completeness_plot <- result %>%
      ggplot2::ggplot(ggplot2::aes({{ sample }}, .data$completeness)) +
      ggplot2::geom_col(fill = "#5680C1", col = "black", size = 1) +
      {
        if (interactive == FALSE) {
          geom_text(
            data = result,
            aes(label = round(.data$completeness, digits = 1)),
            position = position_stack(vjust = 0.5)
          )
        }
      }
  }
  
  # Add common plot elements
  completeness_plot <- completeness_plot +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    {
      if (!missing(digestion)) {
        ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(digestion)), scales = "free", ncol = 2)
      }
    } +
    ggplot2::labs(
      title = "Data completeness per sample",
      x = "",
      y = "Data Completeness [%]"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
      axis.title.y = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(fill = NA),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_rect(fill = "white"),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 12)
    )
  
  if (interactive == FALSE) {
    return(completeness_plot)
  }
  
  suppressWarnings(plotly::ggplotly(completeness_plot))
} 
