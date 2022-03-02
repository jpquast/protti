#' Peak width over retention time
#'
#' Plots one minute binned median precursor elution peak width over retention time for each sample.
#'
#' @param data a data frame containing at least sample names and protein IDs.
#' @param sample a character column in the \code{data} data frame that contains the sample names.
#' @param intensity a numeric column in the \code{data} data frame that contains intensities. If
#' \code{remove_na_intensities = FALSE}, this argument is not required.
#' @param retention_time a numeric column in the \code{data} data frame that contains retention
#' times of precursors.
#' @param peak_width a numeric column in the \code{data} data frame that contains peak width
#' information. It is not required if \code{retention_time_start} and \code{retention_time_end}
#' columns are provided.
#' @param retention_time_start a numeric column in the \code{data} data frame that contains the
#' start time of the precursor elution peak. It is not required if the \code{peak_width} column
#' is provided.
#' @param retention_time_end a numeric column in the \code{data} data frame that contains the end
#' time of the precursor elution peak. It is not required if the \code{peak_width} column is
#' provided.
#' @param remove_na_intensities a logical value that specifies if sample/grouping combinations
#' with intensities that are NA (not quantified IDs) should be dropped from the data frame.
#' Default is TRUE since we are usually interested in the peak width of quantifiable data.
#' @param interactive a logical value that specifies whether the plot should be interactive
#' (default is FALSE).
#'
#' @return A line plot displaying one minute binned median precursor elution peak width over
#' retention time for each sample.
#' @import ggplot2
#' @importFrom stats median
#' @importFrom dplyr distinct rename mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data as_label enquo
#' @importFrom plotly ggplotly
#' @importFrom stringr str_sort
#' @importFrom utils data
#' @export
#'
#' @examples
#'
#' data <- data.frame(
#'   r_file_name = c(rep("sample_1", 10), rep("sample2", 10)),
#'   fg_quantity = c(rep(2000, 20)),
#'   eg_mean_apex_rt = c(rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 2)),
#'   eg_start_rt = c(0.5, 1, 3, 4, 5, 6, 7, 7.5, 8, 9, 1, 2, 2, 3, 4, 5, 5, 8, 9, 9),
#'   eg_end_rt = c(1.5, 2, 3.1, 4.5, 5.8, 6.6, 8, 8, 8.4, 
#'   9.1, 3, 2.2, 4, 3.4, 4.5, 5.5, 5.6, 8.3, 10, 12)
#' )
#' qc_peak_width(
#'   data,
#'   sample = r_file_name,
#'   intensity = fg_quantity,
#'   retention_time = eg_mean_apex_rt,
#'   retention_time_start = eg_start_rt,
#'   retention_time_end = eg_end_rt
#' )
qc_peak_width <- function(data,
                          sample,
                          intensity,
                          retention_time,
                          peak_width = NULL,
                          retention_time_start = NULL,
                          retention_time_end = NULL,
                          remove_na_intensities = TRUE,
                          interactive = FALSE) {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  if (remove_na_intensities == TRUE) {
    if (missing(intensity)) {
      stop(strwrap("Please provide a column containing
intensities or set remove_na_intensities to FALSE",
        prefix = "\n", initial = ""
      ))
    }

    data <- data %>%
      tidyr::drop_na({{ intensity }})
  }
  if (missing(peak_width) & missing(retention_time_start)) {
    stop(strwrap("Please provide either a peak_width column or
retention_time_start and retention_time_end columns.",
      prefix = "\n", initial = ""
    ))
  }
  if (!missing(peak_width)) {
    result <- data %>%
      dplyr::distinct(
        {{ sample }},
        {{ retention_time }},
        {{ peak_width }}
      ) %>%
      dplyr::rename(peak_width = {{ peak_width }})
  }
  if (!missing(retention_time_start)) {
    result <- data %>%
      dplyr::distinct(
        {{ sample }},
        {{ retention_time }},
        {{ retention_time_start }},
        {{ retention_time_end }}
      ) %>%
      dplyr::mutate(peak_width = {{ retention_time_end }} - {{ retention_time_start }})
  }

  peak_width_plot <- result %>%
    dplyr::mutate({{ sample }} := factor({{ sample }},
      levels = unique(stringr::str_sort({{ sample }}, numeric = TRUE))
    )) %>%
    ggplot2::ggplot(ggplot2::aes({{ retention_time }}, .data$peak_width)) +
    ggplot2::stat_summary_bin(aes(col = {{ sample }}),
      size = 1,
      geom = "line",
      binwidth = 1,
      fun = median
    ) +
    ggplot2::labs(
      title = "Median peak width over retention time",
      x = "Retention time [min]",
      y = "Median peak width [min]",
      color = "Sample"
    ) +
    ggplot2::scale_color_manual(values = rep(protti_colours, 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15)
    )

  if (interactive == FALSE) {
    return(peak_width_plot)
  }
  if (interactive == TRUE) plotly::ggplotly(peak_width_plot)
}
