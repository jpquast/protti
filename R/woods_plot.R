#' Woods' plot
#'
#' Creates a Woods' plot that plots log2 fold change of peptides or precursors along the protein
#' sequence. The peptides or precursors are located on the x-axis based on their start and end
#' positions. The position on the y-axis displays the fold change. The vertical size (y-axis) of
#' the box representing the peptides or precursors do not have any meaning.
#'
#' @param data a data frame that contains differential abundance, start and end peptide or
#' precursor positions, protein length and optionally a variable based on which peptides or
#' precursors should be coloured.
#' @param fold_change a numeric column in the \code{data} data frame that contains log2 fold
#' changes.
#' @param start_position a numeric column in the \code{data} data frame that contains the start
#' positions for each peptide or precursor.
#' @param end_position a numeric column in the \code{data} data frame that contains the end
#' positions for each peptide or precursor.
#' @param protein_length a numeric column in the \code{data} data frame that contains the length
#' of the protein.
#' @param coverage optional, a numeric column in the \code{data} data frame that contains coverage
#' in percent. Will appear in the title of the Woods' plot if provided.
#' @param protein_id a character column in the \code{data} data frame that contains protein
#' identifiers.
#' @param targets a character vector that specifies the identifiers of the proteins (depending on
#' \code{protein_id}) that should be plotted. This can also be \code{"all"} if plots for all
#' proteins should be created. Default is \code{"all"}.
#' @param facet a logical value that indicates if plots should be summarised into facets of 20
#' plots. This is recommended for many plots. Default is \code{facet = TRUE}.
#' @param colouring optional, a character or numeric (discrete or continous) column in the data
#' frame containing information by which peptide or precursors should be coloured.
#' @param fold_change_cutoff optional, a numeric value that specifies the log2 fold change cutoff
#' used in the plot. The default value is 2.
#' @param highlight optional, a logical column that specifies whether specific peptides or
#' precursors should be highlighted with an asterisk.
#' @param export a logical value that indicates if plots should be exported as PDF. The output
#' directory will be the current working directory. The name of the file can be chosen using the
#' \code{export_name} argument. Default is \code{export = FALSE}.
#' @param export_name a character vector that provides the name of the exported file if
#' \code{export = TRUE}. Default is \code{export_name = "woods_plots"}
#'
#' @return A list containing Woods' plots is returned. Plotting peptide or precursor log2 fold
#' changes along the protein sequence.
#' @import ggplot2
#' @import tidyr
#' @import progress
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @importFrom utils data
#' @export
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   fold_change = c(2.3, 0.3, -0.4, -4, 1),
#'   pval = c(0.001, 0.7, 0.9, 0.003, 0.03),
#'   start = c(20, 30, 45, 90, 140),
#'   end = c(33, 40, 64, 100, 145),
#'   protein_length = c(rep(150, 5)),
#'   protein_id = c(rep("P1", 5))
#' )
#'
#' # Plot Woods' plot
#' woods_plot(
#'   data = data,
#'   fold_change = fold_change,
#'   start_position = start,
#'   end_position = end,
#'   protein_length = protein_length,
#'   protein_id = protein_id,
#'   colouring = pval
#' )
woods_plot <- function(data,
                       fold_change,
                       start_position,
                       end_position,
                       protein_length,
                       coverage = NULL,
                       protein_id,
                       targets = "all",
                       facet = TRUE,
                       colouring = NULL,
                       fold_change_cutoff = 1,
                       highlight = NULL,
                       export = FALSE,
                       export_name = "woods_plots") {
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  . <- NULL
  # create variable that contains information about the missingness of
  # certain variables so it can be used in the mapping through which the plot is created.
  colouring_missing <- missing(colouring)
  highlight_missing <- missing(highlight)

  if (!missing(highlight) && !all(is.logical(dplyr::pull(data, {{ highlight }})))) {
    stop("Please only provide logicals (i.e. TRUE or FALSE) in the 'highlight' column.")
  }

  # early filter to speed up function
  if (!"all" %in% targets) {
    data <- data %>%
      dplyr::ungroup() %>%
      dplyr::filter({{ protein_id }} %in% targets)
    if (nrow(data) == 0) stop("Target not found in data!")
  }

  data <- data %>%
    dplyr::ungroup() %>%
    tidyr::drop_na({{ protein_id }}) %>%
    dplyr::mutate({{ protein_id }} := factor({{ protein_id }},
      levels = unique(dplyr::pull(data, {{ protein_id }}))
    )) %>%
    dplyr::group_by({{ protein_id }}) %>%
    dplyr::mutate(
      group_number = dplyr::cur_group_id(),
      name = {{ protein_id }}
    ) %>%
    # we do this in preparation for faceting later.
    dplyr::ungroup() %>%
    dplyr::mutate(group_number = ceiling(.data$group_number / 20))

  # Add coverage to protein ID name if present.
  if (!missing(coverage)) {
    data <- data %>%
      dplyr::mutate(name = paste0({{ protein_id }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }

  if (facet == TRUE) {
    data_facet <- data %>%
      split(.$group_number)
  } else {
    data_facet <- data %>%
      split(.$name)
  }

  pb <- progress::progress_bar$new(
    total = length(data_facet),
    format = " Creating Woods' plots [:bar] :current/:total (:percent) :eta"
  )

  plots <- purrr::map2(.x = data_facet, .y = names(data_facet), function(x, y) {
    pb$tick()
    ggplot2::ggplot(data = x) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = 0,
          xmax = {{ protein_length }},
          ymin = -0.01,
          ymax = 0.01
        ),
        fill = "black"
      ) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = {{ start_position }},
          xmax = {{ end_position }},
          ymin = {{ fold_change }} - 0.2,
          ymax = {{ fold_change }} + 0.2,
          fill = {{ colouring }}
        ),
        col = "black",
        size = 0.7,
        alpha = 0.8
      ) +
      {
        if (!highlight_missing) {
          ggplot2::geom_point(
            data = dplyr::filter(x, {{ highlight }} == TRUE),
            ggplot2::aes(
              x = (({{ start_position }} + {{ end_position }}) / 2),
              y = ({{ fold_change }} - 0.3)
            ),
            shape = 8,
            col = "black",
            size = 3
          )
        }
      } +
      ggplot2::geom_hline(
        yintercept = -{{ fold_change_cutoff }},
        col = "blue",
        alpha = .8,
        linewidth = 0.7
      ) +
      ggplot2::geom_hline(
        yintercept = {{ fold_change_cutoff }},
        col = "blue",
        alpha = .8,
        linewidth = 0.7
      ) +
      ggplot2::ylim(
        min(-2.5, dplyr::pull(data, {{ fold_change }})) - 0.5,
        max(2.5, dplyr::pull(data, {{ fold_change }})) + 0.5
      ) +
      ggplot2::scale_x_continuous(limits = NULL, expand = c(0, 0)) +
      {
        if (facet == FALSE) {
          ggplot2::labs(
            title = y,
            x = "Protein Sequence",
            y = "log2(fold change)"
          )
        }
      } +
      {
        if (facet == TRUE) {
          ggplot2::labs(
            title = "Wood's plots",
            x = "Protein Sequence",
            y = "log2(fold change)"
          )
        }
      } +
      {
        if (facet == TRUE) ggplot2::facet_wrap(~ .data$name, scales = "free", ncol = 4)
      } +
      ggplot2::guides(size = "none") +
      {
        if (!colouring_missing && !is.numeric(dplyr::pull(data, {{ colouring }}))) {
          ggplot2::scale_fill_manual(values = protti_colours)
        }
      } +
      {
        if (!colouring_missing && is.numeric(dplyr::pull(data, {{ colouring }}))) {
          ggplot2::scale_fill_gradient(low = protti_colours[1], high = protti_colours[2])
        }
      } +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 20
        ),
        axis.text.x = ggplot2::element_text(
          size = 15
        ),
        axis.text.y = ggplot2::element_text(
          size = 15
        ),
        axis.title.y = ggplot2::element_text(
          size = 15
        ),
        axis.title.x = ggplot2::element_text(
          size = 15
        ),
        legend.title = ggplot2::element_text(
          size = 15
        ),
        legend.text = ggplot2::element_text(
          size = 15
        ),
        strip.text.x = ggplot2::element_text(
          size = 15
        ),
        strip.text = ggplot2::element_text(
          size = 15
        ),
        strip.background = element_blank()
      )
  })

  if (export == FALSE) {
    plots
  } else {
    if (facet == TRUE) {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = 45,
        height = 37.5
      )
      suppressWarnings(print(plots))
      grDevices::dev.off()
    } else {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = 8,
        height = 6
      )
      suppressWarnings(print(plots))
      grDevices::dev.off()
    }
  }
}
