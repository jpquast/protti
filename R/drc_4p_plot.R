#' Perform gene ontology enrichment analysis
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `drc_4p_plot()`.
#'
#' @return If \code{targets = "all"} a list containing plots for every unique identifier in the
#' \code{grouping} variable is created. Otherwise a plot for the specified targets is created with
#' maximally 20 facets.
#' @keywords internal
#' @export
plot_drc_4p <- function(...) {
  # This function has been renamed and is therefore deprecated.
  lifecycle::deprecate_warn("0.2.0",
    "plot_drc_4p()",
    "drc_4p_plot()",
    details = "This function has been renamed."
  )

  drc_4p_plot(...)
}
#' Plotting of four-parameter dose response curves
#'
#' Function for plotting four-parameter dose response curves for each group (precursor, peptide or
#' protein), based on output from \code{fit_drc_4p} function.
#'
#' @param data a data frame that is obtained by calling the \code{fit_drc_4p} function.
#' @param grouping a character column in the \code{data} data frame that contains the precursor,
#' peptide or protein identifiers.
#' @param response a numeric column in a nested data frame called \code{plot_points} that is part
#' of the \code{data} data frame. This column contains the response values, e.g. log2 transformed
#' intensities.
#' @param dose a numeric column in a nested data frame called \code{plot_points} that is part
#' of the \code{data} data frame. This column contains the dose values, e.g. the treatment
#' concentrations.
#' @param targets a character vector that specifies the names of the precursors, peptides or
#' proteins (depending on \code{grouping}) that should be plotted. This can also be \code{"all"}
#' if plots for all curve fits should be created.
#' @param unit a character value specifying the unit of the concentration.
#' @param y_axis_name a character value specifying the name of the y-axis of the plot.
#' @param facet_title_size a numeric value that specifies the size of the facet title. Default is 15.
#' @param facet a logical value that indicates if plots should be summarised into facets of 20
#' plots. This is recommended for many plots.
#' @param scales a character value that specifies if the scales in faceted plots (if more than one
#' target was provided) should be \code{"free"} or \code{"fixed"}.
#' @param x_axis_scale_log10 a logical value that indicates if the x-axis scale should be log10
#' transformed.
#' @param export a logical value that indicates if plots should be exported as PDF. The output
#' directory will be the current working directory. The name of the file can be chosen using the
#' \code{export_name} argument. If only one target is selected and \code{export = TRUE},
#' the plot is exported and in addition returned in R.
#' @param export_height a numeric value that specifies the plot height in inches for an exported plot.
#' The default is `37.5`. For a non-facet plot we recommend using 6.
#' @param export_width a numeric value that specifies the plot height in inches for an exported plot.
#' The default is `45`. For a non-facet plot we recommend using 8.
#' @param export_name a character value providing the name of the exported file if
#' \code{export = TRUE}.
#'
#' @return If \code{targets = "all"} a list containing plots for every unique identifier in the
#' \code{grouping} variable is created. Otherwise a plot for the specified targets is created with
#' maximally 20 facets.
#' @import dplyr
#' @import tidyr
#' @import progress
#' @import ggplot2
#' @importFrom forcats fct_reorder
#' @importFrom purrr pmap
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off pdf
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 2,
#'   frac_change = 1,
#'   n_replicates = 3,
#'   n_conditions = 8,
#'   method = "dose_response",
#'   concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'   additional_metadata = FALSE
#' )
#'
#' # Perform dose response curve fit
#' drc_fit <- fit_drc_4p(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   response = peptide_intensity_missing,
#'   dose = concentration,
#'   retain_columns = c(protein)
#' )
#'
#' str(drc_fit)
#'
#' # Plot dose response curves
#' if (!is.null(drc_fit)) {
#'   drc_4p_plot(
#'     data = drc_fit,
#'     grouping = peptide,
#'     response = peptide_intensity_missing,
#'     dose = concentration,
#'     targets = c("peptide_2_1", "peptide_2_3"),
#'     unit = "pM"
#'   )
#' }
#' }
drc_4p_plot <- function(data,
                        grouping,
                        response,
                        dose,
                        targets,
                        unit = "uM",
                        y_axis_name = "Response",
                        facet_title_size = 15,
                        facet = TRUE,
                        scales = "free",
                        x_axis_scale_log10 = TRUE,
                        export = FALSE,
                        export_height = 37.5,
                        export_width = 45,
                        export_name = "dose-response_curves") {
  . <- NULL

  # early filter to speed up function
  if (!"all" %in% targets) {
    data <- data %>%
      dplyr::filter({{ grouping }} %in% targets)
    if (nrow(data) == 0) stop("Target not found in data!")
  }

  data <- data %>%
    dplyr::ungroup() %>%
    dplyr::mutate(name = paste0(
      {{ grouping }},
      " (correlation = ",
      round(.data$correlation, digits = 2),
      ", EC50 = ",
      round(.data$ec_50),
      ")"
    )) %>%
    dplyr::mutate(correlation = ifelse(is.na(.data$correlation), 0, .data$correlation)) %>%
    dplyr::mutate(name = forcats::fct_reorder(.data$name, dplyr::desc(.data$correlation))) %>%
    dplyr::mutate(
      group_number = ceiling(1:dplyr::n() / 20)
    ) # we do this in preparation for faceting later.

  input_points <- data %>%
    dplyr::select({{ grouping }}, .data$name, .data$group_number, .data$plot_points) %>%
    tidyr::unnest(.data$plot_points)

  input_curve <- data %>%
    dplyr::select({{ grouping }}, .data$name, .data$group_number, .data$plot_curve) %>%
    tidyr::unnest(.data$plot_curve)

  if (!"all" %in% targets) {
    if (length(targets) == 1) {
      input_points_plot <- input_points %>%
        dplyr::filter({{ grouping }} == targets)

      input_curve_plot <- input_curve %>%
        dplyr::filter({{ grouping }} == targets)

      plot <- ggplot2::ggplot(
        data = input_points_plot,
        ggplot2::aes(
          x = {{ dose }},
          y = {{ response }}
        )
      ) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        {
          if (nrow(input_curve_plot) != 1) {
            ggplot2::geom_ribbon(
              data = input_curve_plot,
              ggplot2::aes(
                x = .data$dose,
                y = .data$Prediction,
                ymin = .data$Lower,
                ymax = .data$Upper
              ),
              alpha = 0.2,
              fill = "#B96DAD"
            )
          }
        } +
        {
          if (nrow(input_curve_plot) != 1) {
            ggplot2::geom_line(
              data = input_curve_plot,
              ggplot2::aes(
                x = .data$dose,
                y = .data$Prediction
              ),
              size = 1.2
            )
          }
        } +
        ggplot2::labs(
          title = unique(input_points_plot$name),
          x = paste0("Concentration [", unit, "]"),
          y = y_axis_name
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 18
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

      if (export == FALSE) {
        return(plot)
      } else {
        grDevices::pdf(
          file = paste0(export_name, ".pdf"),
          width = 8,
          height = 6
        )
        suppressWarnings(print(plot))
        grDevices::dev.off()
        return(plot)
      }
    } else {
      input_points <- input_points %>%
        dplyr::filter({{ grouping }} %in% targets)

      input_curve <- input_curve %>%
        dplyr::filter({{ grouping }} %in% targets)
    }
  }

  if (facet == TRUE) {
    input_points_plot <- input_points %>%
      split(.$group_number)

    input_curve_plot <- input_curve %>%
      split(.$group_number)
  } else {
    input_points_plot <- input_points %>%
      split(.$name)

    input_curve_plot <- input_curve %>%
      split(.$name)
  }
  pb <- progress::progress_bar$new(
    total = length(input_points_plot),
    format = " Plotting curves [:bar] :current/:total (:percent) :eta"
  )
  plots <- purrr::pmap(
    list(
      x = input_points_plot,
      y = input_curve_plot,
      z = names(input_points_plot)
    ),
    function(x, y, z) {
      pb$tick()
      ggplot2::ggplot(data = x, ggplot2::aes(x = {{ dose }}, y = {{ response }})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        {
          if (nrow(y) != 1) {
            ggplot2::geom_ribbon(
              data = y,
              ggplot2::aes(
                x = .data$dose,
                y = .data$Prediction,
                ymin = .data$Lower,
                ymax = .data$Upper
              ),
              alpha = 0.2,
              fill = "#B96DAD"
            )
          }
        } +
        {
          if (nrow(y) != 1) {
            ggplot2::geom_line(
              data = y,
              ggplot2::aes(
                x = .data$dose,
                y = .data$Prediction
              ),
              size = 1.2
            )
          }
        } +
        {
          if (facet == FALSE) {
            ggplot2::labs(
              title = z,
              x = paste0("Concentration [", unit, "]"),
              y = y_axis_name
            )
          }
        } +
        {
          if (facet == TRUE) {
            ggplot2::labs(
              title = "Dose-response curves",
              x = paste0("Concentration [", unit, "]"),
              y = y_axis_name
            )
          }
        } +
        {
          if (x_axis_scale_log10 == TRUE) ggplot2::scale_x_log10()
        } +
        {
          if (facet == TRUE) ggplot2::facet_wrap(~ .data$name, scales = scales, ncol = 4)
        } +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 18
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
            size = facet_title_size
          ),
          strip.background = element_blank()
        )
    }
  )
  if (export == FALSE) {
    plots
  } else {
    if (facet == TRUE) {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = export_width,
        height = export_height
      )
      suppressWarnings(print(plots))
      grDevices::dev.off()
    } else {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = export_width,
        height = export_height
      )
      suppressWarnings(print(plots))
      grDevices::dev.off()
    }
  }
}
