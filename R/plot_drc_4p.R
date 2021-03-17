#' Plotting of four-parameter dose response curves
#'
#' Function for plotting four-parameter dose response curves for each group (precursor, peptide or protein), based on output from \code{fit_drc_4p} function.
#'
#' @param data A data frame that is obtained by calling the \code{fit_drc_4p} function.
#' @param grouping The name of the column containing precursor, peptide or protein identifiers.
#' @param response The name of the column containing response values, eg. log2 transformed intensities.
#' @param dose The name of the column containing dose values, eg. the treatment concentrations.
#' @param targets A character vector that specifies the names of the precursors, peptides or proteins (depending on \code{grouping}) 
#' that should be plotted. This can also be \code{"all"} if plots for all curve fits should be created. 
#' @param unit A character vector specifying the unit of the concentration.
#' @param y_axis_name A character vector specifying the name of the y-axis of the plot.
#' @param facet A logical indicating if plots should be summarised into facets of 20 plots. This is recommended for many plots.
#' @param scales A character vector that specifies if the scales in faceted plots (if more than one target was provided) should be \code{"free"} or \code{"fixed"}.
#' @param x_axis_scale_log10 A logical indicating if the x-axis scale should be log10 transformed. 
#' @param export A logical indicating if plots should be exported as PDF. The output directory will be the current working directory. The 
#' name of the file can be chosen using the \code{export_name} argument. If only one target is selected and \code{export = TRUE}, 
#' the plot is exported and in addition returned in R. 
#' @param export_name A character vector providing the name of the exported file if \code{export = TRUE}. 
#' 
#' @return If \code{targets = "all"} a list containing plots for every unique identifier in the \code{grouping} variable is created. Otherwise a plot for the specified targets is created with maximally 20 facets. 
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
#' \dontrun{
#' plot_drc_4p(
#' data,
#' grouping = eg_precursor_id,
#' response = intensity,
#' dose = concentration,
#' targets = c("ABCDEFK")
#' )
#' }
plot_drc_4p <- function(data, grouping, response, dose, targets, unit = "uM", y_axis_name = "Response", facet = TRUE, scales = "free", x_axis_scale_log10 = TRUE, export = FALSE, export_name = "dose-response_curves"){
  . = NULL
  
  #early filter to speed up function 
  if(!"all" %in% targets){
    data <- data %>% 
      dplyr::filter({{grouping}} %in% targets)
    if(nrow(data) == 0) stop("Target not found in data!")
  }
  
  data <- data %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(name = paste0({{grouping}}, " (correlation = ", round(.data$correlation, digits = 2), ", EC50 = ", round(.data$ec_50), ")")) %>%
    dplyr::mutate(name = forcats::fct_reorder(.data$name, desc(.data$correlation))) %>% 
    dplyr::mutate(group_number = 1,
                  group_number = ceiling(cumsum(.data$group_number)/20)) # we do this in preparation for faceting later.
    
  input_points <- data %>% 
    dplyr::select({{grouping}}, .data$name, .data$group_number, .data$plot_points) %>% 
    tidyr::unnest(.data$plot_points)
  
  input_curve <- data %>% 
    dplyr::select({{grouping}}, .data$name, .data$group_number, .data$plot_curve) %>% 
    tidyr::unnest(.data$plot_curve)
  
  if(!"all" %in% targets){
    if(length(targets) == 1){
      input_points_plot <- input_points %>% 
        dplyr::filter({{grouping}} == targets)
      
      input_curve_plot <- input_curve %>% 
        dplyr::filter({{grouping}} == targets)
      
      plot <- ggplot2::ggplot(data = input_points_plot, ggplot2::aes(x = {{dose}}, y = {{response}})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        ggplot2::geom_ribbon(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2, fill = "#B96DAD") +
        ggplot2::geom_line(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction), size = 1.2) +
        ggplot2::labs(title = unique(input_points_plot$name), x = paste0("Concentration [", unit, "]"), y = y_axis_name) +
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
      
      if(export == FALSE){
      return(plot)
      } else {
        grDevices::pdf(file = paste0(export_name, ".pdf"),
            width = 8,
            height = 6)
        suppressWarnings(print(plot))
        grDevices::dev.off()
        return(plot)
      }
      
    } else {
      input_points <- input_points %>% 
        dplyr::filter({{grouping}} %in% targets) 
      
      input_curve <- input_curve %>% 
        dplyr::filter({{grouping}} %in% targets) 
    }
  }

  if(facet == TRUE){
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
    pb <- progress::progress_bar$new(total = length(input_points_plot), format = " Plotting curves [:bar] :current/:total (:percent) :eta")
    plots <- purrr::pmap(list(x = input_points_plot, y = input_curve_plot, z = names(input_points_plot)), function(x, y, z){
      pb$tick()
      ggplot2::ggplot(data = x, ggplot2::aes(x = {{dose}}, y = {{response}})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        ggplot2::geom_ribbon(data = y, ggplot2::aes(x = .data$dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2, fill = "#B96DAD") +
        ggplot2::geom_line(data = y, ggplot2::aes(x=dose, y = .data$Prediction), size = 1.2) +
        {if(facet == FALSE) ggplot2::labs(title = z, x = paste0("Concentration [", unit, "]"), y = y_axis_name)} +
        {if(facet == TRUE) ggplot2::labs(title = "Dose-response curves", x = paste0("Concentration [", unit, "]"), y = y_axis_name)} +
        {if(x_axis_scale_log10 == TRUE) ggplot2::scale_x_log10()} +
        {if(facet == TRUE) ggplot2::facet_wrap(~ .data$name, scales = scales, ncol = 4)} +
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
    }
    )
    if(export == FALSE){
    plots
    } else {
      if(facet == TRUE){
        grDevices::pdf(file = paste0(export_name, ".pdf"),
          width = 45,
          height = 37.5)
      suppressWarnings(print(plots))
      grDevices::dev.off()
      } else {
        grDevices::pdf(file = paste0(export_name, ".pdf"),
            width = 8,
            height = 6)
        suppressWarnings(print(plots))
        grDevices::dev.off()
      }
    }
}