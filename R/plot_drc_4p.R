#' plotting four-parameter dose response curves
#'
#' Function for plotting four-parameter dose response curves for each group (precursor, peptide or protein), based on output from \code{fit_drc_4p} function.
#'
#' @param data A data frame that is obtained by calling the \code{fit_drc_4p} function.
#' @param grouping The name of the column containing precursor, peptide or protein identifiers.
#' @param response The name of the column containing response values, eg. log2 transformed intensities.
#' @param dose The name of the column containing dose values, eg. the treatment concentrations.
#' @param targets A character vector that specifies the names of the precursors, peptides or proteins (depending on \code{grouping}) that should be plottet. This can also be \code{"all"} if plots for all curve fits
#' should be created. If names are provided maximally 20 proteins are plotted at a time, the rest is ignored. If more should be plotted, a mapper over asubsetted data frame should be created.
#' @param unit A character vector specifiying the unit of the concentration.
#' @param y_axis_name A character vector specifiying the name of the y-axis of the plot.
#' @param scales A character vector that specifies if the scales in facetted plots (if more than one target was provided) should be \code{"free"} or \code{"fixed"}.
#' 
#' @return If \code{targets = "all"} a list containing plots for every unique identifier in the \code{grouping} variable is created. Otherwise a plot for the specified targets is created with maximally 20 facets. 
#' @import dplyr
#' @import tidyr
#' @import progress
#' @import ggplot2
#' @importFrom purrr pmap
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
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
plot_drc_4p <- function(data, grouping, response, dose, targets, unit = paste0("\U03BC","M"), y_axis_name = "Response", scales = "free"){
  . = NULL
  
  #early filter to speed up function 
  if(!"all" %in% targets){
    data <- data %>% 
      dplyr::filter({{grouping}} %in% targets)
  }
  
  data <- data %>% 
    dplyr::mutate(name = paste0({{grouping}}, " (correlation = ", round(.data$correlation, digits = 2), ", Kd = ", round(.data$`ec_50:(Intercept)`), ")"))
  
  input_points <- data %>% 
    dplyr::select({{grouping}}, .data$name, .data$plot_points) %>% 
    tidyr::unnest(.data$plot_points)
  
  input_curve <- data %>% 
    dplyr::select({{grouping}}, .data$name, .data$plot_curve) %>% 
    tidyr::unnest(.data$plot_curve)
  
  if("all" %in% targets){
    input_points_plot <- input_points %>% 
      split(.$name)
    
    input_curve_plot <- input_curve %>% 
      split(.$name)
    
    pb <- progress::progress_bar$new(total = length(input_points_plot), format = " Plotting curves [:bar] :current/:total (:percent) :eta")
    plots <- purrr::pmap(list(x = input_points_plot, y = input_curve_plot, z = names(input_points_plot)), function(x, y, z){
      pb$tick()
      ggplot2::ggplot(data = x, ggplot2::aes(x = {{dose}}, y = {{response}})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        suppressWarnings(ggplot2::geom_ribbon(data = y, ggplot2::aes(x = .data$dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2, fill = "#B96DAD")) +
        ggplot2::geom_line(data = y, ggplot2::aes(x=dose, y = .data$Prediction), size = 1.2) +
        ggplot2::labs(title = z, x = paste0("log10 Concentration [", unit, "]"), y = y_axis_name) +
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
    }
    )
    
    return(plots)
  }
  
  if(!"all" %in% targets){
    if(length(targets) == 1){
      input_points_plot <- input_points %>% 
        dplyr::filter({{grouping}} == targets)
      
      input_curve_plot <- input_curve %>% 
        dplyr::filter({{grouping}} == targets)
      
      plot <- ggplot2::ggplot(data = input_points_plot, ggplot2::aes(x = {{dose}}, y = {{response}})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        suppressWarnings(ggplot2::geom_ribbon(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2, fill = "#B96DAD")) +
        ggplot2::geom_line(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction), size = 1.2) +
        ggplot2::labs(title = unique(input_points_plot$name), x = paste0("log10 Concentration [", unit, "]"), y = y_axis_name) +
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
      return(plot)
    } else {
      if(length(targets) > 20){
        n_targets <- length(targets)
        twenty_targets <- unique(dplyr::pull(input_points, {{grouping}}))[1:20]
        
        input_points_plot <- input_points %>% 
          dplyr::filter({{grouping}} %in% twenty_targets)
        
        input_curve_plot <- input_curve %>% 
          dplyr::filter({{grouping}} %in% twenty_targets)
        
        warning(paste("Only the first 20 targets from", rlang::as_name(enquo(grouping)),
                      "have been used for plotting since there are", n_targets,
                      "targets. Consider mapping over subsetted datasets." ))
      } else {
        input_points_plot <- input_points %>% 
          dplyr::filter({{grouping}} %in% targets)
        
        input_curve_plot <- input_curve %>% 
          dplyr::filter({{grouping}} %in% targets)
      }
      
      plot <- ggplot2::ggplot(data = input_points_plot, ggplot2::aes(x = {{dose}}, y = {{response}})) +
        ggplot2::geom_point(size = 2, col = "#5680C1") +
        suppressWarnings(ggplot2::geom_ribbon(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction, ymin = .data$Lower, ymax = .data$Upper), alpha = 0.2, fill = "#B96DAD")) +
        ggplot2::geom_line(data = input_curve_plot, ggplot2::aes(x = .data$dose, y = .data$Prediction), size = 1.2) +
        ggplot2::labs(title = "Dose response model curve fits", x = paste0("log10 Concentration [", unit, "]"), y = y_axis_name) +
        ggplot2::scale_x_log10() +
        ggplot2::facet_wrap(~ .data$name, scales = scales, ncol = 4) +
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
      return(plot)
    }
  }
}