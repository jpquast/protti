#' Volcano plot
#'
#' `r lifecycle::badge('deprecated')`
#' This function was deprecated due to its name changing to `volcano_plot()`.
#'
#' @return Depending on the method used a volcano plot with either highlighted targets
#' (\code{method = "target"}) or highlighted significant proteins (\code{method = "significant"})
#' is returned.
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
#' or adjusted p-value for the corresponding fold changes. The values in this column will be
#' transformed using the -log10 and displayed on the y-axis of the plot.
#' @param method a character value that specifies the method used for the plot.
#' \code{method = "target"} highlights your protein, proteins or any other entities of interest
#' (specified in the `target` argument) in the volcano plot. \code{method = "significant"}
#' highlights all significantly changing entities.
#' @param target_column optional, a column required for \code{method = "target"}, can contain for
#' example protein identifiers or a logical that marks certain proteins such as proteins that are
#' known to interact with the treatment. Can also be provided if \code{method = "significant"}
#' @param target optional, a vector required for \code{method = "target"}. It
#' can contain one or more specific entities of the column provided in \code{target_column}. This
#' can be for example a protein ID if \code{target_column} contains protein IDs or TRUE or FALSE
#' for a logical column.
#' @param split_by optional, an unquoted column name in the \code{data} data frame that contains information by which the data should
#' be split into multiple plots.
#' @param plot_ncol optional, a numeric value that specifies the number of columns in the layout
#' of multiple plots.
#' @param title optional, a character value that specifies the title of the volcano plot. Default
#' is "Volcano plot".
#' @param x_axis_label optional, a character value that specifies the x-axis label. Default is
#' "log2(fold change)".
#' @param y_axis_label optional, a character value that specifies the y-axis label. Default is
#' "-log10(q-value)".
#' @param legend_label optional, a character value that specifies the legend label. Default is
#' "Target".
#' @param colour optional, a character vector containing colours that should be used to colour
#' points according to the selected method. IMPORTANT: the first value in the vector is the
#' default point colour, the additional values specify colouring of target or significant points.
#' E.g. `c("grey60", "#5680C1")` to achieve the same colouring as the default for the "significant"
#' method. For the "significant" method, if two colors are provided after the background color,
#' they will be used for up and down regulated points respectively.
#' @param show_counts logical, whether to show the counts of up and down regulated significant points
#' in the legend. Default is FALSE.
#' @param label_top_n numeric, number of top significant entities to label in each direction (up/down)
#' based on adjusted differential values. Default is 0 (no labeling).
#' @param log2FC_cutoff optional, a numeric value that specifies the log2 transformed fold change
#' cutoff used for the vertical lines, which can be used to assess the significance of changes.
#' Default value is 1.
#' @param significance_cutoff optional, a character vector that specifies the p-value cutoff used
#' for the horizontal cutoff line, which can be used to assess the significance of changes. The
#' vector can consist solely of one element, which is the cutoff value. In that case the cutoff
#' will be applied directly to the plot. Alternatively, a second element can be provided to the
#' vector that specifies a column in the \code{data} data frame which contains e.g. adjusted
#' p-values. In that case the y-axis of the plot could display p-values that are provided to the
#' \code{significance} argument, while the horizontal cutoff line is on the scale of adjusted
#' p-values transformed to the scale of p-values. The provided vector can be e.g.
#' \code{c(0.05, "adj_pval")}. In that case the function looks for the closest adjusted p-value
#' above and below 0.05 and takes the mean of the corresponding p-values as the cutoff line. If
#' there is no adjusted p-value in the data that is below 0.05 no line is displayed. This allows
#' the user to display volcano plots using p-values while using adjusted p-values for the cutoff
#' criteria. This is often preferred because adjusted p-values are related to unadjusted p-values
#' often in a complex way that makes them hard to be interpret when plotted. Default is \code{c(0.01)}.
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
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork patchwork
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
#'   n_conditions = 3,
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
#'   method = "t-test",
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
#'   split_by = comparison,
#'   significance_cutoff = c(0.05, "adj_pval"),
#'   plot_ncol = 2
#' )
volcano_plot <- function(data,
                         grouping,
                         log2FC,
                         significance,
                         method,
                         target_column = NULL,
                         target = NULL,
                         split_by = NULL,
                         plot_ncol = NULL,
                         title = "Volcano plot",
                         x_axis_label = "log2(fold change)",
                         y_axis_label = "-log10(p-value)",
                         legend_label = "Target",
                         colour = NULL,
                         log2FC_cutoff = 1,
                         significance_cutoff = 0.01,
                         show_counts = FALSE,
                         label_top_n = 0) {
  
  protti_colours <- "placeholder"
  utils::data("protti_colours", envir = environment())
  
  if (!missing(colour)) {
    if (length(colour) < 2) {
      stop("Please provide more colours!")
    }
    background <- colour[1]
    additional_colour <- colour[-1]
  } else {
    background <- "grey60"
    additional_colour <- protti_colours
  }
  
  # Function to create a single volcano plot
  create_volcano_plot <- function(data_subset, subtitle = NULL) {
    data_subset <- data_subset %>%
      tidyr::drop_na({{ log2FC }}, {{ significance }})
    
    # Process significance cutoff for this subset
    if (length(significance_cutoff) > 1) {
      adjusted_significance <- significance_cutoff[2]
      data_subset <- data_subset %>%
        dplyr::mutate(
          centered_cutoff = !!rlang::ensym(adjusted_significance) - as.numeric(significance_cutoff[1]),
          positive = .data$centered_cutoff > 0
        ) %>%
        dplyr::mutate(is_closest_to_cutoff = abs(.data$centered_cutoff) == min(abs(.data$centered_cutoff))) %>%
        dplyr::mutate(is_closest_to_cutoff = dplyr::case_when(
          .data$positive == FALSE & .data$is_closest_to_cutoff == TRUE ~ max({{ significance }}) == {{ significance }},
          .data$positive == TRUE & .data$is_closest_to_cutoff == TRUE ~ min({{ significance }}) == {{ significance }}
        )) %>%
        dplyr::mutate(mean_adjusted_cutoff = mean({{ significance }})) %>%
        dplyr::mutate(mean_adjusted_cutoff = ifelse(all(.data$positive), NA, .data$mean_adjusted_cutoff))
    } else {
      data_subset <- data_subset %>%
        dplyr::mutate(mean_adjusted_cutoff = as.numeric(significance_cutoff[1]))
    }
    
    # Create cutoff line information
    cutoff_line <- data_subset %>%
      dplyr::mutate(mean_adjusted_cutoff = -log10(.data$mean_adjusted_cutoff)) %>%
      dplyr::distinct(.data$mean_adjusted_cutoff) %>%
      tidyr::drop_na(.data$mean_adjusted_cutoff)
    
    if (method == "target") {
      data_subset <- data_subset %>%
        dplyr::mutate(target = ifelse({{ target_column }} %in% target, TRUE, FALSE))
      
      plot <- data_subset %>%
        dplyr::filter(.data$target == FALSE) %>%
        ggplot2::ggplot(aes(
          label1 = {{ target_column }},
          label2 = {{ grouping }}
        )) +
        geom_point(
          aes(
            x = {{ log2FC }},
            y = -1 * log10({{ significance }})
          ),
          colour = background
        ) +
        geom_point(
          data = dplyr::filter(data_subset, .data$target == TRUE),
          aes(
            x = {{ log2FC }},
            y = -1 * log10({{ significance }}),
            color = {{ target_column }}
          ),
          size = 3
        )
      
      plot <- plot + 
        ggplot2::scale_color_manual(values = additional_colour)
      
    } else if (method == "significant") {
      data_subset <- data_subset %>%
        dplyr::mutate(regulation = dplyr::case_when(
          ({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < .data$mean_adjusted_cutoff) ~ "Up",
          ({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < .data$mean_adjusted_cutoff) ~ "Down",
          TRUE ~ "Not Significant"
        ))
      
      # Count significant points for this subset
      up_count <- sum(data_subset$regulation == "Up")
      down_count <- sum(data_subset$regulation == "Down")
      
      # Reverse to make color correct
      data_subset <- data_subset %>%
        dplyr::mutate(regulation = dplyr::case_when(
          ({{ log2FC }} > log2FC_cutoff) & ({{ significance }} < .data$mean_adjusted_cutoff) ~ "Down",
          ({{ log2FC }} < -log2FC_cutoff) & ({{ significance }} < .data$mean_adjusted_cutoff) ~ "Up",
          TRUE ~ "Not Significant"
        ))
      
      plot <- data_subset %>%
        dplyr::filter(.data$regulation == "Not Significant") %>%
        ggplot2::ggplot(aes(
          x = {{ log2FC }},
          y = -log10({{ significance }})
        )) +
        geom_point(colour = background) +
        geom_point(
          data = dplyr::filter(data_subset, .data$regulation != "Not Significant"),
          aes(
            x = {{ log2FC }},
            y = -log10({{ significance }}),
            color = .data$regulation
          ),
          size = 3
        )
      
      # Add labels for top entities if requested
      if (label_top_n > 0) {
        label_data <- data_subset %>%
          dplyr::filter(.data$regulation != "Not Significant") %>%
          dplyr::group_by(.data$regulation) %>%
          dplyr::arrange(dplyr::desc(abs({{ log2FC }}))) %>%
          dplyr::slice_head(n = label_top_n) %>%
          dplyr::ungroup()
        
        plot <- plot +
          ggrepel::geom_text_repel(
            data = label_data,
            aes(
              x = {{ log2FC }},
              y = -log10({{ significance }}),
              label = {{ grouping }}
            ),
            size = 3,
            box.padding = 0.5,
            point.padding = 0.1,
            force = 2,
            max.overlaps = Inf,
            min.segment.length = 0
          )
      }
      
      # Set colors and legend
      if (length(additional_colour) >= 2) {
        up_color <- additional_colour[1]
        down_color <- additional_colour[2]
      } else {
        up_color <- "#FF4136"  # Red for up-regulated
        down_color <- "#0074D9"  # Blue for down-regulated
      }
      
      # Scale color assignment is correct but should be verified with the previous definitions
      plot <- plot +
        scale_color_manual(
          values = c("Up" = up_color, "Down" = down_color),
          labels = if(show_counts) {
            c(
              paste0("Up (n=", up_count, ")"),
              paste0("Down (n=", down_count, ")")
            )
          } else {
            c("Up", "Down")
          },
          name = "Regulation"
        )
    }
    
    # Add common elements
    plot <- plot +
      labs(
        title = if (!is.null(subtitle)) subtitle else title,
        x = x_axis_label,
        y = y_axis_label
      ) +
      {
        if (nrow(cutoff_line) != 0) {
          geom_hline(yintercept = cutoff_line$mean_adjusted_cutoff, linetype = "dashed")
        }
      } +
      geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
      theme_bw() +
      theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 15),
        legend.title = ggplot2::element_text(size = 15),
        legend.text = ggplot2::element_text(size = 12),
        legend.position = if(method == "significant" && !show_counts) "none" else "right"
      ) +
      scale_x_continuous(breaks = seq(
        round(-1 * max(abs(dplyr::pull(data_subset, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data_subset, {{ log2FC }})), na.rm = TRUE) + 0.5, 0), 1
      )) +
      coord_cartesian(xlim = c(
        round(-1 * max(abs(dplyr::pull(data_subset, {{ log2FC }})), na.rm = TRUE) - 0.5, 0),
        round(max(abs(dplyr::pull(data_subset, {{ log2FC }})), na.rm = TRUE) + 0.5, 0)
      ))
    
    return(plot)
  }
  
  # If split_by is provided, create multiple plots and combine them
  if (!missing(split_by)) {
    # Get unique values from split_by column using enquo to handle unquoted variables
    split_by_quo <- rlang::enquo(split_by)
    split_values <- unique(dplyr::pull(data, !!split_by_quo))
    
    # Create list of plots
    plot_list <- list()
    for (value in split_values) {
      data_subset <- data %>%
        dplyr::filter(!!split_by_quo == value)
      
      plot_list[[as.character(value)]] <- create_volcano_plot(
        data_subset,
        subtitle = paste(title, "-", value)
      )
    }
    
    # Combine plots using patchwork
    if (length(plot_list) > 0) {
      combined_plot <- patchwork::wrap_plots(plot_list, ncol = plot_ncol %||% 1)
      return(combined_plot)
    }
  } else {
    # Create single plot
    return(create_volcano_plot(data))
  }
}
