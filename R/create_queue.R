#' Creates a mass spectrometer queue for Xcalibur
#'
#' This function creates a measurement queue for sample acquisition for the software Xcalibur. All possible combinations of the
#' provided information will be created to make file and sample names.
#'
#' @param date Optional, start date of the measurements.
#' @param instrument Optional, instrument initials.
#' @param user Optional, user name.
#' @param measurement_type Optional, the measurement type of the samples (e.g "DIA", "DDA", "library" etc.).
#' @param experiment_name Optional, name of the experiment.
#' @param digestion Optional, the digestion types used in this experiment (e.g "LiP" and/or "tryptic control").
#' @param treatment_type_1 Optional, name of the treatment.
#' @param treatment_type_2 Optional, name of a second treatment that was combined with the first treatment.
#' @param treatment_dose_1 Optional, doses used for treatment 1. These can be concentrations or times etc.
#' @param treatment_dose_2 Optional, doses used for treatment 2. These can be concentrations or times etc.
#' @param treatment_unit_1 Optional, unit of the doses for treatment 1 (e.g min, mM, etc.).
#' @param treatment_unit_2 Optional, unit of the doses for treatment 2 (e.g min, mM, etc.).
#' @param n_replicates Optional, number of replicates used per sample.
#' @param number_runs Specifies if file names should be numbered from 1:n instead of adding experiment information. Default is FALSE.
#' @param organism Optional, name of the organism used.
#' @param exclude_combinations Optional, list of lists that contains vectors treatment types and treatment doses whos combinations should be
#' excluded from the final queue.
#' @param inj_vol The volume used for injection in microliter. Will be \code{NA} if not specified. Then it needs to be manually specified before the queue can be used.
#' @param data_path The file path where the MS raw data should be saved. Backslashes should be escaped by another backslash. Will be \code{NA} if not specified, but needs to be specified later on then.
#' @param method_path The file path of the MS acquisition method. Backslashes should be escaped by another backslash. Will be \code{NA} if not specified, but needs to be specified later on then.
#' @param position_row The row positions that can be used for the samples (e.g c("A", "B")). If the number of specified rows and columns does not equal the total number
#' of samples, positions will be repeated.
#' @param position_column The column positions that can be used for the samples (e.g 8). If the number of specified rows and columns does not equal the total number
#' of samples, positions will be repeated.
#' @param blank_every_n Optional, specifies in which intervals a blank sample should be inserted.
#' @param blank_position The plate position of the blank. Will be \code{NA} if not specified, but needs to be specified later on then.
#' @param blank_method_path The file path of the MS acquisition method of the blank. Backslashes should be escaped by another backslash.
#' Will be \code{NA} if not specified, but needs to be specified later on then.
#' @param blank_inj_vol Injection volume of the blank sample. Will be \code{NA} if not specified, but needs to be specified later on then.
#' @param export Logical, specifying if queue should be exported from R and saved as a .csv file. Default is TRUE. Further options for
#' export can be adjusted with the \code{export_to_queue} and \code{queue_path} arguments.
#' @param export_to_queue Logical, specifying if the resulting queue should be appended to an already existing queue. If false result will be saved as \code{queue.csv}.
#' @param queue_path Optional, file path to a queue file to which the generated queue should be appended if \code{export_to_queue = TRUE}.
#' If not specified queue file can be chosen interactively.
#'
#' @return If \code{export_to_queue = FALSE} a file named \code{queue.csv} will be returned that contains the generated queue. If \code{export_to_queue = TRUE}, the resulting
#' generated queue will be appended to an already existing queue that needs to be specified either interactively or through the argument \code{queue_path}.
#' @importFrom tidyr crossing unite
#' @importFrom dplyr mutate select bind_cols
#' @importFrom data.table fread fwrite
#' @importFrom purrr map_dfr map2_dfr
#' @importFrom tibble rowid_to_column
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' create_queue(
#'   date = c("200722"),
#'   instrument = c("EX1"),
#'   user = c("jquast"),
#'   measurement_type = c("DIA"),
#'   experiment_name = c("JPQ031"),
#'   digestion = c("LiP", "tryptic control"),
#'   treatment_type_1 = c("EDTA", "H2O"),
#'   treatment_type_2 = c("Zeba", "unfiltered"),
#'   treatment_dose_1 = c(10, 30, 60),
#'   treatment_unit_1 = c("min"),
#'   n_replicates = 4,
#'   number_runs = FALSE,
#'   organism = c("E. coli"),
#'   exclude_combinations = list(list(
#'     treatment_type_1 = c("H2O"),
#'     treatment_type_2 = c("Zeba", "unfiltered"),
#'     treatment_dose_1 = c(10, 30)
#'   )),
#'   inj_vol = c(2),
#'   data_path = "D:\\2007_Data",
#'   method_path = "C:\\Xcalibur\\methods\\jquast\\DIA_120min_41var_AGC200",
#'   position_row = c("A", "B", "C", "D", "E", "F"),
#'   position_column = 8,
#'   blank_every_n = 4,
#'   blank_position = "1-V1",
#'   blank_method_path = "C:\\Xcalibur\\methods\\blank"
#' )
#' }
create_queue <-
  function(date = NULL,
           instrument = NULL,
           user = NULL,
           measurement_type = NULL,
           experiment_name = NULL,
           digestion = NULL,
           treatment_type_1 = NULL,
           treatment_type_2 = NULL,
           treatment_dose_1 = NULL,
           treatment_dose_2 = NULL,
           treatment_unit_1 = NULL,
           treatment_unit_2 = NULL,
           n_replicates = NULL,
           number_runs = FALSE,
           organism = NULL,
           exclude_combinations = NULL,
           inj_vol = NA,
           data_path = NA,
           method_path = NA,
           position_row = NA,
           position_column = NA,
           blank_every_n = NULL,
           blank_position = NA,
           blank_method_path = NA,
           blank_inj_vol = 1,
           export = TRUE,
           export_to_queue = FALSE,
           queue_path = NULL) {
    data <-
      tidyr::crossing(
        date,
        instrument,
        user,
        measurement_type,
        experiment_name,
        digestion,
        treatment_type_2,
        treatment_type_1,
        treatment_dose_1,
        treatment_dose_2,
        treatment_unit_1,
        treatment_unit_2,
        1:n_replicates
      )

    sample_name <-
      tidyr::crossing(
        organism,
        digestion,
        treatment_type_2,
        treatment_type_1,
        treatment_dose_1,
        treatment_dose_2,
        treatment_unit_1,
        treatment_unit_2,
        1:n_replicates
      )

    if (!is.null(exclude_combinations)) {
      exclude <- purrr::map_dfr(
        .x = exclude_combinations,
        .f = ~ tidyr::crossing(
          .x$treatment_type_2,
          .x$treatment_type_1,
          .x$treatment_dose_1,
          .x$treatment_dose_2
        )
      )

      data <- data %>%
        dplyr::mutate(
          t1_null = is.null(treatment_type_1),
          t2_null = is.null(treatment_type_2),
          d1_null = is.null(treatment_dose_1),
          d2_null = is.null(treatment_dose_2)
        ) %>%
        dplyr::filter(!(
          ifelse(
            !.data$t2_null,
            treatment_type_2 %in% exclude$`.x$treatment_type_2`,
            TRUE
          ) &
            ifelse(
              !.data$t1_null,
              treatment_type_1 %in% exclude$`.x$treatment_type_1`,
              TRUE
            ) &
            ifelse(
              !.data$d2_null,
              treatment_dose_2 %in% exclude$`.x$treatment_dose_2`,
              TRUE
            ) &
            ifelse(
              !.data$d1_null,
              treatment_dose_1 %in% exclude$`.x$treatment_dose_1`,
              TRUE
            )
        )) %>%
        dplyr::select(-c(.data$t1_null, .data$t2_null, .data$d1_null, .data$d2_null))

      sample_name <- sample_name %>%
        dplyr::mutate(
          t1_null = is.null(treatment_type_1),
          t2_null = is.null(treatment_type_2),
          d1_null = is.null(treatment_dose_1),
          d2_null = is.null(treatment_dose_2)
        ) %>%
        dplyr::filter(!(
          ifelse(
            !.data$t2_null,
            treatment_type_2 %in% exclude$`.x$treatment_type_2`,
            TRUE
          ) &
            ifelse(
              !.data$t1_null,
              treatment_type_1 %in% exclude$`.x$treatment_type_1`,
              TRUE
            ) &
            ifelse(
              !.data$d2_null,
              treatment_dose_2 %in% exclude$`.x$treatment_dose_2`,
              TRUE
            ) &
            ifelse(
              !.data$d1_null,
              treatment_dose_1 %in% exclude$`.x$treatment_dose_1`,
              TRUE
            )
        )) %>%
        dplyr::select(-c(.data$t1_null, .data$t2_null, .data$d1_null, .data$d2_null))
    }

    if (number_runs == TRUE) {
      data <- data %>%
        tibble::rowid_to_column("number") %>%
        dplyr::select(
          date,
          instrument,
          user,
          measurement_type,
          experiment_name,
          .data$number
        )
    }

    data <- data %>%
      tidyr::unite("File Name", sep = "_")

    sample_name <- sample_name %>%
      tidyr::unite("Sample Name", sep = " ")

    nrow_data <- nrow(data)

    position <- tidyr::crossing(position_row, 1:position_column) %>%
      tidyr::unite("Position", sep = "")

    nrow_position <- nrow(position)

    ratio <- ceiling(nrow_data / nrow_position)
    position_final <- NULL

    for (i in 1:ratio) {
      position_final <- position_final %>%
        dplyr::bind_rows(position)
    }

    position_final <- position_final %>%
      dplyr::slice(1:nrow_data)

    result <- data %>%
      dplyr::mutate(`Sample Type` = "Unknown") %>%
      dplyr::select(.data$`Sample Type`, .data$`File Name`) %>%
      dplyr::mutate(`Sample ID` = 1) %>%
      dplyr::mutate(Path = data_path) %>%
      dplyr::mutate(`Instrument Method` = method_path) %>%
      dplyr::mutate(`Process Method` = NA) %>%
      dplyr::mutate(`Calibration File` = NA) %>%
      cbind(position_final) %>%
      dplyr::mutate(`Inj Vol` = inj_vol) %>%
      dplyr::mutate(`Level` = NA) %>%
      dplyr::mutate(`Sample Wt` = 0) %>%
      dplyr::mutate(`Sample Vol` = 0) %>%
      dplyr::mutate(`ISTD Amt` = 0) %>%
      dplyr::mutate(`Dil Factor` = 1) %>%
      dplyr::mutate(`L1 Study` = NA) %>%
      dplyr::mutate(`L2 Client` = NA) %>%
      dplyr::mutate(`L3 Laboratory` = NA) %>%
      dplyr::mutate(`L4 Company` = NA) %>%
      dplyr::mutate(`L5 Phone` = NA) %>%
      dplyr::mutate(`Comment` = NA) %>%
      dplyr::bind_cols(sample_name)

    if (!is.null(blank_every_n)) {
      n_blanks <- nrow(result) / blank_every_n

      blank_names <-
        tidyr::crossing(date, instrument, user, "blank", 1:n_blanks)

      blank <- tidyr::unite(blank_names, "File Name", sep = "_") %>%
        dplyr::mutate(`Sample Type` = "QC") %>%
        dplyr::select(.data$`Sample Type`, .data$`File Name`) %>%
        dplyr::mutate(`Sample ID` = 1) %>%
        dplyr::mutate(Path = data_path) %>%
        dplyr::mutate(`Instrument Method` = blank_method_path) %>%
        dplyr::mutate(`Process Method` = NA) %>%
        dplyr::mutate(`Calibration File` = NA) %>%
        dplyr::mutate(`Position` = blank_position) %>%
        dplyr::mutate(`Inj Vol` = blank_inj_vol) %>%
        dplyr::mutate(`Level` = NA) %>%
        dplyr::mutate(`Sample Wt` = 0) %>%
        dplyr::mutate(`Sample Vol` = 0) %>%
        dplyr::mutate(`ISTD Amt` = 0) %>%
        dplyr::mutate(`Dil Factor` = 1) %>%
        dplyr::mutate(`L1 Study` = NA) %>%
        dplyr::mutate(`L2 Client` = NA) %>%
        dplyr::mutate(`L3 Laboratory` = NA) %>%
        dplyr::mutate(`L4 Company` = NA) %>%
        dplyr::mutate(`L5 Phone` = NA) %>%
        dplyr::mutate(`Comment` = NA) %>%
        dplyr::mutate(`Sample Name` = "blank") %>%
        split(1:n_blanks)

      result <- result %>%
        split((as.numeric(rownames(result)) - 1) %/% 4) %>%
        purrr::map2_dfr(
          .y = blank,
          .f = ~ .y %>%
            dplyr::bind_rows(.x)
        )
    }
    if (export == FALSE) {
      return(result)
    }
    if (export_to_queue == FALSE) {
      # a file named queue.csv is exported
      write("Bracket Type=4", file = "queue.csv")
      data.table::fwrite(result,
        "queue.csv",
        append = TRUE,
        col.names = TRUE
      )
    }
    if (export_to_queue == TRUE) {
      # load queue interactively if no path is provided in the queue_path argument
      if (is.null(queue_path)) {
        queue_path <- file.choose(".")
      }
      queue <- data.table::fread(queue_path, skip = 1)

      queue <- queue %>%
        dplyr::bind_rows(result)

      write("Bracket Type=4", file = queue_path)
      data.table::fwrite(queue,
        file = queue_path,
        append = TRUE,
        col.names = TRUE
      )
    }
  }
