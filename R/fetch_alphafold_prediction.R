#' Fetch AlphaFold prediction
#'
#' Fetches atom level data for AlphaFold predictions either for selected proteins or whole
#' organisms.
#'
#' @param uniprot_ids optional, a character vector of UniProt identifiers for which predictions
#' should be fetched. This argument is mutually exclusive to the \code{organism_name} argument.
#' @param organism_name optional, a character value providing the name of an organism for which
#' all available AlphaFold predictions should be retreived. The name should be the capitalised
#' scientific species name (e.g. "Homo sapiens"). **Note:** Some organisms contain a lot of
#' predictions which might take a considerable amount of time and memory to fetch. Therefore, you
#' should be sure that your system can handle fetching predictions for these organisms. This
#' argument is mutually exclusive to the \code{uniprot_ids} argument.
#' @param timeout a numeric value specifying the time in seconds until the download of an organism
#' archive times out. The default is 3600 seconds.
#' @param return_data_frame a logical value that specifies if true, a data frame instead of a list
#' is returned. It is recommended to only use this if not many pdb structures are retrieved.
#' Default is FALSE.
#' @param show_progress a logical value that specifies if true, a progress bar will be shown.
#' Default is TRUE.
#'
#' @return A list that contains atom level data for AlphaFold predictions. If return_data_frame is
#' TRUE, a data frame with this information is returned instead. The data frame contains the
#' following columns:
#' \itemize{
#' \item{label_id: }{Uniquely identifies every atom in the prediction following the standardised
#' convention for mmCIF files.}
#' \item{type_symbol: }{The code used to identify the atom species representing this atom type.
#' This code is the element symbol.}
#' \item{label_atom_id: }{Uniquely identifies every atom for the given residue following the
#' standardised convention for mmCIF files.}
#' \item{label_comp_id: }{A chemical identifier for the residue. This is the three- letter code
#' for the amino acid.}
#' \item{label_asym_id: }{Chain identifier following the standardised convention for mmCIF files.
#' Since every prediction only contains one protein this is always "A".}
#' \item{label_seq_id: }{Uniquely and sequentially identifies residues for each protein. The
#' numbering corresponds to the UniProt amino acid positions.}
#' \item{x: }{The x coordinate of the atom.}
#' \item{y: }{The y coordinate of the atom.}
#' \item{z: }{The z coordinate of the atom.}
#' \item{prediction_score: }{Contains the prediction score for each residue.}
#' \item{auth_seq_id: }{Same as \code{label_seq_id}. But of type character.}
#' \item{auth_comp_id: }{Same as \code{label_comp_id}.}
#' \item{auth_asym_id: }{Same as \code{label_asym_id}.}
#' \item{uniprot_id: }{The UniProt identifier of the predicted protein.}
#' \item{score_quality: }{Score annotations.}
#' }
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom tibble tibble
#' @importFrom utils download.file untar
#' @importFrom readr read_tsv
#' @importFrom stringr str_replace_all str_detect
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' alphafold <- fetch_alphafold_prediction(
#'   uniprot_ids = c("F4HVG8", "O15552"),
#'   return_data_frame = TRUE
#' )
#'
#' head(alphafold, n = 10)
#' }
fetch_alphafold_prediction <- function(uniprot_ids = NULL,
                                       organism_name = NULL,
                                       timeout = 3600,
                                       return_data_frame = FALSE,
                                       show_progress = TRUE) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  if (!missing(uniprot_ids) & !missing(organism_name)) {
    stop(strwrap("Please only provide either a list of UniProt identifiers or one organism name!",
      prefix = "\n", initial = ""
    ))
  }
  # if organism name is provided fetch all information about that organism
  if (!missing(organism_name)) {
    organism_name <- match.arg(organism_name, c(
      "Arabidopsis thaliana",
      "Caenorhabditis elegans",
      "Candida albicans",
      "Danio rerio",
      "Dictyostelium discoideum",
      "Drosophila melanogaster",
      "Escherichia coli",
      "Glycine max",
      "Homo sapiens",
      "Leishmania infantum",
      "Methanocaldococcus jannaschii",
      "Mus musculus",
      "Mycobacterium tuberculosis",
      "Oryza sativa",
      "Plasmodium falciparum",
      "Rattus norvegicus",
      "Saccharomyces cerevisiae",
      "Schizosaccharomyces pombe",
      "Staphylococcus aureus",
      "Trypanosoma cruzi",
      "Zea mays"
    ))

    organism_file <- switch(organism_name,
      "Arabidopsis thaliana" = "UP000006548_3702_ARATH.tar",
      "Caenorhabditis elegans" = "UP000001940_6239_CAEEL.tar",
      "Candida albicans" = "UP000000559_237561_CANAL.tar",
      "Danio rerio" = "UP000000437_7955_DANRE.tar",
      "Dictyostelium discoideum" = "UP000002195_44689_DICDI.tar",
      "Drosophila melanogaster" = "UP000000803_7227_DROME.tar",
      "Escherichia coli" = "UP000000625_83333_ECOLI.tar",
      "Glycine max" = "UP000008827_3847_SOYBN.tar",
      "Homo sapiens" = "UP000005640_9606_HUMAN.tar",
      "Leishmania infantum" = "UP000008153_5671_LEIIN.tar",
      "Methanocaldococcus jannaschii" = "UP000000805_243232_METJA.tar",
      "Mus musculus" = "UP000000589_10090_MOUSE.tar",
      "Mycobacterium tuberculosis" = "UP000001584_83332_MYCTU.tar",
      "Oryza sativa" = "UP000059680_39947_ORYSJ.tar",
      "Plasmodium falciparum" = "UP000001450_36329_PLAF7.tar",
      "Rattus norvegicus" = "UP000002494_10116_RAT.tar",
      "Saccharomyces cerevisiae" = "UP000002311_559292_YEAST.tar",
      "Schizosaccharomyces pombe" = "UP000002485_284812_SCHPO.tar",
      "Staphylococcus aureus" = "UP000008816_93061_STAA8.tar",
      "Trypanosoma cruzi" = "UP000002296_353153_TRYCC.tar",
      "Zea mays" = "UP000007305_4577_MAIZE.tar"
    )

    url <- paste0("https://ftp.ebi.ac.uk/pub/databases/alphafold/", organism_file)

    # set new longer timeout and reset to standard once function exits.
    old <- options(timeout = timeout)
    on.exit(options(old))

    utils::download.file(url, destfile = paste0(tempdir(), "/alphafold.tar"))

    utils::untar(
      tarfile = paste0(tempdir(), "/alphafold.tar"),
      exdir = paste0(tempdir(), "/alphafold")
    )

    all_files <- paste0(
      tempdir(),
      "/alphafold/",
      list.files(
        path = paste0(
          tempdir(),
          "/alphafold"
        ),
        pattern = ".cif.gz"
      )
    )

    all_protein_ids <- str_extract(all_files,
      pattern = "(?<=AF-).+(?=-F1)"
    )

    names(all_files) <- all_protein_ids

    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(
        total = length(all_files),
        format = paste0(
          "Importing AlphaFold predictions for ",
          organism_name,
          "[:bar] :current/:total (:percent) :eta"
        )
      )
    }

    query_result <- purrr::map(
      .x = all_files,
      .f = ~ {
        if (show_progress == TRUE) {
          pb$tick()
        }
        readr::read_tsv(.x, col_names = FALSE, quote = "", show_col_types = FALSE, progress = FALSE) %>%
          dplyr::filter(stringr::str_detect(X1, pattern = "^ATOM\\s+\\d|^HETATM\\s+\\d")) %>%
          dplyr::mutate(X2 = stringr::str_replace_all(X1, pattern = "\\s+", replacement = " ")) %>%
          tidyr::separate(X2,
            sep = " ",
            into = c(
              "x1",
              "label_id",
              "type_symbol",
              "label_atom_id",
              "x2",
              "label_comp_id",
              "label_asym_id",
              "entity_id",
              "label_seq_id",
              "x3",
              "x",
              "y",
              "z",
              "site_occupancy",
              "prediction_score",
              "formal_charge",
              "auth_seq_id",
              "auth_comp_id",
              "auth_asym_id",
              "x4",
              "pdb_model_number",
              "uniprot_id",
              "x5",
              "x6",
              "x7"
            )
          ) %>%
          dplyr::select(-c(
            .data$X1,
            .data$x1,
            .data$x2,
            .data$x3,
            .data$x4,
            .data$x5,
            .data$x6,
            .data$x7,
            .data$formal_charge,
            .data$site_occupancy,
            .data$entity_id,
            .data$pdb_model_number
          )) %>%
          dplyr::mutate(
            label_id = as.numeric(.data$label_id),
            label_seq_id = as.numeric(.data$label_seq_id),
            x = as.numeric(.data$x),
            y = as.numeric(.data$y),
            z = as.numeric(.data$z),
            prediction_score = as.numeric(.data$prediction_score),
            auth_seq_id = .data$auth_seq_id
          ) %>%
          dplyr::mutate(score_quality = dplyr::case_when(
            .data$prediction_score > 90 ~ "very_good",
            .data$prediction_score > 70 ~ "confident",
            .data$prediction_score > 50 ~ "low",
            .data$prediction_score <= 50 ~ "very_low"
          ))
      }
    )
    # delete files after everything is done
    unlink(paste0(tempdir(), "/alphafold"), recursive = TRUE)
    unlink(paste0(tempdir(), "/alphafold.tar"))
  } else {
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
    }
    # remove NAs from UniProt IDs
    uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]

    batches <- purrr::map(
      .x = uniprot_ids,
      .f = ~ paste0("https://alphafold.ebi.ac.uk/files/AF-", .x, "-F1-model_v1.cif")
    )

    names(batches) <- uniprot_ids

    if (show_progress == TRUE) {
      pb <- progress::progress_bar$new(
        total = length(batches),
        format = "  Fetching AlphaFold predictions [:bar] :current/:total (:percent) :eta"
      )
    }

    query_result <- purrr::map(
      .x = batches,
      .f = ~ {
        # query information from database
          query <- try_query(.x,
            type = "text/tab-separated-values",
            col_names = FALSE,
            quote = "",
            show_col_types = FALSE,
            progress = FALSE
          )
          
          if (show_progress == TRUE) {
          pb$tick()
          }
          
        # only proceed with data if it was correctly retrieved
        if ("tbl" %in% class(query)) {
          query %>%
            dplyr::filter(stringr::str_detect(X1, pattern = "^ATOM\\s+\\d|^HETATM\\s+\\d")) %>%
            dplyr::mutate(X2 = stringr::str_replace_all(X1, pattern = "\\s+", replacement = " ")) %>%
            tidyr::separate(X2,
              sep = " ",
              into = c(
                "x1",
                "label_id",
                "type_symbol",
                "label_atom_id",
                "x2",
                "label_comp_id",
                "label_asym_id",
                "entity_id",
                "label_seq_id",
                "x3",
                "x",
                "y",
                "z",
                "site_occupancy",
                "prediction_score",
                "formal_charge",
                "auth_seq_id",
                "auth_comp_id",
                "auth_asym_id",
                "x4",
                "pdb_model_number",
                "uniprot_id",
                "x5",
                "x6",
                "x7"
              )
            ) %>%
            dplyr::select(-c(
              .data$X1,
              .data$x1,
              .data$x2,
              .data$x3,
              .data$x4,
              .data$x5,
              .data$x6,
              .data$x7,
              .data$formal_charge,
              .data$site_occupancy,
              .data$entity_id,
              .data$pdb_model_number
            )) %>%
            dplyr::mutate(
              label_id = as.numeric(.data$label_id),
              label_seq_id = as.numeric(.data$label_seq_id),
              x = as.numeric(.data$x),
              y = as.numeric(.data$y),
              z = as.numeric(.data$z),
              prediction_score = as.numeric(.data$prediction_score),
              auth_seq_id = .data$auth_seq_id
            ) %>%
            dplyr::mutate(score_quality = dplyr::case_when(
              .data$prediction_score > 90 ~ "very_good",
              .data$prediction_score > 70 ~ "confident",
              .data$prediction_score > 50 ~ "low",
              .data$prediction_score <= 50 ~ "very_low"
            ))
        } else {
          query
        }
      }
    )
  }
  
  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>% 
    purrr::keep(.p = ~ is.character(.x))
  
  error_table <- tibble::tibble(id = names(error_list),
                                error = unlist(error_list)) %>% 
    dplyr::distinct()
  
  if(nrow(error_table) != 0){
    message("The following IDs have not be retrieved correctly.")
    message(paste0(capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output
  
  query_result <- query_result %>% 
    purrr::keep(.p = ~ !is.character(.x))
  
  if (return_data_frame == FALSE) {
    return(query_result)
  } else {
    query_result_df <- purrr::map_dfr(
      .x = query_result,
      .f = ~.x
    )
    return(query_result_df)
  }
}
