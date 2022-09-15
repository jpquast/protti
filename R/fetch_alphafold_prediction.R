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
#' is returned. It is recommended to only use this if information for few proteins is retrieved.
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
#' @importFrom utils capture.output
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
      "Zea mays",
      "Ajellomyces capsulatus",
      "Brugia malayi",
      "Campylobacter jejuni",
      "Cladophialophora carrionii",
      "Dracunculus medinensis",
      "Enterococcus faecium",
      "Fonsecaea pedrosoi",
      "Haemophilus influenzae",
      "Helicobacter pylori",
      "Klebsiella pneumoniae",
      "Leishmania infantum",
      "Madurella mycetomatis",
      "Mycobacterium leprae",
      "Mycobacterium tuberculosis",
      "Mycobacterium ulcerans",
      "Neisseria gonorrhoeae",
      "Nocardia brasiliensis",
      "Onchocerca volvulus",
      "Paracoccidioides lutzii",
      "Plasmodium falciparum",
      "Pseudomonas aeruginosa",
      "Salmonella typhimurium",
      "Schistosoma mansoni",
      "Shigella dysenteriae",
      "Sporothrix schenckii",
      "Staphylococcus aureus",
      "Streptococcus pneumoniae",
      "Strongyloides stercoralis",
      "Trichuris trichiura",
      "Trypanosoma brucei",
      "Trypanosoma cruzi",
      "Wuchereria bancrofti"
    ))

    organism_file <- switch(organism_name,
      "Arabidopsis thaliana" = "UP000006548_3702_ARATH_v3.tar",
      "Caenorhabditis elegans" = "UP000001940_6239_CAEEL_v3.tar",
      "Candida albicans" = "UP000000559_237561_CANAL_v3.tar",
      "Danio rerio" = "UP000000437_7955_DANRE_v3.tar",
      "Dictyostelium discoideum" = "UP000002195_44689_DICDI_v3.tar",
      "Drosophila melanogaster" = "UP000000803_7227_DROME_v3.tar",
      "Escherichia coli" = "UP000000625_83333_ECOLI_v3.tar",
      "Glycine max" = "UP000008827_3847_SOYBN_v3.tar",
      "Homo sapiens" = "UP000005640_9606_HUMAN_v3.tar",
      "Leishmania infantum" = "UP000008153_5671_LEIIN_v3.tar",
      "Methanocaldococcus jannaschii" = "UP000000805_243232_METJA_v3.tar",
      "Mus musculus" = "UP000000589_10090_MOUSE_v3.tar",
      "Mycobacterium tuberculosis" = "UP000001584_83332_MYCTU_v3.tar",
      "Oryza sativa" = "UP000059680_39947_ORYSJ_v3.tar",
      "Plasmodium falciparum" = "UP000001450_36329_PLAF7_v3.tar",
      "Rattus norvegicus" = "UP000002494_10116_RAT_v3.tar",
      "Saccharomyces cerevisiae" = "UP000002311_559292_YEAST_v3.tar",
      "Schizosaccharomyces pombe" = "UP000002485_284812_SCHPO_v3.tar",
      "Staphylococcus aureus" = "UP000008816_93061_STAA8_v3.tar",
      "Trypanosoma cruzi" = "UP000002296_353153_TRYCC_v3.tar",
      "Zea mays" = "UP000007305_4577_MAIZE_v3.tar",
      "Ajellomyces capsulatus" = "UP000001631_447093_AJECG_v3.tar",
      "Brugia malayi" = "UP000006672_6279_BRUMA_v3.tar",
      "Campylobacter jejuni" = "UP000000799_192222_CAMJE_v3.tar",
      "Cladophialophora carrionii" = "UP000094526_86049_9EURO1_v3.tar",
      "Dracunculus medinensis" = "UP000274756_318479_DRAME_v3.tar",
      "Enterococcus faecium" = "UP000325664_1352_ENTFC_v3.tar",
      "Fonsecaea pedrosoi" = "UP000053029_1442368_9EURO2_v3.tar",
      "Haemophilus influenzae" = "UP000000579_71421_HAEIN_v3.tar",
      "Helicobacter pylori" = "UP000000429_85962_HELPY_v3.tar",
      "Klebsiella pneumoniae" = "UP000007841_1125630_KLEPH_v3.tar",
      "Leishmania infantum" = "UP000008153_5671_LEIIN_v3.tar",
      "Madurella mycetomatis" = "UP000078237_100816_9PEZI1_v3.tar",
      "Mycobacterium leprae" = "UP000000806_272631_MYCLE_v3.tar",
      "Mycobacterium tuberculosis" = "UP000001584_83332_MYCTU_v3.tar",
      "Mycobacterium ulcerans" = "UP000020681_1299332_MYCUL_v3.tar",
      "Neisseria gonorrhoeae" = "UP000000535_242231_NEIG1_v3.tar",
      "Nocardia brasiliensis" = "UP000006304_1133849_9NOCA1_v3.tar",
      "Onchocerca volvulus" = "UP000024404_6282_ONCVO_v3.tar",
      "Paracoccidioides lutzii" = "UP000002059_502779_PARBA_v3.tar",
      "Plasmodium falciparum" = "UP000001450_36329_PLAF7_v3.tar",
      "Pseudomonas aeruginosa" = "UP000002438_208964_PSEAE_v3.tar",
      "Salmonella typhimurium" = "UP000001014_99287_SALTY_v3.tar",
      "Schistosoma mansoni" = "UP000008854_6183_SCHMA_v3.tar",
      "Shigella dysenteriae" = "UP000002716_300267_SHIDS_v3.tar",
      "Sporothrix schenckii" = "UP000018087_1391915_SPOS1_v3.tar",
      "Staphylococcus aureus" = "UP000008816_93061_STAA8_v3.tar",
      "Streptococcus pneumoniae" = "UP000000586_171101_STRR6_v3.tar",
      "Strongyloides stercoralis" = "UP000035681_6248_STRER_v3.tar",
      "Trichuris trichiura" = "UP000030665_36087_TRITR_v3.tar",
      "Trypanosoma brucei" = "UP000008524_185431_TRYB2_v3.tar",
      "Trypanosoma cruzi" = "UP000002296_353153_TRYCC_v3.tar",
      "Wuchereria bancrofti" = "UP000270924_6293_WUCBA_v3.tar"
    )

    url <- paste0("https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/", organism_file)

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
    # remove NAs from UniProt IDs
    uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]

    batches <- purrr::map(
      .x = uniprot_ids,
      .f = ~ paste0("https://alphafold.ebi.ac.uk/files/AF-", .x, "-F1-model_v3.cif")
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

  if (length(error_list) != 0) {
    error_table <- tibble::tibble(
      id = names(error_list),
      error = unlist(error_list)
    ) %>%
      dplyr::distinct()

    message("The following IDs have not been retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output

  query_result <- query_result %>%
    purrr::keep(.p = ~ !is.character(.x))

  if (length(query_result) == 0) {
    message("No valid information could be retrieved!")
    return(invisible(NULL))
  }

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
