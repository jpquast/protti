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
#' @param version a character value that specifies the alphafold version that should be used. This
#' is regularly updated by the database. We always try to make the current version the default version.
#' Available version can be found here: https://ftp.ebi.ac.uk/pub/databases/alphafold/
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
                                       version = "v4",
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
      "Arabidopsis thaliana" = "UP000006548_3702_ARATH",
      "Caenorhabditis elegans" = "UP000001940_6239_CAEEL",
      "Candida albicans" = "UP000000559_237561_CANAL",
      "Danio rerio" = "UP000000437_7955_DANRE",
      "Dictyostelium discoideum" = "UP000002195_44689_DICDI",
      "Drosophila melanogaster" = "UP000000803_7227_DROME",
      "Escherichia coli" = "UP000000625_83333_ECOLI",
      "Glycine max" = "UP000008827_3847_SOYBN",
      "Homo sapiens" = "UP000005640_9606_HUMAN",
      "Leishmania infantum" = "UP000008153_5671_LEIIN",
      "Methanocaldococcus jannaschii" = "UP000000805_243232_METJA",
      "Mus musculus" = "UP000000589_10090_MOUSE",
      "Mycobacterium tuberculosis" = "UP000001584_83332_MYCTU",
      "Oryza sativa" = "UP000059680_39947_ORYSJ",
      "Plasmodium falciparum" = "UP000001450_36329_PLAF7",
      "Rattus norvegicus" = "UP000002494_10116_RAT",
      "Saccharomyces cerevisiae" = "UP000002311_559292_YEAST",
      "Schizosaccharomyces pombe" = "UP000002485_284812_SCHPO",
      "Staphylococcus aureus" = "UP000008816_93061_STAA8",
      "Trypanosoma cruzi" = "UP000002296_353153_TRYCC",
      "Zea mays" = "UP000007305_4577_MAIZE",
      "Ajellomyces capsulatus" = "UP000001631_447093_AJECG",
      "Brugia malayi" = "UP000006672_6279_BRUMA",
      "Campylobacter jejuni" = "UP000000799_192222_CAMJE",
      "Cladophialophora carrionii" = "UP000094526_86049_9EURO1",
      "Dracunculus medinensis" = "UP000274756_318479_DRAME",
      "Enterococcus faecium" = "UP000325664_1352_ENTFC",
      "Fonsecaea pedrosoi" = "UP000053029_1442368_9EURO2",
      "Haemophilus influenzae" = "UP000000579_71421_HAEIN",
      "Helicobacter pylori" = "UP000000429_85962_HELPY",
      "Klebsiella pneumoniae" = "UP000007841_1125630_KLEPH",
      "Leishmania infantum" = "UP000008153_5671_LEIIN",
      "Madurella mycetomatis" = "UP000078237_100816_9PEZI1",
      "Mycobacterium leprae" = "UP000000806_272631_MYCLE",
      "Mycobacterium tuberculosis" = "UP000001584_83332_MYCTU",
      "Mycobacterium ulcerans" = "UP000020681_1299332_MYCUL",
      "Neisseria gonorrhoeae" = "UP000000535_242231_NEIG1",
      "Nocardia brasiliensis" = "UP000006304_1133849_9NOCA1",
      "Onchocerca volvulus" = "UP000024404_6282_ONCVO",
      "Paracoccidioides lutzii" = "UP000002059_502779_PARBA",
      "Plasmodium falciparum" = "UP000001450_36329_PLAF7",
      "Pseudomonas aeruginosa" = "UP000002438_208964_PSEAE",
      "Salmonella typhimurium" = "UP000001014_99287_SALTY",
      "Schistosoma mansoni" = "UP000008854_6183_SCHMA",
      "Shigella dysenteriae" = "UP000002716_300267_SHIDS",
      "Sporothrix schenckii" = "UP000018087_1391915_SPOS1",
      "Staphylococcus aureus" = "UP000008816_93061_STAA8",
      "Streptococcus pneumoniae" = "UP000000586_171101_STRR6",
      "Strongyloides stercoralis" = "UP000035681_6248_STRER",
      "Trichuris trichiura" = "UP000030665_36087_TRITR",
      "Trypanosoma brucei" = "UP000008524_185431_TRYB2",
      "Trypanosoma cruzi" = "UP000002296_353153_TRYCC",
      "Wuchereria bancrofti" = "UP000270924_6293_WUCBA"
    )

    url <- paste0("https://ftp.ebi.ac.uk/pub/databases/alphafold/", version, "/", organism_file, "_", version, ".tar")

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
            "X1",
            "x1",
            "x2",
            "x3",
            "x4",
            "x5",
            "x6",
            "x7",
            "formal_charge",
            "site_occupancy",
            "entity_id",
            "pdb_model_number"
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
      .f = ~ paste0("https://alphafold.ebi.ac.uk/files/AF-", .x, "-F1-model_", version, ".cif")
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
