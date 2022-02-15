#' Fetch structural information about protein-metal binding from MetalPDB
#'
#' Fetches information about protein-metal binding sites from the
#' MetalPDB database. A complete list of different possible search
#' queries can be found on their website.
#'
#' @param id_type a character value that specifies the type of the IDs provided to \code{id_value}.
#' Default is "uniprot". Possible options include: "uniprot", "pdb", "ec_number", "molecule" and
#' "organism".
#' @param id_value a character vector supplying IDs that are of the ID type that was specified in
#' \code{id_type}. E.g. UniProt IDs. Information for these IDs will be retreived.
#' @param site_type optional, a character value that specifies a nuclearity for which information
#' should be retrieved. The specific nuclearity can be supplied as e.g. "tetranuclear".
#' @param pfam optional, a character value that specifies a Pfam domain for which information
#' should be retrieved. The domain can be specified as e.g. "Carb_anhydrase".
#' @param cath optional, a character value that specifies a CATH ID for which information
#' should be retrieved. The ID can be specified as e.g. "3.10.200.10".
#' @param scop optional, a character value that specifies a SCOP ID for which information
#' should be retrieved. The ID can be specified as e.g. "b.74.1.1".
#' @param representative optional, a logical that indicates if only information of representative
#' sites of a family should be retrieved it can be specified here. A representative site is a
#' site selected to represent a cluster of equivalent sites. The selection is done by choosing
#' the PDB structure with the best X-ray resolution among those containing the sites in the
#' cluster. NMR structures are generally discarded in favor of X-ray structures, unless all the
#' sites in the cluster are found in NMR structures. If it is \code{TRUE}, only representative
#' sites are retrieved, if it is \code{FALSE}, all sites are retrieved.
#' @param metal optional, a character value that specifies a metal for which information
#' should be retrieved. The metal can be specified as e.g. "Zn".
#' @param ligands optional, a character value that specifies a metal ligand residue for which
#' information should be retrieved. The ligand can be specified as e.g. "His".
#' @param geometry optional, a character value that specifies a metal site geometry for which
#' information should be retrieved. The geometry can be specified here based on the three letter
#' code for geometries provided on their website.
#' @param coordination optional, a character value that specifies a coordination number for which
#' information should be retrieved. The number can be specified as e.g. "3".
#' @param donors optional, a character value that specifies a metal ligand atom for which
#' information should be retrieved. The atom can be specified as e.g. "S" for sulfur.
#' @param columns optional, a character vector that specifies specific columns that should be
#' retrieved based on the MetalPDB website. If
#' nothing is supplied here, all possible columns will be retrieved.
#' @param show_progress logical, if true, a progress bar will be shown. Default is TRUE.
#'
#' @return A data frame that contains information about protein-metal binding sites. The data
#' frame contains some columns that might not be self explanatory.
#' \itemize{
#' \item{auth_id_metal: }{Unique structure atom identifier of the metal, which is provided by
#' the author of the structure in order to match the identification used in the publication
#' that describes the structure.}
#' \item{auth_seq_id_metal: }{Residue identifier of the metal, which is provided by the author of
#' the structure in order to match the identification used in the publication that describes the
#' structure.}
#' \item{pattern: }{Metal pattern for each metal bound by the structure.}
#' \item{is_representative: }{A representative site is a site selected to represent a cluster of
#' equivalent sites. The selection is done by choosing the PDB structure with the best X-ray
#' resolution among those containing the sites in the cluster. NMR structures are generally
#' discarded in favor of X-ray structures, unless all the sites in the cluster are found in NMR
#' structures.}
#' \item{auth_asym_id_ligand: }{Chain identifier of the metal-coordinating ligand residues, which
#' is provided by the author of the structure in order to match the identification used in the
#' publication that describes the structure.}
#' \item{auth_seq_id_ligand: }{Residue identifier of the metal-coordinating ligand residues, which
#' is provided by the author of the structure in order to match the identification used in the
#' publication that describes the structure.}
#' \item{auth_id_ligand: }{Unique structure atom identifier of the metal-coordinating ligand r
#' esidues, which is provided by the author of the structure in order to match the identification
#' used in the publication that describes the structure.}
#' \item{auth_atom_id_ligand: }{Unique residue specific atom identifier of the metal-coordinating
#' ligand residues, which is provided by the author of the structure in order to match the
#' identification used in the publication that describes the structure.}
#' }
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom stringr str_detect
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_metal_pdb(id_value = c("P42345", "P00918")))
#'
#' fetch_metal_pdb(id_type = "pdb", id_value = c("1g54"), metal = "Zn")
#' }
#'
fetch_metal_pdb <- function(id_type = "uniprot",
                            id_value,
                            site_type = NULL,
                            pfam = NULL,
                            cath = NULL,
                            scop = NULL,
                            representative = NULL,
                            metal = NULL,
                            ligands = NULL,
                            geometry = NULL,
                            coordination = NULL,
                            donors = NULL,
                            columns = NULL,
                            show_progress = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package \"httr\" is needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }

  . <- NULL
  # Make sure arguments are provided correctly
  if (!missing(site_type)) {
    if (length(site_type) > 1) {
      stop('Please only provide one value for the "site_type" argument')
    }
    site_type <- paste0(",site_type:", site_type)
  }

  if (!missing(pfam)) {
    if (length(pfam) > 1) {
      stop('Please only provide one value for the "pfam" argument')
    }
    pfam <- paste0(",pfam:", pfam)
  }

  if (!missing(cath)) {
    if (length(cath) > 1) {
      stop('Please only provide one value for the "cath" argument')
    }
    cath <- paste0(",cath:", cath)
  }

  if (!missing(scop)) {
    if (length(scop) > 1) {
      stop('Please only provide one value for the "scop" argument')
    }
    scop <- paste0(",scop:", scop)
  }

  if (!missing(representative)) {
    if (representative != TRUE | representative != FALSE) {
      stop('Please only provide TRUE or FALSE for the "representative" argument')
    }
    if (representative == TRUE) {
      representative <- paste0(",representative:", representative)
    }
  }

  if (!missing(metal)) {
    if (length(metal) > 1) {
      stop('Please only provide one value for the "metal" argument')
    }
    metal <- paste0(",metal:", metal)
  }

  if (!missing(ligands)) {
    if (length(ligands) > 1) {
      stop('Please only provide one value for the "ligands" argument')
    }
    ligands <- paste0(",ligands:", ligands)
  }

  if (!missing(geometry)) {
    if (length(geometry) > 1) {
      stop('Please only provide one value for the "geometry" argument')
    }
    geometry <- paste0(",geometry:", geometry)
  }

  if (!missing(coordination)) {
    if (length(coordination) > 1) {
      stop('Please only provide one value for the "coordination" argument')
    }
    coordination <- paste0(",coordination:", coordination)
  }

  if (!missing(donors)) {
    if (length(donors) > 1) {
      stop('Please only provide one value for the "donors" argument')
    }
    donors <- paste0(",donors:", donors)
  }

  if (!missing(columns)) {
    columns <- paste0("&columns=", paste0(columns, collapse = ","))
  }

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(total = length(id_value), show_after = 0)
    pb$tick(0)
  }

  url <- "http://metalpdb.cerm.unifi.it/api?query="
  names(id_value) <- id_value

  query_result <- purrr::map(
    .x = id_value,
    .f = function(x) {
      query_url <- paste0(
        url,
        id_type, ":", x,
        site_type,
        pfam,
        cath,
        scop,
        representative,
        metal,
        ligands,
        geometry,
        coordination,
        donors,
        columns
      )
      query <- try_query(query_url, type = "application/json", simplifyDataFrame = TRUE)

      if (show_progress == TRUE) {
        pb$tick()
      }
      query
    }
  )

  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>%
    purrr::map(.f = ~ {
      if (length(.x) == 0) {
        paste("Not found in database")
      } else {
        .x
      }
    }) %>%
    purrr::keep(.p = ~ is.character(.x))

  error_table <- tibble::tibble(
    id = names(error_list),
    error = unlist(error_list)
  ) %>%
    dplyr::distinct()

  if (nrow(error_table) == length(id_value) && str_detect(error_table$error, pattern = "502")) {
    message("502 Bad Gateway. The server cannot be reached right now. Try again later!")
    return(invisible(NULL))
  }

  if (nrow(error_table) != 0) {
    message("The following IDs have not be retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output

  content <- query_result %>%
    purrr::map(.f = ~ {
      if (length(.x) == 0) {
        paste("Not found in database")
      } else {
        .x
      }
    }) %>%
    purrr::keep(.p = ~ !is.character(.x)) %>%
    purrr::map_dfr(
      .f = ~.x
    )

  if (nrow(content) == 0) {
    return(content)
  }

  # make sure that the data is complete even if there are columns missing
  columns <- c(
    "site",
    "organism",
    "scop",
    "site_type",
    "ec_number",
    "pfam",
    "metals",
    "molecule",
    "cath",
    "uniprot",
    "is_representative",
    "pdb"
  )
  should_be_here <- columns[!columns %in% colnames(content)]

  content_metal <- content %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here),
        nrow = nrow(content)
      )),
      should_be_here
    )) %>%
    tidyr::unnest(.data$metals)

  if ("metals" %in% colnames(content_metal)) {
    content_metal <- content_metal %>%
      dplyr::select(-c(.data$metals))
  }

  columns_metal <- c(
    "residue_pdb_number",
    "atom_pdb_number",
    "symbol",
    "ligands",
    "name",
    "pattern",
    "geometry",
    "coordination"
  )
  should_be_here_metal <- columns_metal[!columns_metal %in% colnames(content_metal)]

  content_ligand <- content_metal %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here_metal),
        nrow = nrow(content_metal)
      )),
      should_be_here_metal
    )) %>%
    dplyr::rename(
      auth_seq_id_metal = .data$residue_pdb_number,
      auth_id_metal = .data$atom_pdb_number,
      symbol_metal = .data$symbol
    ) %>%
    dplyr::group_by(.data$site) %>%
    dplyr::mutate(check = length(.data$ligands[[1]])) %>%
    dplyr::mutate(ligands = ifelse(.data$check == 0, NA, .data$ligands)) %>%
    tidyr::unnest_longer(.data$ligands) %>%
    dplyr::bind_cols(donors = .$ligands)

  columns_ligand <- c("check", "chain", "donors", "residue_pdb_number", "residue")
  should_be_here_ligand <- columns_ligand[!columns_ligand %in% colnames(content_ligand)]

  content_donor <- content_ligand %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here_ligand),
        nrow = nrow(content_ligand)
      )),
      should_be_here_ligand
    )) %>%
    tidyr::unnest_longer(.data$donors) %>%
    dplyr::bind_cols(atom = .$donors) %>%
    dplyr::select(-c(.data$ligands, .data$donors, .data$check))

  columns_donor <- c("atom", "symbol", "atom_pdb_number", "distance")
  should_be_here_donor <- columns_donor[!columns_donor %in% colnames(content_donor)]

  result <- content_donor %>%
    dplyr::bind_cols(stats::setNames(
      data.frame(matrix(
        ncol = length(should_be_here_donor),
        nrow = nrow(content_donor)
      )),
      should_be_here_donor
    )) %>%
    dplyr::rename(
      auth_asym_id_ligand = .data$chain,
      auth_seq_id_ligand = .data$residue_pdb_number,
      auth_id_ligand = .data$atom_pdb_number,
      auth_atom_id_ligand = .data$atom
    ) %>%
    dplyr::ungroup()

  return(result)
}
