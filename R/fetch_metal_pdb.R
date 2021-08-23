#' Fetch structural information about protein-metal binding from MetalPDB
#'
#' Fetches information about protein-metal binding sites from the \href{https://metalpdb.cerm.unifi.it}{MetalPDB} database.
#' A complete list of different search queries possible can be found \href{https://metalpdb.cerm.unifi.it/api_help}{here}.
#'
#' @param id_type character, specifying the type of the IDs provided to \code{id_value}. Default is "uniprot". Possible options include:
#' "uniprot", "pdb", "ec_number", "molecule" and "organism"
#' @param id_value a character vector supplying IDs that have the ID type that was specified in \code{id_type}. E.g. UniProt IDs.
#' Information for these IDs will be retreived.
#' @param site_type optional character, if only information about structures with a certain nuclearity should be retrieved, the
#' specific nuclearity can be supplied here. E.g. "tetranuclear".
#' @param pfam optional character, if only information about structures with a certain Pfam domain should be retrieved, the domain
#' can be specified here. E.g. "Carb_anhydrase".
#' @param cath optional character, if only information about structures with a certain CATH ID should be retreived, the ID can be
#' specified here. E.g. "3.10.200.10".
#' @param scop optional character, if only information about structures with a certain SCOP ID should be retreived, the ID can be
#' specified here. E.g. "b.74.1.1".
#' @param representative optional logical, if only information of representative sites of a family should be retrieved it can
#' be specified here. A representative site is a site selected to represent a cluster of equivalent sites. The
#' selection is done by choosing the PDB structure with the best X-ray resolution among those containing the sites in the cluster.
#' NMR structures are generally discarded in favor of X-ray structures, unless all the sites in the cluster are found in NMR structures.
#' If it is \code{TRUE}, only representative sites are retrieved, if it is \code{FALSE}, all sites are retrieved.
#' @param metal optional character, if only information about structures with a certain metal should be retrieved, the metal can be
#' specified here. E.g. "Zn".
#' @param ligands optional character, if only information about structures with a certain metal ligand residue should be retrieved,
#' the ligand can be specified here. E.g. "His".
#' @param geometry optional character, if only information about structures with a certain metal site geometry should be retrieved,
#' the geometry can be specified here based on the three letter code for geometries provided on \href{https://metalpdb.cerm.unifi.it/perGeometry}{MetalPDB}.
#' @param coordination optional character, if only information about structures with a certain coordination number should be retrieved,
#' the number can be specified here. E.g. "3".
#' @param donors optional character, if only information about structures with a certain metal ligand atom should be retrieved,
#' the atom can be specified here. E.g. "S" for sulfur.
#' @param columns optional character vector, if only specific columns should be retrieved these can be specified here based on the MetalPDB
#' \href{https://metalpdb.cerm.unifi.it/api_help}{website}. If nothing is supplied here, all possible columns will be retrieved.
#' @param show_progress logical, if true, a progress bar will be shown. Default is TRUE.
#'
#' @return A data frame that contains information about protein-metal binding sites. The data frame contains some columns
#' that might not be self explanatory.
#' \itemize{
#' \item{auth_id_metal: }{Unique structure atom identifier of the metal, which is provided by the author of the structure in order to match the identification
#' used in the publication that describes the structure.}
#' \item{auth_seq_id_metal: }{Residue identifier of the metal, which is provided by the author of the structure in order to match the identification
#' used in the publication that describes the structure.}
#' \item{pattern: }{Metal pattern for each metal bound by the structure.}
#' \item{is_representative: }{A representative site is a site selected to represent a cluster of equivalent sites. The selection is done
#' by choosing the PDB structure with the best X-ray resolution among those containing the sites in the cluster. NMR structures are generally
#' discarded in favor of X-ray structures, unless all the sites in the cluster are found in NMR structures.}
#' \item{auth_asym_id_ligand: }{Chain identifier of the metal-coordinating ligand residues, which is provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure.}
#' \item{auth_seq_id_ligand: }{Residue identifier of the metal-coordinating ligand residues, which is provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure.}
#' \item{auth_id_ligand: }{Unique structure atom identifier of the metal-coordinating ligand residues, which is provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure.}
#' \item{auth_atom_id_ligand: }{Unique residue specific atom identifier of the metal-coordinating ligand residues, which is provided by the author of the structure in order
#' to match the identification used in the publication that describes the structure.}
#' }
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
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

  content <- purrr::map_dfr(
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
      if (!is.null(id_value)) {
        query <- try_query(query_url, type = "application/json", simplifyDataFrame = TRUE)
      }

      if (show_progress == TRUE & class(query) %in% c("data.frame", "list")) {
        pb$tick()
      }
      # if previous ID had a connection problem change IDs to NULL, which breaks the mapping.
      if (!class(query) %in% c("data.frame", "list")) {
        id_value <<- NULL
      }
      query
    }
  )

  if (nrow(content) == 0) {
    return(content)
  }

  # make sure that the data is complete even if there are columns missing
  columns <- c("site", "organism", "scop", "site_type", "ec_number", "pfam", "metals", "molecule", "cath", "uniprot", "is_representative", "pdb")
  should_be_here <- columns[!columns %in% colnames(content)]

  content_metal <- content %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here), nrow = nrow(content))), should_be_here)) %>%
    tidyr::unnest(.data$metals)

  if ("metals" %in% colnames(content_metal)) {
    content_metal <- content_metal %>%
      dplyr::select(-c(.data$metals))
  }

  columns_metal <- c("residue_pdb_number", "atom_pdb_number", "symbol", "ligands", "name", "pattern", "geometry", "coordination")
  should_be_here_metal <- columns_metal[!columns_metal %in% colnames(content_metal)]

  content_ligand <- content_metal %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here_metal), nrow = nrow(content_metal))), should_be_here_metal)) %>%
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
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here_ligand), nrow = nrow(content_ligand))), should_be_here_ligand)) %>%
    tidyr::unnest_longer(.data$donors) %>%
    dplyr::bind_cols(atom = .$donors) %>%
    dplyr::select(-c(.data$ligands, .data$donors, .data$check))

  columns_donor <- c("atom", "symbol", "atom_pdb_number", "distance")
  should_be_here_donor <- columns_donor[!columns_donor %in% colnames(content_donor)]

  result <- content_donor %>%
    dplyr::bind_cols(stats::setNames(data.frame(matrix(ncol = length(should_be_here_donor), nrow = nrow(content_donor))), should_be_here_donor)) %>%
    dplyr::rename(
      auth_asym_id_ligand = .data$chain,
      auth_seq_id_ligand = .data$residue_pdb_number,
      auth_id_ligand = .data$atom_pdb_number,
      auth_atom_id_ligand = .data$atom
    ) %>%
    dplyr::ungroup()

  return(result)
}
