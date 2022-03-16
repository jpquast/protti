#' Colour scheme for protti
#'
#' A colour scheme for protti that contains 100 colours.
#'
#' @format A vector containing 100 colours
#' @source Dina's imagination.
"protti_colours"

#' Viridis colour scheme
#'
#' A colour scheme by the viridis colour scheme from the viridis R package.
#'
#' @format A vector containing 256 colours
#' @source viridis R package
"viridis_colours"

#' Rapamycin 10 uM example data
#'
#' Rapamycin example data used for the vignette about binary control/treated data. The data was
#' obtained from \href{https://www.nature.com/articles/s41467-020-18071-x}{Piazza 2020}
#' and corresponds to experiment 18. FKBP1A the rapamycin binding protein and 49 other randomly
#' sampled proteins were used for this example dataset. Furthermore, only the DMSO control and the
#' 10 uM condition were used.
#'
#' @format A data frame containing peptide level data from a Spectronaut report.
#' @source Piazza, I., Beaton, N., Bruderer, R. et al. A machine learning-based chemoproteomic
#' approach to identify drug targets and binding sites in complex proteomes. Nat Commun 11, 4200
#' (2020). https://doi.org/10.1038/s41467-020-18071-x
"rapamycin_10uM"

#' Rapamycin dose response example data
#'
#' Rapamycin example data used for the vignette about dose response data. The data was obtained
#' from \href{https://www.nature.com/articles/s41467-020-18071-x}{Piazza 2020} and corresponds
#' to experiment 18. FKBP1A the rapamycin binding protein and 39 other randomly sampled proteins
#' were used for this example dataset. The concentration range includes the following points:
#' 0 (DMSO control), 10 pM, 100 pM, 1 nM, 10 nM, 100 nM, 1 uM, 10 uM and 100 uM.
#'
#' @format A data frame containing peptide level data from a Spectronaut report.
#' @source Piazza, I., Beaton, N., Bruderer, R. et al. A machine learning-based chemoproteomic
#' approach to identify drug targets and binding sites in complex proteomes. Nat Commun 11, 4200
#' (2020). https://doi.org/10.1038/s41467-020-18071-x
"rapamycin_dose_response"

#' Structural analysis example data
#'
#' Example data used for the vignette about structural analysis. The data was obtained from
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867420316913}{Cappelletti 2021}
#' and corresponds to two separate experiments. Both experiments were limited proteolyis coupled to
#' mass spectrometry (LiP-MS) experiments conducted on purified proteins. The first protein is
#' phosphoglycerate kinase 1 (pgk) and it was treated with 25mM 3-phosphoglyceric acid (3PG).
#' The second protein is phosphoenolpyruvate-protein phosphotransferase (ptsI) and it was treated
#' with 25mM fructose 1,6-bisphosphatase (FBP). From both experiments only peptides belonging to
#' either protein were used for this data set. The ptsI data set contains precursor level data
#' while the pgk data set contains peptide level data. The pgk data can be obtained from
#' supplementary table 3 from the tab named "pgk+3PG". The ptsI data is only included as raw data
#' and was analysed using the functions of this package.
#'
#' @format A data frame containing differential abundances and adjusted p-values for
#' peptides/precursors of two proteins.
#' @source Cappelletti V, Hauser T, Piazza I, Pepelnjak M, Malinovska L, Fuhrer T, Li Y, Dörig C,
#' Boersema P, Gillet L, Grossbach J, Dugourd A, Saez-Rodriguez J, Beyer A, Zamboni N, Caflisch A,
#' de Souza N, Picotti P. Dynamic 3D proteomes reveal protein functional alterations at high
#' resolution in situ. Cell. 2021 Jan 21;184(2):545-559.e22. doi: 10.1016/j.cell.2020.12.021.
#' Epub 2020 Dec 23. PMID: 33357446; PMCID: PMC7836100.
"ptsi_pgk"
