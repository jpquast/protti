# library(tidyverse)
# library(protti)
#
# # Source: Cappelletti V, Hauser T, Piazza I, Pepelnjak M, Malinovska L, Fuhrer T, Li Y, Dörig C, Boersema P, Gillet L, Grossbach J, Dugourd A, Saez-Rodriguez J, Beyer A, Zamboni N, Caflisch A, de Souza N, Picotti P. Dynamic 3D proteomes reveal protein functional alterations at high resolution in situ. Cell. 2021 Jan 21;184(2):545-559.e22. doi: 10.1016/j.cell.2020.12.021. Epub 2020 Dec 23. PMID: 33357446; PMCID: PMC7836100.
#
# # The pgk data set is from supplementary table 3, the tab is called "pgk+3PG". The data does not contain precursor level data since charge states are
# # missing from peptides.
# pgk <- read_protti("pgk.csv")
#
# # The ptsI data set is not part of the supplementary tables. The raw data is included in the pride repository. We exported the Spectronaut report
# # and analysed that data using prottis standard pipeline.
# ptsi <- read_protti("ptsi.csv")
#
# # pgk data tidying
#
# pgk_tidy <- pgk %>%
#   filter(concentration == "25mM") %>% # filter to only retain the 25 mM concentration
#   rename(eg_precursor_id = peptide_sequence,
#          pg_protein_accessions = uniprot_id,
#          diff = log2fc,
#          adj_pval = qvalue) %>%
#   distinct(eg_precursor_id,
#            diff,
#            adj_pval,
#            pg_protein_accessions) %>%
#   mutate(pep_stripped_sequence = str_remove_all(eg_precursor_id, pattern = "(?<=\\[)[\\w\\(\\)\\s\\-]+(?=\\])")) %>% # removes "[Carbamidomethyl]" from peptides.
#   mutate(pep_stripped_sequence = str_remove_all(pep_stripped_sequence, pattern = "[\\[\\]]"))
#
# # ptsi data tidying
#
# ptsi_tidy <- ptsi %>%
#   rename(eg_precursor_id = precursor_id)
#
# # combining data
#
# ptsi_pgk <- pgk_tidy %>%
#   bind_rows(ptsi_tidy)
#
# usethis::use_data(ptsi_pgk, overwrite = TRUE)
