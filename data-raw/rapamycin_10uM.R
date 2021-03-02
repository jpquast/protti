# library(tidyverse)
# library(protti)
# 
# set.seed(1234)
#
# # Source: Piazza, I., Beaton, N., Bruderer, R. et al. A machine learning-based chemoproteomic approach to identify drug targets and binding sites in complex proteomes. Nat Commun 11, 4200 (2020). https://doi.org/10.1038/s41467-020-18071-x
# 
# rapa <- read_protti("rapamycin_dose_response.csv")
# 
# # filter to only retain DMSO control and 10 uM concentration
# 
# rapa_filtered <- rapa %>%
#   distinct(r_file_name, r_condition, pep_stripped_sequence, eg_precursor_id, pg_protein_accessions, fg_quantity, pep_is_proteotypic, eg_is_decoy) %>%
#   filter(r_condition == 0 | r_condition == 7) %>%
#   mutate(r_condition = ifelse(r_condition == 0, "control", "rapamycin")) %>%
#   mutate(r_file_name = paste0(r_condition, "_", str_sub(r_file_name, start = 35, end = 36)))
# 
# all_proteins <- unique(rapa_filter$pg_protein_accessions)
# 
# all_proteins_wo_FKBP1A <- all_proteins[all_proteins != "P62942"]
# 
# sampled_bg <- sample(all_proteins_wo_FKBP1A, size = 49)
# 
# rapamycin_10uM <- rapa_filtered %>%
#   filter(pg_protein_accessions %in% c(sampled_bg, "P62942"))
# 
# usethis::use_data(rapamycin_10uM, overwrite = TRUE)
