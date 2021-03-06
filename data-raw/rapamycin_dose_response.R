# library(tidyverse)
# library(protti)
# 
# set.seed(123)
# 
# # Source: Piazza, I., Beaton, N., Bruderer, R. et al. A machine learning-based chemoproteomic approach to identify drug targets and binding sites in complex proteomes. Nat Commun 11, 4200 (2020). https://doi.org/10.1038/s41467-020-18071-x
# 
rapa <- read_protti("rapamycin_dose_response.csv")

# Filter to only contain necessary columns. Simplify file names. Annotate conditions with concentrations in pM.

rapa_filtered <- rapa %>%
  distinct(r_file_name, r_condition, eg_precursor_id, pg_protein_accessions, fg_quantity, pep_is_proteotypic, eg_is_decoy) %>%
  mutate(r_file_name = paste0("sample_", str_sub(r_file_name, start = 35, end = 36))) %>%
  mutate(r_condition = case_when(
    r_condition == 0 ~ 0,
    r_condition == 1 ~ 10,
    r_condition == 2 ~ 100,
    r_condition == 3 ~ 1000,
    r_condition == 4 ~ 10000,
    r_condition == 5 ~ 100000,
    r_condition == 6 ~ 1000000,
    r_condition == 7 ~ 10000000,
    r_condition == 8 ~ 100000000,
  ))

all_proteins <- unique(rapa_filtered$pg_protein_accessions)

all_proteins_wo_FKBP1A <- all_proteins[all_proteins != "P62942"]

sampled_bg <- sample(all_proteins_wo_FKBP1A, size = 39)

rapamycin_dose_response <- rapa_filtered %>%
  filter(pg_protein_accessions %in% c(sampled_bg, "P62942"))

usethis::use_data(rapamycin_dose_response, overwrite = TRUE)
