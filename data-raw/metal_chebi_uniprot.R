# This is version 1, which was created on 22/06/03

# This code specificially creates a list of ChEBI IDs that appear in UniProt and that 
# are related to metals.
# When updated also update the documentation!

library(tidyverse)
library(protti)

# First we retrieve all reviewed annotations from UniProt
# The relevant columns are: cc_cofactor, cc_catalytic_activity

url_chebi_uniprot <- utils::URLencode("http://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=tsv&fields=accession,cc_cofactor,cc_catalytic_activity")
input_chebi_uniprot <- protti:::try_query(url_chebi_uniprot, timeout = 1000, progress = FALSE, show_col_types = FALSE)
colnames(input_chebi_uniprot) <- janitor::make_clean_names(c("accession", "cc_cofactor", "cc_catalytic_activity"))

# pivot columns longer for extraction of ChEBI IDs

extracted_chebi_uniprot <- input_chebi_uniprot %>% 
  tidyr::pivot_longer(-accession,
                      names_to = "source",
                      values_to = "ids"
  ) %>%
  dplyr::filter(ids != "") %>% 
  # extract all chebi IDs
  dplyr::mutate(chebi_id = stringr::str_extract_all(
                               ids,
                               pattern = "(?<=CHEBI:)[:digit:]+")) %>% 
  tidyr::unnest(chebi_id) %>% 
  dplyr::distinct(accession, chebi_id) %>% 
  dplyr::mutate(chebi_id = as.numeric(chebi_id))

# Retrieve information from the ChEBI database

chebi <- protti::fetch_chebi()

# There are ChEBI IDs in UniProt without a formula but that are metal related
# If there is a new version of this reference created, check if these are still valid or new have been added.
# Run the code below to check for them:
# without_formula <- chebi %>% 
#   dplyr::filter(id %in% extracted_chebi_uniprot$chebi_id) %>% 
#   dplyr::filter(is.na(formula))

metal_chebi_ids_wo_formula <- c(25213, # a metal cation
                                30408, # iron-sulfur cluster
                                30413, # heme
                                38201, # a bacteriochlorophyll
                                60240, # a divalent metal cation
                                60242, # a monovalent cation (this can also be non-metal so should be handled with care)
                                60400 # [Ni-Fe-S] cluster
                                )

# Extract all ChEBI IDs that appear in UniProt and that contain a matal in their formula or that are in the metal_chebi_ids_wo_formula vector above

metal_chebi_uniprot <- chebi %>% 
  dplyr::filter(id %in% extracted_chebi_uniprot$chebi_id) %>% 
  dplyr::distinct(id, chebi_accession, definition, star, type_name, name, formula, mass, charge, monoisotopic_mass) %>% 
  dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) | 
                  id %in% metal_chebi_ids_wo_formula)

# In version 1 there were 147 metal related ChEBI IDs
length(unique(metal_chebi_uniprot$id))

usethis::use_data(metal_chebi_uniprot, overwrite = TRUE)
