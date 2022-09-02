# This is version 3, which was created on 22/08/12

# This code specificially creates a list of ChEBI IDs that appear in UniProt and that 
# are related to metals.
# When updated also update the documentation!

library(tidyverse)
library(protti)

# First we retrieve all reviewed annotations from UniProt
# The relevant columns are: cc_cofactor, cc_catalytic_activity

url_chebi_uniprot <- utils::URLencode("http://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=tsv&fields=accession,cc_cofactor,cc_catalytic_activity,ft_binding")
input_chebi_uniprot <- protti:::try_query(url_chebi_uniprot, timeout = 1000, progress = FALSE, show_col_types = FALSE)
colnames(input_chebi_uniprot) <- janitor::make_clean_names(c("accession", "cc_cofactor", "cc_catalytic_activity", "ft_binding"))

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

chebi <- protti::fetch_chebi(stars = c(2, 3))

# There are ChEBI IDs in UniProt without a formula that are still metal related
# If there is a new version of this reference created, check if these are still valid or new have been added.
# Run the code below to check for them:
# without_formula <- chebi %>%
#   dplyr::filter(id %in% extracted_chebi_uniprot$chebi_id) %>%
#   dplyr::filter(is.na(formula)) %>%
#   dplyr::filter(.data$type_name == "STANDARD") %>%
#   dplyr::distinct(.data$id, .data$chebi_accession, .data$star, .data$definition, .data$name)

metal_chebi_ids_wo_formula <- c("25213" = "25213", # a metal cation
                                "30408" = "18248", # iron-sulfur cluster
                                "30413" = "18248", # heme
                                "38201" = "25107", # a bacteriochlorophyll
                                "60240" = "60240", # a divalent metal cation
                                "60242" = "60242", # a monovalent cation (this can also be non-metal so should be handled with care)
                                "60400" = "28112,18248", # [Ni-Fe-S] cluster
                                "190135" = "18248" # [2Fe-2S] cluster
                                )

# Create annotation vector that indicates what ChEBI ID a certain element symbol should be annotated with

metal_chebi_annotation <- setNames(metal_list$chebi_id, metal_list$symbol)

# Extract all ChEBI IDs that appear in UniProt and that contain a metal in their formula or that are in the metal_chebi_ids_wo_formula vector above

metal_chebi_uniprot <- chebi %>% 
  dplyr::filter(.data$type_name == "STANDARD") %>% 
  dplyr::filter(id %in% extracted_chebi_uniprot$chebi_id) %>% 
  dplyr::distinct(id, chebi_accession, definition, star, type_name, name, formula, charge, monoisotopic_mass) %>% 
  dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) | 
                  id %in% as.numeric(names(metal_chebi_ids_wo_formula))) %>% 
  dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>% 
  tidyr::unnest(extract_formula) %>% 
  dplyr::mutate(metal_atom_id = ifelse(is.na(extract_formula), 
                                       metal_chebi_ids_wo_formula[as.character(id)],
                                       metal_chebi_annotation[extract_formula])) %>% 
  dplyr::select(-.data$extract_formula) %>% 
  dplyr::group_by(.data$id) %>% 
  dplyr::mutate(metal_atom_id = paste0(.data$metal_atom_id, collapse = ",")) %>% 
  dplyr::distinct()

# In version 1 there were 147 metal related ChEBI IDs
# In version 2 there were 188 metal related ChEBI IDs
# In version 3 there were 188 metal related ChEBI IDs
length(unique(metal_chebi_uniprot$id))

usethis::use_data(metal_chebi_uniprot, overwrite = TRUE)
