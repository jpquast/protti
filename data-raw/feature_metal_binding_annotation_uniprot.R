# This is version 1, which was created on 22/06/03

# This code annotates all possible names that can appear in the
# ft_metal column in UniProt with ChEBI IDs that match.
# It also includes annotations for complexed metals that appear 
# in cc_cofactor and cc_catalytic_activity UniProt columns. In 
# addition we also we also have annotations for complexed metals
# from the QuickGO database.
# As complexed metals we define any metal of which the formula 
# is not simply the element symbol.
# When updated also update the documentation!

library(tidyverse)
library(protti)

# First we retrieve all reviewed annotations from UniProt
# The relevant columns are: ft_metal

url_fmb_uniprot <- utils::URLencode("http://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=tsv&fields=accession,ft_metal")
input_fmb_uniprot <- protti:::try_query(url_fmb_uniprot, timeout = 600, progress = FALSE, show_col_types = FALSE)
colnames(input_fmb_uniprot) <- janitor::make_clean_names(c("accession", "ft_metal"))

# Extract all possible names

extracted_fmb_uniprot <- input_fmb_uniprot %>%
  tidyr::drop_na(ft_metal) %>%
  dplyr::mutate(fmb_name = stringr::str_extract_all(
    ft_metal,
    pattern = "METAL.+?(?=METAL)|METAL.+$"
  )) %>%
  tidyr::unnest(fmb_name) %>%
  dplyr::mutate(fmb_name = str_extract(fmb_name, pattern = '(?<=/note=\\")[^\\";]+(?=[\\";])')) %>%
  mutate(fmb_name = stringr::str_trim(stringr::str_replace_all(fmb_name, pattern = "[:space:]\\d{1,2}(?=[:space:]|$)", replacement = ""))) %>%
  dplyr::distinct(fmb_name) %>%
  dplyr::arrange(fmb_name) %>%
  tidyr::drop_na()

# retrieve ChEBI relations in order to find all possible sub IDs

chebi_relations <- fetch_chebi(relation = TRUE)

# create data frame containing annotations
# I used the command bellow to generate a list of names in the console
# This list was copied and pasted below and information was added to it manually.
# I compared names to often occurring ChEBI IDs in protti::metal_chebi_uniprot.
# I generally used fetch_chebi() for finding matching ChEBI IDs as well.
# message(paste0('"', paste0(extracted_fmb_uniprot$fmb_name, collapse = '",\n"'), '"'))

# fmb_name column is pasted here
fmb_annotation_uniprot_manuall <- data.frame(fmb_name = c(
  paste0(c(
    "Arsenite",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "22633",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "22633"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Ca in calcium-manganese-oxide [Ca-4Mn-5O]",
    "39123",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "39123"
    )), collapse = ","),
    "calcium-manganese-oxide [Ca-4Mn-5O]",
    "",
    "There is no ChEBI ID yet refering to the calcium-manganese-oxide cluster."
  ), collapse = ";"),
  paste0(c(
    "Cadmium",
    "63063",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "63063"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Calcium",
    "39123",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "39123"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Chloride",
    "",
    "",
    "",
    "",
    "Chloride is not a metal and is therefore not annotated with a ChEBI ID"
  ), collapse = ";"),
  paste0(c(
    "Co(2+)",
    "48828",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48828"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Co(3+)",
    "49415",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49415"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt (adenosylcob(III)alamin axial ligand)",
    "49415",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49415"
    )), collapse = ","),
    "18408",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18408"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt (adenosylcobalamin axial ligand)",
    "49415",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49415"
    )), collapse = ","),
    "18408",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18408"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt (cob(II)alamin axial ligand)",
    "48828",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48828"
    )), collapse = ","),
    "16304",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "16304"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt (Cob(II)alamin axial ligand)",
    "48828",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48828"
    )), collapse = ","),
    "16304",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "16304"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Cobalt (methylcob(III)alamin axial ligand)",
    "49415",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49415"
    )), collapse = ","),
    "28115",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "28115"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper A",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "47357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "47357"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper A1",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "47357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "47357"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper A2",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "47357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "47357"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper B",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "47357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "47357"
    )), collapse = ","),
    "Copper B is a copper B center if there is no Copper A annotation for the protein. There is no ChEBI ID for that though."
  ), collapse = ";"),
  paste0(c(
    "Copper Z1",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "33730",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33730"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper Z2",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "33730",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33730"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper Z3",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "33730",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33730"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper Z4",
    "23378",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23378"
    )), collapse = ","),
    "33730",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33730"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper(1+)",
    "49552",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49552"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Copper(2+)",
    "29036",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29036"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Cu(+)",
    "49552",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49552"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Cu(2+)",
    "29036",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29036"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Divalent cations",
    "60240",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60240"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Divalent metal cation",
    "60240",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60240"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Fe(2+)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Fe(3+)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "FeMo cofactor",
    "24875|37239",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ",")),
    collapse = "|"),
    "30409",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30409"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-oxo-sulfur (4Fe-2O-2S)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "60519",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60519"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "30408",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30408"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur-carbon (8Fe-9S-C-homocitryl)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "60504",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60504"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur (2Fe-2S)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "49601",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49601"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur (3Fe-4S)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "21137",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "21137"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur (4Fe-4S-S-AdoMet)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "49883",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49883"
    )), collapse = ","),
    "Should be the same as a normal Iron-sulfur (4Fe-4S) cluster, but has an exchangable SAM."
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur (4Fe-4S)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "49883",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49883"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron-sulfur (8Fe-7S)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "21143",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "21143"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (Fe-coproporphyrin III axial ligand)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "68438",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "68438"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme A axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "24479",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24479"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme A3 axial ligand)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "83282",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83282"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "30413",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30413"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme B axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b distal ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b proximal ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b558 axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b562 axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b566 axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme b595 axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme bD axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme bP axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme c axial ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "61717",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61717"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme C axial ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "61717",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61717"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme d1 axial ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60549",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60549"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme d1 distal ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60549",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60549"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme d1 proximal ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60549",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60549"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme distal ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "30413",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30413"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme o axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "24480",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24480"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme O axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "24480",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24480"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (heme proximal ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "30413",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30413"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (high-spin heme A3 axial ligand)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "83282",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83282"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (high-spin heme b axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (low-spin heme A axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "24479",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24479"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (low-spin heme b axial ligand)",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "26355",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26355"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron (siroheme axial ligand)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60052",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60052"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Iron B",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Li(+)",
    "49713",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49713"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll a axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "61720",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61720"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "38201",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "38201"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll b axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "22686",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "22686"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll c axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "60197",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60197"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll d axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90955",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90955"
    )), collapse = ","),
    "This should be bacteriochlorophyll d, which does not exist in ChEBI. 
    The only protein containing this binds to becteriochlorophyllide d, 
    which is the precursor to bacteriochlorophyl d. Therefore, it is likely a false annotation."
  ), collapse = ";"),
  paste0(c(
    "Magnesium (bacteriochlorophyll e axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "189438",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "189438"
    )), collapse = ","),
    "This is only a 2-star entry in ChEBI so technically it does not exist properly."
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a' A1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "189419",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "189419"
    )), collapse = ","),
    "This is only a 2-star entry in ChEBI so technically it does not exist properly."
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a A3 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a B1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a B3 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a ChlzD1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a ChlzD2 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a PD1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-a PD2 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll-b axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "61721",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61721"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll a axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58416",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58416"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (chlorophyll axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "28966",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "28966"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (divinyl chlorophyll-a' A1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "189420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "189420"
    )), collapse = ","),
    "This is only a 2-star entry in ChEBI so technically it does not exist properly."
  ), collapse = ";"),
  paste0(c(
    "Magnesium (divinyl chlorophyll-a A3 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "73095",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73095"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (divinyl chlorophyll-a B1 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "73095",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73095"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Magnesium (divinyl chlorophyll-a B3 axial ligand)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "73095",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73095"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Manganese",
    "25155",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25155"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Mercury",
    "25197",
    paste0(c(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25197"
    )), "16170"), collapse = ","),
    "",
    "",
    "Sub ID also contains elemental mercury"
  ), collapse = ";"),
  paste0(c(
    "Metal cation",
    "25213",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25213"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Mn in calcium-manganese-oxide [Ca-4Mn-5O]",
    "25155",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25155"
    )), collapse = ","),
    "calcium-manganese-oxide [Ca-4Mn-5O]",
    "",
    "There is no ChEBI ID yet refering to the calcium-manganese-oxide cluster."
  ), collapse = ";"),
  paste0(c(
    "Molybdate",
    "37239",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ","),
    "36264",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "36264"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum-iron-sulfur-carbon (7Fe-Mo-9S-C-homocitryl)",
    "24875|37239",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "37239"
      )), collapse = ",")),
      collapse = "|"),
    "30409",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30409"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum (Mo-bis(molybdopterin guanine dinucleotide) metal ligand)",
    "37239",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ","),
    "60539",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60539"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum (Mo-molybdopterin cytosine dinucleotide metal ligand)",
    "49414",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49414"
    )), collapse = ","),
    "71308",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "71308"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum (Mo-molybdopterin guanine dinucleotide metal ligand)",
    "49414",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49414"
    )), collapse = ","),
    "71310",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "71310"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum (Mo-molybdopterin metal ligand)",
    "37239",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ","),
    "71302",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "71302"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Monovalent cation",
    "60242",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60242"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Ni(2+)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel",
    "25516",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25516"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel-iron-sulfur",
    "25516|24875",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "25516"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
      )), collapse = ",")),
      collapse = "|"),
    "60400",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60400"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel-iron-sulfur (Ni-4Fe-4S)",
    "25516|24875",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "25516"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
      )), collapse = ",")),
      collapse = "|"),
    "47739",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "47739"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel-iron-sulfur (Ni-4Fe-5S)",
    "25516|24875",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "25516"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
      )), collapse = ",")),
      collapse = "|"),
    "177874",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "177874"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel (coenzyme f430 axial ligand)",
    "25516",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25516"
    )), collapse = ","),
    "28265",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "28265"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel (coenzyme F430 axial ligand)",
    "25516",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25516"
    )), collapse = ","),
    "28265",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "28265"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Nickel (Ni(II)-pyridinium-3,5-bisthiocarboxylate mononucleotide metal ligand)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "137373",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "137373"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Potassium",
    "29103",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29103"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Sodium",
    "29101",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29101"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Tungsten (W-bis(molybdopterin guanine dinucleotide) metal ligand)",
    "60401",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60401"
    )), collapse = ","),
    "60537",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60537"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Vanadium",
    "35172",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35172"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Vanadium-iron-sulfur-carbon (7Fe-V-9S-C-homocitryl)",
    "35172|24875",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "35172"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
      )), collapse = ",")),
      collapse = "|"),
    "60357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60357"
    )), collapse = ","),
    ""
  ), collapse = ";"),
  paste0(c(
    "Zinc",
    "63056",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "63056"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Zn(2+)",
    "29105",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29105"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "ferroheme b(2-)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60344",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60344"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "coenzyme F430(5-)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "60540",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60540"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "heme d cis-diol(2-)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "62814",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "62814"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "5-hydroxybenzimidazolylcob(I)amide(1-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60494",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60494"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "dihydrogenvanadate",
    "35172",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35172"
    )), collapse = ","),
    "35169",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35169"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "Mo(VI)-oxido Se-molybdopterin cytosine dinucleotide(4-)",
    "37239",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ","),
    "73094",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73094"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "cob(I)alamin(1-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60488",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60488"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "divinyl chlorophyll b(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "73096",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73096"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "bis(molybdopterin)tungsten cofactor",
    "60401|18420",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "60401"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "18420"
      )), collapse = ",")),
      collapse = "|"),
    "30402",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30402"
    )), collapse = ","),
    "From cc_cofactor"
  ), collapse = ";"),
  paste0(c(
    "tetrahydroxoborate(1-)",
    "22909",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "22909"
    )), collapse = ","),
    "41132",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "41132"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "tetrahydroxoborate(1-)",
    "22909",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "22909"
    )), collapse = ","),
    "41132",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "41132"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "methylarsonous acid",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "17826",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "17826"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "methylarsonate(1-)",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "33409",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33409"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "chlorophyll(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "139291",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "139291"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "chlorophyllide(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "139292",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "139292"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-sirohydrochlorin(8-)",
    "48828",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48828"
    )), collapse = ","),
    "60049",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60049"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "methyl-Co(2+)",
    "48828",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48828"
    )), collapse = ","),
    "85035",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "85035"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "protochlorophyllide(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "83350",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83350"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "arsenate(2-)",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "48597",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "48597"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "iron(III) oxide-hydroxide(1-)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "78619",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "78619"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "ferroheme o(2-)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "60530",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60530"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "ferroheme a(2-)",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "61715",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61715"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "ferrienterobactin(3-)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "28199",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "28199"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "magnesium protoporphyrin(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "60492",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60492"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "tellurite",
    "60271",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60271"
    )), collapse = ","),
    "30477",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30477"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "methanetelluronate(1-)",
    "60271",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60271"
    )), collapse = ","),
    "71624",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "71624"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcobinamide",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "2480",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "2480"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcobinamide phosphate(1-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58502",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58502"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcobinamide guanosyl diphosphate(1-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60487",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60487"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cob(II)yrinic acid a,c diamide(4-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58537",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58537"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cob(II)yrinate(6-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58894",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58894"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "ferroheme c di-L-cysteine(2-) residue",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "83739",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83739"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcobalamin 5'-phosphate(2-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60493",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60493"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "2,4-divinyl protochlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "58632",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58632"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "aquacob(III)alamin",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "15852",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "15852"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cyanocob(III)alamin",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "17439",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "17439"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "7(1)-hydroxychlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "83357",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83357"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "dimethylarsinate",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "16223",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "16223"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "7(1)-hydroxychlorophyll a(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "83377",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83377"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "beta-hematin",
    "29033",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29033"
    )), collapse = ","),
    "55377",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "55377"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "bacteriochlorophyllide c(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90965",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90965"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "Mo(=O)(=S)-molybdopterin cofactor(2-)",
    "37239",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "37239"
    )), collapse = ","),
    "82685",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "82685"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "magnesium protoporphyrin 13-monomethyl ester(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "60491",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60491"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "nickel-sirohydrochlorin(8-)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "136841",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "136841"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "nickel-sirohydrochlorin a,c-diamide(6-)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "136887",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "136887"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "bacteriochlorophyllide e(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "136512",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "136512"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "bacteriochlorophyllide f(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "139237",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "139237"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-7(7-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "70791",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "70791"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-6B(7-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "72780",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "72780"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-4(7-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60061",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60061"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt(II)-factor III(8-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "73299",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "73299"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcob(III)yrinate a,c-diamide(4-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58503",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58503"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "3-vinylbacteriochlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "83373",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83373"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "3-acetylchlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90794",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90794"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "bacteriochlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90795",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90795"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "3-(1-hydroxyethyl)bacteriochlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90791",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90791"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "3-(1-hydroxyethyl)chlorophyllide a(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90792",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90792"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "Co-methyl-Co-5-hydroxybenzimidazolylcob(III)amide",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "16379",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "16379"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "8,12-diethyl-3-vinylbacteriochlorophyllide d(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90964",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90964"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "12-ethyl-8-propyl-3-vinylbacteriochlorophyllide d(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90966",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90966"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "12-ethyl-8-isobutyl-3-vinylbacteriochlorophyllide d(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90967",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90967"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "3-vinylbacteriochlorophyllide d(1-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "90963",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "90963"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-8(6-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "70792",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "70792"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt(II)-factor IV(6-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "85471",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "85471"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "tungstate",
    "49955",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49955"
    )), collapse = ","),
    "46502",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "46502"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-5B(8-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60063",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60063"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-6A(7-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60064",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60064"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-5A(7-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60062",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60062"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "arseno-mycothiol(2-)",
    "35827",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35827"
    )), collapse = ","),
    "59655",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "59655"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-2(8-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60053",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60053"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "cobalt-precorrin-3(8-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "60060",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60060"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "15,17(3)-seco-F430-17(3)-acid(6-)",
    "49786",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49786"
    )), collapse = ","),
    "136888",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "136888"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "primary alkymercury(1+)",
    "25197",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25197"
    )), collapse = ","),
    "83725",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "83725"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "adenosylcobyrate",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58504",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58504"
    )), collapse = ","),
    "From cc_catalytic_activity"
  ), collapse = ";"),
  paste0(c(
    "chrysobactin",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "61345",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61345"
    )), collapse = ","),
    "From QuickGO, ChEBI does not contain a metal but is able to bind one. Complexes are usually not 1:1 but multiple ligands are required to bind one Fe3+ (https://doi.org/10.1007/BF01079695)."
  ), collapse = ";"),
  paste0(c(
    "achromobactin",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "61346",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "61346"
    )), collapse = ","),
    "From QuickGO, ChEBI does not contain a metal but is able to bind one."
  ), collapse = ";"),
  paste0(c(
    "enterobactin(1-)",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "77805",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "77805"
    )), collapse = ","),
    "From QuickGO, ChEBI does not contain a metal but is able to bind one."
  ), collapse = ";"),
  paste0(c(
    "Fe3S4 iron-sulfur cluster",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "64606",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "64606"
    )), collapse = ","),
    "From QuickGO, ChEBI has no formula."
  ), collapse = ";"),
  paste0(c(
    "MgATP(2-)",
    "18420",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "18420"
    )), collapse = ","),
    "30617",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30617"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "cluster",
    "25213",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "25213"
    )), collapse = ","),
    "33731",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33731"
    )), collapse = ","),
    "From QuickGO, ChEBI has no formula."
  ), collapse = ";"),
  paste0(c(
    "antimonous acid",
    "49867",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49867"
    )), collapse = ","),
    "49870",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "49870"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "silicic acid",
    "30584",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "30584"
    )), collapse = ","),
    "26675",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26675"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "chromate(2-)",
    "33516",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33516"
    )), collapse = ","),
    "35404",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "35404"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "iron chelate",
    "24875",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "24875"
    )), collapse = ","),
    "5975",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "5975"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "N',N'',N'''-triacetylfusarinine C",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "60481",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60481"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
    paste0(c(
      "iron(III) hydroxamate",
      "29034",
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "29034"
      )), collapse = ","),
      "28163",
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "28163"
      )), collapse = ","),
      "From QuickGO."
    ), collapse = ";"),
      paste0(c(
        "achromobactin",
        "29034",
        paste0(unlist(protti:::find_all_subs(
          data = chebi_relations,
          ids = "29034"
        )), collapse = ","),
        "61346",
        paste0(unlist(protti:::find_all_subs(
          data = chebi_relations,
          ids = "61346"
        )), collapse = ","),
        "From QuickGO."
      ), collapse = ";"),
  paste0(c(
    "calcium oxalate",
    "29108",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29108"
    )), collapse = ","),
    "60579",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60579"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "hydroxylapatite",
    "29108",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29108"
    )), collapse = ","),
    "52255",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "52255"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "selenate",
    "60250",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "60250"
    )), collapse = ","),
    "15075",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "15075"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "siderophore",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "26672",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "26672"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "cob(I)yrinate a,c diamide(5-)",
    "23336",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "23336"
    )), collapse = ","),
    "58575",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "58575"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "ferric enterobactin",
    "29034",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "29034"
    )), collapse = ","),
    "144426",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "144426"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";"),
  paste0(c(
    "boric acid",
    "33610",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33610"
    )), collapse = ","),
    "33118",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "33118"
    )), collapse = ","),
    "From QuickGO."
  ), collapse = ";")
))

fmb_annotation_uniprot <- fmb_annotation_uniprot_manuall %>%
  separate(fmb_name, into = c(
    "fmb_name",
    "metal_id",
    "metal_sub_id",
    "complex_id",
    "complex_sub_id",
    "comment"
  ), sep = ";")

not_yet_in_data_frame <- extracted_fmb_uniprot$fmb_name[!extracted_fmb_uniprot$fmb_name %in% fmb_annotation_uniprot$fmb_name]

if(length(not_yet_in_data_frame) > 0){
  message("The following IDs have not yet been manually annotated. Please add them to the list:\n", 
          paste0(not_yet_in_data_frame, collapse = ", "))
}

usethis::use_data(fmb_annotation_uniprot, overwrite = TRUE)