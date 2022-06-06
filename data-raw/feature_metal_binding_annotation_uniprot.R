# This is version 1, which was created on 22/06/03

# This code annotates all possible names that can appear in the
# feature(METAL BINDING) column in UniProt with ChEBI IDs that match.
# When updated also update the documentation!

library(tidyverse)
library(protti)

# First we retrieve all reviewed annotations from UniProt
# The relevant columns are: feature(METAL BINDING)

url_fmb_uniprot <- utils::URLencode("http://www.uniprot.org/uniprot/?query=reviewed:yes&format=tab&columns=id,feature(METAL BINDING)")
input_fmb_uniprot <- protti:::try_query(url_fmb_uniprot, progress = FALSE, show_col_types = FALSE)
colnames(input_fmb_uniprot) <- janitor::make_clean_names(c("id", "feature(METAL BINDING)"))

# Extract all possible names

extracted_fmb_uniprot <- input_fmb_uniprot %>%
  tidyr::drop_na(feature_metal_binding) %>%
  dplyr::mutate(fmb_name = stringr::str_extract_all(
    feature_metal_binding,
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
    "22633",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "22633"
    )), collapse = ","),
    "",
    "",
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
    "36264",
    paste0(unlist(protti:::find_all_subs(
      data = chebi_relations,
      ids = "36264"
    )), collapse = ","),
    "",
    "",
    ""
  ), collapse = ";"),
  paste0(c(
    "Molybdenum-iron-sulfur-carbon (7Fe-Mo-9S-C-homocitryl)",
    "37239|24875",
    paste0(
      c(paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "37239"
      )), collapse = ","),
      paste0(unlist(protti:::find_all_subs(
        data = chebi_relations,
        ids = "24875"
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
  ), collapse = ";")
))

fmb_annotation_uniprot <- fmb_annotation_uniprot_manuall %>%
  separate(fmb_name, into = c(
    "fmb_name",
    "metal_id",
    "metal_sub_id",
    "prosthetic_group_id",
    "prosthetic_group_sub_id",
    "comment"
  ), sep = ";")

not_yet_in_data_frame <- extracted_fmb_uniprot$fmb_name[!extracted_fmb_uniprot$fmb_name %in% fmb_annotation_uniprot$fmb_name]

if(length(not_yet_in_data_frame) > 0){
  message("The following IDs have not yet been manually annotated. Please add them to the list:\n", 
          paste0(not_yet_in_data_frame, collapse = ", "))
}

usethis::use_data(fmb_annotation_uniprot, overwrite = TRUE)