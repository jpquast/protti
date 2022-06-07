# Metal gene ontology slim subset
# This is version 1, which was created on 22/06/03

# This code generates a gene ontology subset that contains 
# metal related entries. The slimming process provided by 
# the QuickGO database is used.

library(tidyverse)
library(protti)
library(igraph)

# First we retrieve all possible GO terms for annotation purposes

terms <- fetch_quickgo(type = "terms")

# Manually define IDs that should be used for the slimming process

parent_metal_ids <- c("GO:0046872", # metal ion binding
                      "GO:0008237", # metallopeptidase activity
                      "GO:0140487", # metal ion sequestering activity
                      "GO:0004023", # alcohol dehydrogenase activity, metal ion-independent
                      "GO:0030946", # protein tyrosine phosphatase activity, metal-dependent
                      "GO:0051002", # ligase activity, forming nitrogen-metal bonds
                      "GO:0016722", # oxidoreductase activity, acting on metal ions
                      "GO:1904517", # MgATP(2-) binding
                      "GO:0051540", # metal cluster binding
                      "GO:0044590", # iron-sulfur-molybdenum cofactor binding
                      "GO:0140132", # iron-sulfur cluster carrier activity
                      "GO:0020037", # heme binding
                      "GO:0140488", # heme receptor activity
                      "GO:0015232", # heme transmembrane transporter activity
                      "GO:0016168", # chlorophyll binding
                      "GO:0140784", # metal ion sensor activity
                      "GO:0016530", # metallochaperone activity
                      "GO:0046873", # metal ion transmembrane transporter activity
                      "GO:0015343", # siderophore transmembrane transporter activity
                      "GO:0015621", # ferric triacetylfusarinine C transmembrane transporter activity
                      "GO:0032523", # silicon efflux transmembrane transporter activity
                      "GO:0080139", # borate efflux transmembrane transporter activity
                      "GO:0046715", # active borate transmembrane transporter activity
                      "GO:0015603", # iron chelate transmembrane transporter activity
                      "GO:0015105", # arsenite transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0015098", # molybdate ion transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0015654", # tellurite transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0015109", # chromate transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0015104", # antimonite transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0015115", # silicate transmembrane transporter activity (metal ion transmembrane transporter activity)
                      "GO:0030973", # molybdate ion binding (part of anion binding)
                      "GO:1901359", # tungstate binding (part of anion binding)
                      "GO:0046714", # borate binding (part of anion binding)
                      "GO:0046904", # calcium oxalate binding
                      "GO:0046848", # hydroxyapatite binding
                      "GO:0015624", # ABC-type ferric-enterobactin transporter activity
                      "GO:0140481", # ABC-type iron-sulfur cluster transporter activity
                      "GO:1901238", # ABC-type tungstate transporter activity
                      "GO:0015625", # ABC-type ferric hydroxamate transporter activity
                      "GO:0015345" # ferric enterobactin:proton symporter activity
                      )

# Retrieve slim dataset 

metal_slim_subset <- fetch_quickgo(type = "slims", go_id_slims = parent_metal_ids)

# Annotate slim subset

metal_slim_subset_annotated <- metal_slim_subset %>% 
  left_join(terms, by = c("slims_from_id" = "main_id")) %>% 
  filter(ontology == "molecular_function") # filter for molecular function ontology

# Find terms that still contain metal-related ChEBI IDs that are not yet covered by the subset dataset
# Retrieve ChEBI IDs in order to find all metal related ChEBI entries

# As defined in metal_chebi_uniprot.R
metal_chebi_ids_wo_formula <- c(25213, # a metal cation
                                30408, # iron-sulfur cluster
                                30413, # heme
                                38201, # a bacteriochlorophyll
                                60240, # a divalent metal cation
                                60242, # a monovalent cation (this can also be non-metal so should be handled with care)
                                60400 # [Ni-Fe-S] cluster
)

metal_chebi <- protti::fetch_chebi() %>% 
  dplyr::distinct(id, chebi_accession, formula) %>% 
  dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) | 
                  id %in% metal_chebi_ids_wo_formula)

# If terms_metal contains any rows add parent GO terms to parent_metal_ids above!
terms_metal <- terms %>% 
  filter(ontology == "molecular_function") %>% # only use molecular function since it only true reports on direct interactions
  filter(chebi_id != "CHEBI:60242") %>% # exclude monovalent cations since it is not specific
  filter(chebi_id %in% unique(metal_chebi$chebi_accession)) %>%
  filter(!main_id %in% metal_slim_subset_annotated$slims_from_id)

# Assign a ChEBI ID to all GO terms.
# If child terms do not have a ChEBI ID they take over the ChEBI ID from their parent

# The igraph package is used for this

go_network <- metal_slim_subset_annotated %>% 
  distinct(slims_from_id, child_id, chebi_id, relations_relation, relations_url, database, relations_term) %>% 
  group_by(slims_from_id) %>% 
  mutate(chebi_id = paste0(chebi_id, collapse = ";"),
         relations_relation = paste0(relations_relation, collapse = ";"),
         relations_url = paste0(relations_url, collapse = ";"),
         database = paste0(database, collapse = ";"),
         relations_term = paste0(relations_term, collapse = ";")) %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(child_id = str_split(child_id, pattern = ";")) %>% 
  unnest(child_id) %>% 
  mutate(chebi_id = ifelse(chebi_id == "NA", NA_character_, chebi_id),
         relations_relation = ifelse(relations_relation == "NA", NA_character_, relations_relation),
         relations_url = ifelse(relations_url == "NA", NA_character_, relations_url),
         database = ifelse(database == "NA", NA_character_, database),
         relations_term = ifelse(relations_term == "NA", NA_character_, relations_term)) 

go_network_connections <- go_network %>% 
  filter(child_id %in% metal_slim_subset_annotated$slims_from_id) %>% 
  distinct(slims_from_id, child_id)

# Custom function to propagate chebi_ids 
# Caution, this function cannot deal with two IDs converging on one unannotated note. It will take over 
# the ID of the parent note with the least children.
# One entry not correct is: GO:0015277 has only sodium annotation even though there is also potassium

propagate_attribute <- function(g, attribute_name = NULL){
  ids_not_na <- V(g)[!is.na(eval(bquote(`$`(V(g), .(as.name(attribute_name))))))]
  
  sub_ids_ids_not_na <- map(.x = ids_not_na,
                            .f = ~{
                              igraph::subcomponent(g, .x, "out")
                            })
  
  sub_ids_ids_not_na
  
  deleted_ids <- list()
  
  for(i in 1:length(sub_ids_ids_not_na)){
    current_id <- sub_ids_ids_not_na[[i]]
    current_id_length <- length(current_id)
    
    for(j in 1:length(sub_ids_ids_not_na)){
      sub_id <- sub_ids_ids_not_na[[j]]
      if(current_id_length > length(sub_id)){
        current_id <- current_id[!current_id %in% sub_id]
      }
    }
    
    deleted_ids[[i]] <- current_id
  }
  
  deleted_ids
  
  to_replace <- V(g)[is.na(eval(bquote(`$`(V(g), .(as.name(attribute_name))))))]
  
  for(i in 1:length(to_replace)){
    print(to_replace[i])
    new_value <- NA
    for(j in 1:length(deleted_ids)){
      if(to_replace[i] %in% deleted_ids[[j]]){
        new_value <- eval(bquote(`$`(V(g)[ids_not_na[j]], .(as.name(attribute_name)))))
      }
    }
    g <- set_vertex_attr(graph = g, name = attribute_name, index = to_replace[i], value = new_value)
  }
  
  return(g)
}

############### Chebi IDs

go_network_vertex_attribute_chebi_id <- go_network %>% 
  distinct(slims_from_id, chebi_id) 

# Generate graph
g_chebi_id <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_chebi_id) 

# Update graph

g_chebi_id_new <- propagate_attribute(g = g_chebi_id, attribute_name = "chebi_id")

propagated_chebi_ids <- as_data_frame(g_chebi_id_new, what="vertices") %>% 
  rename(slims_from_id = name)

############### relations_relation

go_network_vertex_attribute_relations_relation <- go_network %>% 
  distinct(slims_from_id, relations_relation) 

# Generate graph
g_relations_relation <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_relation) 

# Update graph

g_relations_relation_new <- propagate_attribute(g = g_relations_relation, attribute_name = "relations_relation")

propagated_relations_relation <- as_data_frame(g_relations_relation_new, what="vertices") %>% 
  rename(slims_from_id = name)

############### relations_url

go_network_vertex_attribute_relations_url <- go_network %>% 
  distinct(slims_from_id, relations_url) 

# Generate graph
g_relations_url <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_url) 

# Update graph

g_relations_url_new <- propagate_attribute(g = g_relations_url, attribute_name = "relations_url")

propagated_relations_urls <- as_data_frame(g_relations_url_new, what="vertices") %>% 
  rename(slims_from_id = name)

############### database

go_network_vertex_attribute_database <- go_network %>% 
  distinct(slims_from_id, database) 

# Generate graph
g_database <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_database) 

# Update graph

g_database_new <- propagate_attribute(g = g_database, attribute_name = "database")

propagated_databases <- as_data_frame(g_database_new, what="vertices") %>% 
  rename(slims_from_id = name)

############### relations_term

go_network_vertex_attribute_relations_term <- go_network %>% 
  distinct(slims_from_id, relations_term) 

# Generate graph
g_relations_term <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_term) 

# Update graph

g_relations_term_new <- propagate_attribute(g = g_relations_term, attribute_name = "relations_term")

propagated_relations_terms <- as_data_frame(g_relations_term_new, what="vertices") %>% 
  rename(slims_from_id = name)

# Create final annotation set 

metal_go_slim_subset_pre <- metal_slim_subset_annotated %>% 
  select(-c(chebi_id, relations_relation, relations_url, database, relations_term)) %>% 
  distinct() %>% 
  left_join(propagated_chebi_ids, by = "slims_from_id") %>% 
  left_join(propagated_relations_relation, by = "slims_from_id") %>% 
  left_join(propagated_relations_urls, by = "slims_from_id") %>% 
  left_join(propagated_databases, by = "slims_from_id") %>% 
  left_join(propagated_relations_terms, by = "slims_from_id") %>% 
  mutate(chebi_id = str_split(chebi_id, pattern = ";"),
         relations_relation = str_split(relations_relation, pattern = ";"),
         relations_url = str_split(relations_url, pattern = ";"),
         database = str_split(database, pattern = ";"),
         relations_term = str_split(relations_term, pattern = ";")) %>% 
  unnest(c(chebi_id, relations_relation, relations_url, database, relations_term))

# Create another list that contains manually annotated ChEBI IDs that are considered metal even though they do not contain a formula
# There are many siderophores. Their formula is often without metal but they are counted as metal associated ChEBI ID.
# First need to create non_metal_chebis then the list below can be created with the commented code.
# message('"', paste0(unique(non_metal_chebis$relations_term), collapse = '",\n"'), '"')
more_metal_ids_wo_formula <- c("cadmium cation",
                               "chrysobactin",
                               "achromobactin",
                               "enterobactin(1-)",
                               "silver cation",
                               "zinc cation",
                               "Fe3S4 iron-sulfur cluster",
                               "Fe4S4 iron-sulfur cluster",
                               "cluster",
                               "magnesium ion",
                               "iron chelate",
                               "transition element cation",
                               "alkali metal cation",
                               "Fe4S3 iron-sulfur cluster",
                               "lead cation",
                               "mercury cation",
                               "selenate")

metal_chebis <- metal_go_slim_subset_pre %>% 
  filter((chebi_id %in% unique(metal_chebi$chebi_accession)) | 
           (as.numeric(str_replace(chebi_id, pattern = "CHEBI:", replacement = "")) %in% metal_chebi_ids_wo_formula) |
           (relations_term %in% more_metal_ids_wo_formula)) 

non_metal_chebis <- metal_go_slim_subset_pre %>% 
  filter(!(chebi_id %in% unique(metal_chebi$chebi_accession)) & 
           !(as.numeric(str_replace(chebi_id, pattern = "CHEBI:", replacement = "")) %in% metal_chebi_ids_wo_formula) &
           !(relations_term %in% more_metal_ids_wo_formula) &
           !(slims_from_id %in% metal_chebis$slims_from_id)) 

# Manual annotation of remaining GO IDs

message(
  paste0(
    paste0('paste0("', 
           unique(non_metal_chebis$slims_from_id), ";", unique(non_metal_chebis$main_name), '", ";", "CHEBI")'), collapse = ',\n'))

manual_go_chebi_annotation <- data.frame(name = c(paste0("GO:0052851;ferric-chelate reductase (NADPH) activity", ";", "CHEBI:29034|CHEBI:29033"),
                                                  paste0("GO:0030946;protein tyrosine phosphatase activity, metal-dependent", ";", "CHEBI:25213"),
                                                  paste0("GO:0061473;murein tripeptide carboxypeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0050453;cob(II)alamin reductase activity", ";", "CHEBI:60488|CHEBI:16304|CHEBI:48828|CHEBI:49415"),
                                                  paste0("GO:0140492;metal-dependent deubiquitinase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0000293;ferric-chelate reductase activity", ";", "GO:0000293"),
                                                  paste0("GO:0140487;metal ion sequestering activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0140486;zinc ion sequestering activity", ";", "CHEBI:29105"),
                                                  paste0("GO:0008823;cupric reductase activity", ";", "	CHEBI:49552|CHEBI:29036"),
                                                  paste0("GO:0004222;metalloendopeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0070006;metalloaminopeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0004322;ferroxidase activity", ";", "CHEBI:29034|CHEBI:29033"),
                                                  paste0("GO:0015344;siderophore uptake transmembrane transporter activity", ";", "CHEBI:26672|CHEBI:29034"),
                                                  paste0("GO:0015343;siderophore transmembrane transporter activity", ";", "CHEBI:26672|CHEBI:29034"),
                                                  paste0("GO:0004023;alcohol dehydrogenase activity, metal ion-independent", ";", "CHEBI:25213"),
                                                  paste0("GO:0051002;ligase activity, forming nitrogen-metal bonds", ";", "CHEBI:25213"),
                                                  paste0("GO:0051003;ligase activity, forming nitrogen-metal bonds, forming coordination complexes", ";", "CHEBI:25213"),
                                                  paste0("GO:0004181;metallocarboxypeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0030586;[methionine synthase] reductase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0016152;mercury (II) reductase activity", ";", "CHEBI:25197"),
                                                  paste0("GO:0051116;cobaltochelatase activity", ";", "CHEBI:58537|CHEBI:48828"),
                                                  paste0("GO:0008235;metalloexopeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0008237;metallopeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0047852;diferric-transferrin reductase activity", ";", "CHEBI:29034|CHEBI:29033"),
                                                  paste0("GO:0043783;oxidoreductase activity, acting on metal ions, flavin as acceptor", ";", "CHEBI:25213"),
                                                  paste0("GO:0043784;cob(II)yrinic acid a,c-diamide reductase activity", ";", "CHEBI:58575|CHEBI:58537|CHEBI:48828|CHEBI:85033"),
                                                  paste0("GO:0140618;ferric-chelate reductase (NADH) activity", ";", "CHEBI:26672|CHEBI:29034|CHEBI:29033"),
                                                  paste0("GO:0033787;cyanocobalamin reductase (cyanide-eliminating) activity", ";", "CHEBI:60488|CHEBI:17439|CHEBI:49415|CHEBI:85033"),
                                                  paste0("GO:0140758;metal-dependent deNEDDylase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0009046;zinc D-Ala-D-Ala carboxypeptidase activity", ";", "CHEBI:29105"),
                                                  paste0("GO:0008191;metalloendopeptidase inhibitor activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0102274;glutathione S-conjugate carboxypeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0070573;metallodipeptidase activity", ";", "CHEBI:25213"),
                                                  paste0("GO:0016722;oxidoreductase activity, acting on metal ions", ";", "CHEBI:25213"),
                                                  paste0("GO:0016723;oxidoreductase activity, acting on metal ions, NAD or NADP as acceptor", ";", "CHEBI:25213"),
                                                  paste0("GO:0016724;oxidoreductase activity, acting on metal ions, oxygen as acceptor", ";", "CHEBI:25213"),
                                                  paste0("GO:0140315;iron ion sequestering activity", ";", "CHEBI:24875"),
                                                  paste0("GO:0140314;calcium ion sequestering activity", ";", "CHEBI:39123"),
                                                  paste0("GO:0016851;magnesium chelatase activity", ";", "CHEBI:18420"),
                                                  paste0("GO:0140571;transmembrane ascorbate ferrireductase activity", ";", "CHEBI:29034|CHEBI:29033"),
                                                  paste0("GO:0015620;ferric-enterobactin transmembrane transporter activity", ";", "CHEBI:144426|CHEBI:29034"),
                                                  paste0("GO:1902945;metalloendopeptidase activity involved in amyloid precursor protein catabolic process", ";", "CHEBI:25213"))) %>% 
  separate(name, into = c("slims_from_id", "main_name", "chebi_id"), sep = ";") %>% 
  mutate(chebi_id = str_split(chebi_id, pattern = "\\|")) %>% 
  unnest(chebi_id) %>% 
  left_join(distinct(non_metal_chebis, slims_from_id, slims_to_ids), by = "slims_from_id")

metal_go_slim_subset <- metal_chebis %>% 
  bind_rows(manual_go_chebi_annotation) %>% 
  distinct(slims_from_id, slims_to_ids, main_name, chebi_id, relations_relation, relations_term, database) %>% 
  mutate(database = ifelse(is.na(database), "manual", database))

usethis::use_data(metal_go_slim_subset, overwrite = TRUE)