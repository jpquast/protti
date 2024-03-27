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

parent_metal_ids <- c(
  "GO:0046872", # metal ion binding
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
  "GO:0015345", # ferric enterobactin:proton symporter activity
  "GO:1903981", # enterobactin binding
  "GO:0004076", # biotin synthase activity
  "GO:0016041", # glutamate synthase (ferredoxin) activity
  "GO:0018695", # 4-cresol dehydrogenase (hydroxylating) activity
  "GO:0018694", # p-cymene methyl hydroxylase activity
  "GO:0018683", # camphor 5-monooxygenase activity
  "GO:0018685", # alkane 1-monooxygenase activity
  "GO:0043757", # adenosylcobinamide-phosphate synthase activity
  "GO:0043756", # adenosylcobinamide hydrolase activity
  "GO:0043779", # cobalt-precorrin-5A acetaldehyde-lyase activity
  "GO:0004130", # cytochrome-c peroxidase activity
  "GO:0043776", # cobalt-precorrin-6B C5-methyltransferase activity
  "GO:0043778", # cobalt-precorrin-8 methylmutase activity
  "GO:0043777", # cobalt-precorrin-7 C15-methyltransferase activity
  "GO:0004128", # cytochrome-b5 reductase activity, acting on NAD(P)H
  "GO:0102842", # 1-18:1-2-16:2-monogalactosyldiacylglycerol desaturase activity (SN2-16:3 forming)
  "GO:0102843", # 1-18:2-2-16:0-monogalactosyldiacylglycerol desaturase activity (SN2-16:1 forming)
  "GO:0102850", # 1-18:1-2-16:0-phosphatidylglycerol omega-6 desaturase activity
  "GO:0102859", # 1-18:1-2-18:2-phosphatidylcholine desaturase activity (SN2-18:3 forming)
  "GO:0102867", # molybdenum cofactor sulfurtransferase activity
  "GO:0102866", # di-homo-gamma-linolenate delta5 desaturase activity
  "GO:0102865", # delta6-acyl-lipid desaturase activity
  "GO:0004507", # steroid 11-beta-monooxygenase activity
  "GO:0018801", # glutaconyl-CoA decarboxylase activity
  "GO:0043805", # indolepyruvate ferredoxin oxidoreductase activity
  "GO:0043807", # 3-methyl-2-oxobutanoate dehydrogenase (ferredoxin) activity
  "GO:0018836", # alkylmercury lyase activity
  "GO:0043884", # CO-methylating acetyl-CoA synthase activity
  "GO:0043823", # spheroidene monooxygenase activity
  "GO:0018522", # benzoyl-CoA reductase activity
  "GO:0018525", # 4-hydroxybenzoyl-CoA reductase activity
  "GO:0018491", # 2-oxobutyrate synthase activity
  "GO:0018493", # formylmethanofuran dehydrogenase activity
  "GO:0018492", # carbon-monoxide dehydrogenase (acceptor) activity
  "GO:0102220", # hydrogenase activity (NAD+, ferredoxin)
  "GO:0004408", # holocytochrome-c synthase activity
  "GO:0004460", # L-lactate dehydrogenase (cytochrome) activity
  "GO:0004458", # D-lactate dehydrogenase (cytochrome) activity
  "GO:0004498", # calcidiol 1-monooxygenase activity
  "GO:0102988", # 9,12-cis-hexadecadienoic acid delta 15 desaturase activity
  "GO:0102985", # delta12-fatty-acid desaturase activity
  "GO:0102987", # palmitoleic acid delta 12 desaturase activity
  "GO:0102993", # linolenate delta15 desaturase activity
  "GO:0016966", # nitric oxide reductase activity
  "GO:0016630", # protochlorophyllide reductase activity
  "GO:0016633", # galactonolactone dehydrogenase activity
  "GO:0016689", # manganese peroxidase activity
  "GO:0043834", # trimethylamine methyltransferase activity
  "GO:0043833", # [methyl-Co(III) methylamine-specific corrinoid protein]:coenzyme M methyltransferase activity
  "GO:0043852", # monomethylamine methyltransferase activity
  "GO:0102771", # sphingolipid very long chain fatty acid alpha-hydroxylase activity
  "GO:0102786", # stearoyl-[acp] desaturase activity
  "GO:0043782", # cobalt-precorrin-3 C17-methyltransferase activity
  "GO:0043791", # dimethylamine methyltransferase activity
  "GO:0043797", # glyceraldehyde-3-phosphate dehydrogenase (ferredoxin) activity
  "GO:0004324", # ferredoxin-NADP+ reductase activity
  "GO:0004325", # ferrochelatase activity
  "GO:0004392", # heme oxygenase (decyclizing) activity
  "GO:0102654", # 1-18:1-2-16:0-phosphatidylglycerol trans-3 desaturase activity
  "GO:0004129", # cytochrome-c oxidase activity
  "GO:0043781", # cobalt-factor II C20-methyltransferase activity
  "GO:0043780", # cobalt-precorrin-5B C1-methyltransferase activity
  "GO:0016163", # nitrogenase activity
  "GO:0090523", # cytochrome-b5 reductase activity, acting on NADPH
  "GO:0004768", # stearoyl-CoA 9-desaturase activity
  "GO:0051073", # adenosylcobinamide-GDP ribazoletransferase activity
  "GO:0047152", # methanol-5-hydroxybenzimidazolylcobamide Co-methyltransferase activity
  "GO:0047111", # formate dehydrogenase (cytochrome-c-553) activity
  "GO:0047103", # 3-alpha,7-alpha,12-alpha-trihydroxycholestan-26-al 26-oxidoreductase activity
  "GO:0051266", # sirohydrochlorin ferrochelatase activity
  "GO:0016852", # sirohydrochlorin cobaltochelatase activity
  "GO:0016992", # lipoate synthase activity
  "GO:0003958", # NADPH-hemoprotein reductase activity
  "GO:0050421", # nitrite reductase (NO-forming) activity
  "GO:0015047", # NADPH-cytochrome-c2 reductase activity
  "GO:0015046", # rubredoxin-NADP+ reductase activity
  "GO:0015044", # rubredoxin-NAD+ reductase activity
  "GO:0015039", # NADPH-adrenodoxin reductase activity
  "GO:0052610", # beta-cryptoxanthin hydroxylase activity
  "GO:0052611", # beta-carotene 3-hydroxylase activity
  "GO:0052631", # sphingolipid delta-8 desaturase activity
  "GO:0052607", # 7-hydroxy-chlorophyllide a oxygenase activity
  "GO:0052606", # chlorophyllide a oxygenase activity
  "GO:0050046", # delta7-sterol 5(6)-desaturase activity
  "GO:0052662", # zeaxanthin epoxidase activity
  "GO:0033738", # methylenetetrahydrofolate reductase (ferredoxin) activity
  "GO:0008794", # arsenate reductase (glutaredoxin) activity
  "GO:0033728", # divinyl chlorophyllide a 8-vinyl-reductase activity
  "GO:0033726", # aldehyde ferredoxin oxidoreductase activity
  "GO:0030059", # aralkylamine dehydrogenase (azurin) activity
  "GO:0030338", # CMP-N-acetylneuraminate monooxygenase activity
  "GO:0030269", # tetrahydromethanopterin S-methyltransferase activity
  "GO:0042279", # nitrite reductase (cytochrome, ammonia-forming) activity
  "GO:0042284", # sphingolipid delta-4 desaturase activity
  "GO:0042242", # cobyrinic acid a,c-diamide synthase activity
  "GO:0140741", # tRNA U4 sulfurtransferase
  "GO:0017105", # acyl-CoA delta11-(Z)-desaturase activity
  "GO:0051743", # red chlorophyll catabolite reductase activity
  "GO:0051744", # 3,8-divinyl protochlorophyllide a 8-vinyl reductase activity
  "GO:0036199", # cholest-4-en-3-one 26-monooxygenase activity
  "GO:0047959", # glycine dehydrogenase (cytochrome) activity
  "GO:0008386", # cholesterol monooxygenase (side-chain-cleaving) activity
  "GO:0047783", # corticosterone 18-monooxygenase activity
  "GO:0047806", # cytochrome-c3 hydrogenase activity
  "GO:0070225", # sulfide dehydrogenase activity
  "GO:0047889", # ferredoxin-nitrate reductase activity
  "GO:0047898", # formate dehydrogenase (cytochrome) activity
  "GO:0008820", # cobinamide phosphate guanylyltransferase activity
  "GO:0008818", # cobalamin 5'-phosphate synthase activity
  "GO:0008901", # ferredoxin hydrogenase activity
  "GO:0102036", # methyltetrahydrofolate:corrinoid/iron-sulfur protein methyltransferase activity
  "GO:0050605", # superoxide reductase activity
  "GO:0050600", # acyl-CoA 11-(E)-desaturase activity
  "GO:0050618", # phycoerythrobilin:ferredoxin oxidoreductase activity
  "GO:0050619", # phytochromobilin:ferredoxin oxidoreductase activity
  "GO:0050617", # 15,16-dihydrobiliverdin:ferredoxin oxidoreductase activity
  "GO:0050612", # arsenate reductase (donor) activity
  "GO:0050610", # methylarsonate reductase activity
  "GO:0050611", # arsenate reductase (azurin) activity
  "GO:0050620", # phycocyanobilin:ferredoxin oxidoreductase activity
  "GO:0050626", # trimethylamine-N-oxide reductase (cytochrome c) activity
  "GO:0047595", # 6-hydroxynicotinate reductase activity
  "GO:0047553", # 2-oxoglutarate synthase activity
  "GO:0047527", # 2,3-dihydroxybenzoate-serine ligase activity
  "GO:0106344", # 4-amino-5-hydroxymethyl-2-methylpyrimidine phosphate synthase activity from histidine and PLP
  "GO:0106364", # 4-hydroxy-3-all-trans-hexaprenylbenzoate oxygenase activity
  "GO:0008860", # ferredoxin-NAD+ reductase activity
  "GO:0008849", # enterochelin esterase activity
  "GO:0033740", # hydroxylamine oxidoreductase activity
  "GO:0102100", # mycothiol-arsenate ligase activity
  "GO:0033797", # selenate reductase activity
  "GO:0102172", # 4alpha-hydroxymethyl,4beta,14alpha-dimethyl-9beta,19-cyclo-5alpha-ergost-24(241)-en-3beta-ol-4alpha-methyl oxidase activity
  "GO:0102178", # 4alpha-formyl-ergosta-7,24(241)-dien-3beta-ol-methyl oxidase activity
  "GO:0102177", # 24-methylenelophenol methyl oxidase activity
  "GO:0102179", # 24-ethylidenelophenol 4alpha-methyl oxidase activity
  "GO:0102174", # 4alpha-formyl,4beta,14alpha-dimethyl-9beta,19-cyclo-5alpha-ergost-24(241)-en-3beta-ol-4alpha-methyl oxidase activity
  "GO:0102173", # 24-methylenecycloartanol 4alpha-methyl oxidase activity
  "GO:0102181", # 4alpha-formyl-stigmasta-7,24(241)-dien-3beta-ol-methyl oxidase activity
  "GO:0102180", # 4alpha-hydroxymethyl-stigmasta-7,24(241)-dien-3beta-ol-methyl oxidase activity
  "GO:0008495", # protoheme IX farnesyltransferase activity
  "GO:0047726", # iron-cytochrome-c reductase activity
  "GO:0008121", # ubiquinol-cytochrome-c reductase activity
  "GO:0047749", # cholestanetriol 26-monooxygenase activity
  "GO:0047748", # cholestanetetraol 26-dehydrogenase activity
  "GO:0047059", # polyvinyl alcohol dehydrogenase (cytochrome) activity
  "GO:0047051", # D-lactate dehydrogenase (cytochrome c-553) activity
  "GO:1990088", # [methyl-Co(III) methanol-specific corrinoid protein]:coenzyme M methyltransferase
  "GO:0015420", # ABC-type vitamin B12 transporter activity
  "GO:0015446", # ATPase-coupled arsenite transmembrane transporter activity
  "GO:0050310", # sulfite dehydrogenase activity
  "GO:0050311", # sulfite reductase (ferredoxin) activity
  "GO:0050338", # thiosulfate dehydrogenase activity
  "GO:0052933", # alcohol dehydrogenase (cytochrome c(L)) activity
  "GO:0052934", # alcohol dehydrogenase (cytochrome c) activity
  "GO:0050304", # nitrous-oxide reductase activity
  "GO:0062185", # secalciferol 1-monooxygenase activity
  "GO:0062181", # 1-alpha,25-dihydroxyvitamin D3 23-hydroxylase activity
  "GO:0062180", # 25-hydroxycholecalciferol-23-hydroxylase activity
  "GO:0050183", # phosphatidylcholine 12-monooxygenase activity
  "GO:0050184", # phosphatidylcholine desaturase activity
  "GO:0050140", # nitrate reductase (cytochrome) activity
  "GO:0050087", # mannitol dehydrogenase (cytochrome) activity
  "GO:0050207", # plasmanylethanolamine desaturase activity
  "GO:0052876", # methylamine dehydrogenase (amicyanin) activity
  "GO:0061603", # molybdenum cofactor guanylyltransferase activity
  "GO:0061602", # molybdenum cofactor cytidylyltransferase activity
  "GO:0010277", # chlorophyllide a oxygenase [overall] activity
  "GO:0036200", # 3-ketosteroid 9-alpha-monooxygenase activity
  "GO:0000254", # C-4 methylsterol oxidase activity
  "GO:0036354", # 2-desacetyl-2-hydroxyethyl bacteriochlorophyllide a dehydrogenase activity
  "GO:0036428", # adenosylcobinamide kinase (GTP-specific) activity
  "GO:0036429", # adenosylcobinamide kinase (ATP-specific) activity
  "GO:0061599", # molybdopterin molybdotransferase activity
  "GO:0048529", # magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase activity
  "GO:0051921", # adenosylcobyric acid synthase (glutamine-hydrolyzing) activity
  "GO:0048307", # ferredoxin-nitrite reductase activity
  "GO:0061896", # all-trans retinol 3,4-desaturase activity
  "GO:0103012", # ferredoxin-thioredoxin reductase activity
  "GO:0019164", # pyruvate synthase activity
  "GO:0019133", # choline monooxygenase activity
  "GO:0046408", # chlorophyll synthetase activity
  "GO:0046406", # magnesium protoporphyrin IX methyltransferase activity
  "GO:0009496", # plastoquinol--plastocyanin reductase activity
  "GO:0046429", # 4-hydroxy-3-methylbut-2-en-1-yl diphosphate synthase activity
  "GO:0032441" # pheophorbide a oxygenase activity
)

# Retrieve slim dataset

metal_slim_subset <- fetch_quickgo(type = "slims", go_id_slims = parent_metal_ids)

# Annotate slim subset

metal_slim_subset_annotated <- metal_slim_subset %>%
  left_join(terms, by = c("slims_from_id" = "main_id")) %>%
  filter(ontology == "molecular_function") # filter for molecular function ontology

# Find terms that still contain metal-related ChEBI IDs that are not yet covered by the subset dataset
# Retrieve ChEBI IDs in order to find all metal related ChEBI entries

# Use directly from metal_chebi_uniprot which is defined in metal_chebi_uniprot.R
metal_chebi_ids_wo_formula <- setNames(
  metal_chebi_uniprot %>%
    dplyr::filter(is.na(formula)) %>%
    dplyr::distinct(metal_atom_id) %>%
    dplyr::pull(),
  metal_chebi_uniprot %>%
    dplyr::filter(is.na(formula)) %>%
    dplyr::distinct(id) %>%
    dplyr::pull()
)

# This vector contains additional ChEBI IDs that are metal related but do not contain a formula.
# This should be updated if promted bellow (around line 520)
more_metal_chebi_ids_wo_formula <- c(
  "63063" = "22977", # cadmium cation
  "61345" = "18248", # chrysobactin
  "61346" = "18248", # achromobactin
  "77805" = "18248", # enterobactin(1-)
  "60253" = "30512", # silver cation
  "63056" = "27363", # zinc cation
  "64606" = "18248", # Fe3S4 iron-sulfur cluster
  "64607" = "18248", # Fe4S4 iron-sulfur cluster
  "33731" = "33521", # cluster
  "5975" = "18248", # iron chelate
  "33515" = "27081", # transition element cation
  "33504" = "22314", # alkali metal cation
  "85620" = "18248", # Fe4S3 iron-sulfur cluster
  "60252" = "25016", # lead cation
  "25197" = "25195", # mercury cation
  "15075" = "27568", # selenate
  "26672" = "18248" # siderophore (was included in the manual annotation)
)

metal_chebi_ids_wo_formula <- c(metal_chebi_ids_wo_formula, more_metal_chebi_ids_wo_formula)

# Create annotation vector that indicates what ChEBI ID a certain element symbol should be annotated with
metal_chebi_annotation <- setNames(metal_list$chebi_id, metal_list$symbol)

metal_chebi <- protti::fetch_chebi(stars = c(2, 3)) %>%
  dplyr::distinct(id, chebi_accession, formula) %>%
  dplyr::filter(stringr::str_detect(formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])")) |
    id %in% as.numeric(names(metal_chebi_ids_wo_formula))) %>%
  dplyr::mutate(extract_formula = stringr::str_extract_all(.data$formula, pattern = paste0(paste0(metal_list$symbol, collapse = "(?![:lower:])|"), "(?![:lower:])"))) %>%
  dplyr::mutate(
    extract_formula = ifelse(
      as.character(.data$extract_formula) == "character(0)",
      NA,
      .data$extract_formula
    )
  ) %>%
  tidyr::unnest(extract_formula) %>%
  dplyr::mutate(metal_atom_id = ifelse(is.na(extract_formula),
    metal_chebi_ids_wo_formula[as.character(id)],
    metal_chebi_annotation[extract_formula]
  ))

# If terms_metal contains any rows add parent GO terms to parent_metal_ids above!
terms_metal <- terms %>%
  filter(ontology == "molecular_function") %>% # only use molecular function since it only true reports on direct interactions
  filter(chebi_id != "CHEBI:60242") %>% # exclude monovalent cations since it is not specific
  filter(chebi_id %in% unique(metal_chebi$chebi_accession)) %>%
  filter(!main_id %in% metal_slim_subset_annotated$slims_from_id)

terms_metal_id_name <- terms_metal %>%
  distinct(main_id, main_name)

terms_metal_paste <- paste0(paste0('"', terms_metal_id_name$main_id, '", \\# ', terms_metal_id_name$main_name), collapse = "\n")
# Assign a ChEBI ID to all GO terms.
# If child terms do not have a ChEBI ID (or ChEBI IDs that are not metals) they take over the ChEBI ID from their parent

# The igraph package is used for this

go_network <- metal_slim_subset_annotated %>%
  mutate(
    chebi_id = ifelse(!(chebi_id %in% metal_chebi$chebi_accession) | is.na(chebi_id), NA, chebi_id),
    relations_relation = ifelse(is.na(chebi_id), NA, relations_relation),
    relations_url = ifelse(is.na(chebi_id), NA, relations_url),
    database = ifelse(is.na(chebi_id), NA, database),
    relations_term = ifelse(is.na(chebi_id), NA, relations_term)
  ) %>%
  group_by(slims_from_id) %>%
  filter(all(is.na(chebi_id)) | !is.na(chebi_id)) %>%
  ungroup() %>%
  # exclude all non-metal IDs according to the above ChEBI list
  # this might however also exclude non-metal ChEBI IDs if they are present.
  # One can run this pipeline once commenting these arguments out and checking for new non-metal IDs
  # Then run it again with the arguments enabled.
  distinct(slims_from_id, child_id, chebi_id, relations_relation, relations_url, database, relations_term) %>%
  group_by(slims_from_id) %>%
  mutate(
    chebi_id = paste0(chebi_id, collapse = ";"),
    relations_relation = paste0(relations_relation, collapse = ";"),
    relations_url = paste0(relations_url, collapse = ";"),
    database = paste0(database, collapse = ";"),
    relations_term = paste0(relations_term, collapse = ";")
  ) %>%
  ungroup() %>%
  distinct() %>%
  mutate(child_id = str_split(child_id, pattern = ";")) %>%
  unnest(child_id) %>%
  mutate(
    chebi_id = ifelse(chebi_id == "NA", NA_character_, chebi_id),
    relations_relation = ifelse(relations_relation == "NA", NA_character_, relations_relation),
    relations_url = ifelse(relations_url == "NA", NA_character_, relations_url),
    database = ifelse(database == "NA", NA_character_, database),
    relations_term = ifelse(relations_term == "NA", NA_character_, relations_term)
  )

go_network_connections <- go_network %>%
  filter(child_id %in% metal_slim_subset_annotated$slims_from_id) %>%
  distinct(slims_from_id, child_id)

# Custom function to propagate chebi_ids
# Caution: this function cannot deal with two IDs converging on one unannotated note. It will take over
# the ID of the parent note with the least children.
# The only incorrect entry is: GO:0015277 has only sodium annotation even though there is also potassium

propagate_attribute <- function(g, attribute_name = NULL) {
  ids_not_na <- V(g)[!is.na(eval(bquote(`$`(V(g), .(as.name(attribute_name))))))]

  sub_ids_ids_not_na <- map(
    .x = ids_not_na,
    .f = ~ {
      igraph::subcomponent(g, .x, "out")
    }
  )

  sub_ids_ids_not_na

  deleted_ids <- list()

  for (i in 1:length(sub_ids_ids_not_na)) {
    current_id <- sub_ids_ids_not_na[[i]]
    current_id_length <- length(current_id)

    for (j in 1:length(sub_ids_ids_not_na)) {
      sub_id <- sub_ids_ids_not_na[[j]]
      if (current_id_length > length(sub_id)) {
        current_id <- current_id[!current_id %in% sub_id]
      }
    }

    deleted_ids[[i]] <- current_id
  }

  deleted_ids

  to_replace <- V(g)[is.na(eval(bquote(`$`(V(g), .(as.name(attribute_name))))))]

  for (i in 1:length(to_replace)) {
    print(to_replace[i])
    new_value <- NA
    for (j in 1:length(deleted_ids)) {
      if (to_replace[i] %in% deleted_ids[[j]]) {
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

propagated_chebi_ids <- igraph::as_data_frame(g_chebi_id_new, what = "vertices") %>%
  rename(slims_from_id = name)

############### relations_relation

go_network_vertex_attribute_relations_relation <- go_network %>%
  distinct(slims_from_id, relations_relation)

# Generate graph
g_relations_relation <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_relation)

# Update graph

g_relations_relation_new <- propagate_attribute(g = g_relations_relation, attribute_name = "relations_relation")

propagated_relations_relation <- igraph::as_data_frame(g_relations_relation_new, what = "vertices") %>%
  rename(slims_from_id = name)

############### relations_url

go_network_vertex_attribute_relations_url <- go_network %>%
  distinct(slims_from_id, relations_url)

# Generate graph
g_relations_url <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_url)

# Update graph

g_relations_url_new <- propagate_attribute(g = g_relations_url, attribute_name = "relations_url")

propagated_relations_urls <- igraph::as_data_frame(g_relations_url_new, what = "vertices") %>%
  rename(slims_from_id = name)

############### database

go_network_vertex_attribute_database <- go_network %>%
  distinct(slims_from_id, database)

# Generate graph
g_database <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_database)

# Update graph

g_database_new <- propagate_attribute(g = g_database, attribute_name = "database")

propagated_databases <- igraph::as_data_frame(g_database_new, what = "vertices") %>%
  rename(slims_from_id = name)

############### relations_term

go_network_vertex_attribute_relations_term <- go_network %>%
  distinct(slims_from_id, relations_term)

# Generate graph
g_relations_term <- igraph::graph_from_data_frame(go_network_connections, directed = TRUE, vertices = go_network_vertex_attribute_relations_term)

# Update graph

g_relations_term_new <- propagate_attribute(g = g_relations_term, attribute_name = "relations_term")

propagated_relations_terms <- igraph::as_data_frame(g_relations_term_new, what = "vertices") %>%
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
  mutate(
    chebi_id = str_split(chebi_id, pattern = ";"),
    relations_relation = str_split(relations_relation, pattern = ";"),
    relations_url = str_split(relations_url, pattern = ";"),
    database = str_split(database, pattern = ";"),
    relations_term = str_split(relations_term, pattern = ";")
  ) %>%
  unnest(c(chebi_id, relations_relation, relations_url, database, relations_term))

# These are all the GO terms with metal containing ChEBI IDs
metal_chebis <- metal_go_slim_subset_pre %>%
  mutate(chebi_id = ifelse(chebi_id == "CHEBI:39128", "CHEBI:18420", chebi_id)) %>% # just do this to correct an entry (magnesium ion to Mg2+)
  filter((chebi_id %in% unique(metal_chebi$chebi_accession)))

# These are all the GO terms with NA ChEBI IDs or potential metal ChEBI IDs
non_metal_chebis <- metal_go_slim_subset_pre %>%
  filter(!(slims_from_id %in% metal_chebis$slims_from_id)) # find all entries that are not in the metal_chebis list

message(
  "If this list is not empty please add to the more_metal_chebi_ids_wo_formula vector in line 273!\n In order for this not to be always empty first comment out the lines around 333.",
  '"', paste0(na.omit(unique(non_metal_chebis$relations_term)), collapse = '",\n"'), '"'
)

# Manual annotation of remaining GO IDs

message(
  paste0(
    paste0(
      'paste0("',
      unique(non_metal_chebis$slims_from_id), ";", unique(non_metal_chebis$main_name), '", ";", "CHEBI")'
    ),
    collapse = ",\n"
  )
)

manual_go_chebi_annotation <- data.frame(name = c(
  paste0("GO:0030946;protein tyrosine phosphatase activity, metal-dependent", ";", "CHEBI:25213"),
  paste0("GO:0061473;murein tripeptide carboxypeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0140492;metal-dependent deubiquitinase activity", ";", "CHEBI:25213"),
  paste0("GO:0000293;ferric-chelate reductase activity", ";", "CHEBI:29034|CHEBI:29033"),
  paste0("GO:0140487;metal ion sequestering activity", ";", "CHEBI:25213"),
  paste0("GO:0140486;zinc ion sequestering activity", ";", "CHEBI:29105"),
  paste0("GO:0008823;cupric reductase activity", ";", "CHEBI:49552|CHEBI:29036"),
  paste0("GO:0004222;metalloendopeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0070006;metalloaminopeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0015344;siderophore uptake transmembrane transporter activity", ";", "CHEBI:26672"),
  paste0("GO:0015343;siderophore transmembrane transporter activity", ";", "CHEBI:26672"),
  paste0("GO:0004023;alcohol dehydrogenase activity, metal ion-independent", ";", "CHEBI:25213"),
  paste0("GO:0051002;ligase activity, forming nitrogen-metal bonds", ";", "CHEBI:25213"),
  paste0("GO:0051003;ligase activity, forming nitrogen-metal bonds, forming coordination complexes", ";", "CHEBI:25213"),
  paste0("GO:0004181;metallocarboxypeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0008235;metalloexopeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0008237;metallopeptidase activity", ";", "CHEBI:25213"),
  paste0("GO:0043783;oxidoreductase activity, acting on metal ions, flavin as acceptor", ";", "CHEBI:25213"),
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
  paste0("GO:0106423;tubulin-tyrosine carboxypeptidase", ";", "CHEBI:25213"),
  paste0("GO:0015620;ferric-enterobactin transmembrane transporter activity", ";", "CHEBI:28199"),
  paste0("GO:1902945;metalloendopeptidase activity involved in amyloid precursor protein catabolic process", ";", "CHEBI:25213")
)) %>%
  separate(name, into = c("slims_from_id", "main_name", "chebi_id"), sep = ";") %>%
  mutate(chebi_id = str_split(chebi_id, pattern = "\\|")) %>%
  unnest(chebi_id) %>%
  left_join(distinct(non_metal_chebis, slims_from_id, slims_to_ids), by = "slims_from_id")

metal_go_slim_subset <- metal_chebis %>%
  bind_rows(manual_go_chebi_annotation) %>%
  distinct(slims_from_id, slims_to_ids, main_name, chebi_id, relations_relation, relations_term, database) %>%
  mutate(database = ifelse(is.na(database), "manual", database)) %>%
  left_join(distinct(metal_chebi, chebi_accession, metal_atom_id), by = c("chebi_id" = "chebi_accession")) %>%
  mutate(chebi_id = str_remove(chebi_id, pattern = "CHEBI:"))

usethis::use_data(metal_go_slim_subset, overwrite = TRUE)
