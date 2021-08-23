# protti 0.2.0.9000

## New features

* `fetch_pdb()` was added. It fetches PDB structure metadata from RCSB.
* `fetch_pdb_structure()` was added. It fetches atom level data for a PDB structure from RCSB.
* `fetch_metal_pdb()` was added. It fetches information about protein-metal binding sites from the MetalPDB database.
* `find_peptide_in_structure()` was added. It returns the positions of a protein region, peptide or amino acid within a protein structure using UniProt positions as input. This is necessary because often amino acid positions in a protein structure vary from their positions in the UniProt protein sequence. 
* `map_peptides_on_structure()` was added. It can map peptides onto PDB structures based on their positions. This is accomplished by replacing the B-factor information in the structure file with values that allow highlighting peptides when the structure is coloured by B-factor.
* `create_structure_contact_map()` was added. It creates a contact map of a subset or of all atom or residue distances in a structure file.
* `calculate_aa_scores()` was added. It calculates a score for each amino acid of a protein based on the product of the -log10(adjusted p-value) and the absolute log2 fold change per peptide.
* `woods_plot()` recieved the `highlight` argument, which allows to highlight peptides in the plot with an asterisk based on a logical column. It also got `export` and `export_name` arguments which now makes it possible to directly export plots from the function. The new `target` argument allows the user to specify one or multiple proteins from their data frame for which a plot should be returned, however the default option is to return plots for all proteins in the provided data frame. Now it is also possible to provide more than 20 proteins at a time while the `facet` argument is `TRUE`. For every 20 proteins a new plot is created and all of them are returned together in a list. 
* `assign_missingness()` and `calculate_diff_abundance()` now also take `"all"` as input to their `ref_condition` argument. This will create all pairwise condition pairs. Previously only one reference condition could be provided.

## Renamed functions
 
 Multiple functions have been renamed. More function follow the convention that they should be verbs. Other functions are more similar to each other in their naming. The old functions still work but they are deprecated and will be removed in the future. Please use the new versions instead.
 
* `diff_abundance()` has been renamed to `calculate_diff_abundance()`. 
* `network_analysis()` has been renamed to `create_functional_network()`. 

## Bug fixes and documentation updates

* `scale_protti()` can now deal with the case that a vector of only equal numbers (e.g. `c(1, 1, 1)`) is provided. In this case it always scales all values to 1 for `method = "01"` and to 0 for `method = "center"`.
* The `header` argument was removed from `try_query()`. It can now be directly supplied to the elipsis (...) of this function depending on which read method is used. For `type = "text/tab-separated-values"` the argument `col_names` of the corresponding `readr` function `read_tsv()` can be supplied.
* `split_metal_name()` removes more non-metal specific words that lead to false identifications. Also many metal names containing a "-" are now considered. A bug was resolved that did not deal with iron-sulfur cluster information correctly due to a wrong variable name.
* `extract_metal_binders()` now also correctly identifies monovalent inorganic cations as metals.
* Improve documentation of "Dose-Response Data Analysis Workflow" vignette. Added code for the creation of the `passed_filter` example column that is used in functions afterwards. 
* Small fixes in some examples.

# protti 0.1.1

* Fixed adjustment for multiple testing in the `diff_abundance()` function. Previously, the p-values in the whole result table were adjusted together. Now adjustments are made by comparisons. This change only affects analyses that have more than one comparison (two conditions). 
* Fixed `fetch_uniprot()`, `fetch_uniprot_proteome()` and `fetch_kegg()`. They now fail gracefully with informative messages if there is no internet connection, the database times out or the data resource has changed. Specifically the underlying `try_query()` function was changed.
* Small fixes in some examples.
* Small changes in Vignettes. 

# protti 0.1.0

* First release version.
