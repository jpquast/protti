# protti 0.2.2.9000

## New features

* `fetch_eco` was added. It fetches evidence & conclusion ontology information from the EBI database.
* The "Protein Structure Analysis Workflow" vignette was added. It contains an example workflow for the analysis of structural differential proteomics data.

## Bug fixes

* Fixed a bug in `map_peptides_on_structure()` that caused an error if the column provided to the `auth_seq_id` argument was called "residue".

## Additional changes

* Improved test coverage for a few functions.

# protti 0.2.2

## New features

* `calculate_go_enrichment()` now has the argument `label`. If `TRUE` labels are added to the plot that specify how many proteins from the specific GO term are among the significant hits. This new argument is by default `TRUE`. 
* The version of STRINGdb that should be used for network analysis can now be provided through the `version` argument in the `analyse_functional_network()` function.
* Examples were added to some additional functions.

## Bug fixes 

* All tests and examples that are run on CRAN servers are only for functions that use packages `protti` depends on and not also packages it suggests.
* Removed the functionalities relying on the `iq` package from `protti` for now since `iq` is currently not available on CRAN. Once it is available again we will add the functionalities back.
* `fetch_metal_pdb()` now gives more informative feedback regarding the reasons resources were not fetched correctly.
* Fixed a bug in `qc_proteome_coverage` that flipped the "Detected" and "Not detected" labels.
* When the `highlight` argument of the `woods_plot()` function was used huge plots were generated when more than 20 proteins were provided to the function. This is fixed and was due to a wrong variable used in a loop.

# protti 0.2.1

## New features

* `create_structure_contact_map()` now takes an optional `data2` argument in which a second data frame of positions and structures can be provided. The positions in this data frame are used for distance calculations relative to the positions provided in the `data` argument. This helps to reduce the size of the map since only comparisons of interest are made. If a selection provided through the `data` argument should be compared to the whole structure the `data2` argument should not be provided, which is the previous default setting. 
* `parallel_create_structure_contact_map()` can create structure contact maps using parallel processing. It is also recommended for sequential processing if a large number of contact maps should be created. The non-parallel function should not be used if more than 50 maps should be created at the same time. In order to reduce contact maps to domains or even only peptides, a search pattern is created that selects only relevant regions. This pattern is shared over all maps and would become too large if too maps are created at once. 
* The "combined" condition in the plot of the `qc_cvs()` function now has a grey colour. This ensures that the other colours for each condition match the colours of other quality control plots.
* The `run_order` argument of `qc_sample_correlation()` now has a gradient colour scheme that is easier to interpret than distinct colours. The colours are inspired by the "plasma" colour scheme of the viridis package.
* If the sample column provided to `qc_intensity_distribution()` is of type factor instead of character the provided factor levels are used for sample ordering in the plot. This allows for custom sample ordering. If you want this functionality in other functions, please let us know in an issue on GitHub.

## Bug fixes and documentation updates

* Fixed a bug in `qc_ids()`, which caused an error when the optional `condition` argument was not provided. Also fixed a bug that did not take the state of the `remove_na_intensities` argument into account.
* Small documentation updates in the "Dose-Response Data Analysis Workflow" vignette in which we correct statements about multiple testing correction and p-value distributions.
* The `auth_seq_id` variable in structure files was handled as a numeric variable even though it sometimes can be character. This was fixed in all functions using it. These include: `fetch_pdb_structure()`, `fetch_alphafold_prediction()` (output is now numeric), `fetch_pdb()` (provides the `auth_seq_id` column not as a list vector but as character vector with semicolons as separators), `find_peptide_in_structure()` (now uses the new `auth_seq_id` format and returns the character vector of all positions of each peptide as `auth_seq_id` in additon to `auth_seq_id_start` and `auth_seq_id_end`), `create_structure_contact_map()`, `map_peptides_on_structure()` (now take the new `auth_seq_id` column instead of start and end positions as input).
* Fixed a bug in `map_peptides_on_structure()`, which caused an error when file names larger than 256 characters were generated. This was the case if the number of proteins that are part of one structure is very large (e.g. ribosome). If the number of proteins is too large to fit into a normal length file name, they are abbreviated as for example "51_proteins".
* `fetch_alphafold_prediction()` and `fetch_pdb_structure()` now return more informative messages if individual IDs have not been retrieved correctly (e.g. internet connection problem or outdated ID).
* When "proDA" was selected as method in `calculate_diff_abundance()`, comparison names starting with numbers were not supported. This has been fixed. There are no more restrictions for condition names.

# protti 0.2.0

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
* `volcano_plot()` can now plot unadjusted p-values while the horizontal cutoff line is based on the adjusted p-value. the `significance_cutoff` argument was modified to also take a second value which specifies the adjusted p-value column name. In that case the y-axis of the plot could display p-values that are provided to the `significance` argument, while the horizontal cutoff line is on the scale of adjusted p-values transformed to the scale of p-values. The provided vector can be e.g. `c(0.05, "adj_pval")`. In that case the function looks for the closest adjusted p-value above and below 0.05 and takes the mean of the corresponding p-values as the cutoff line. If there is no adjusted p-value in the data that is below 0.05 no line is displayed. This allows the user to display volcano plots using p-values while using adjusted p-values for the cutoff criteria. This is often preferred because adjusted p-values are related to unadjusted p-values often in a complex way that makes them hard to be interpret when plotted.

## Renamed functions
 
 Multiple functions have been renamed. More function follow the convention that they should be verbs. Other functions are more similar to each other in their naming. The old functions still work but they are deprecated and will be removed in the future. Please use the new versions instead.
 
* `diff_abundance()` has been renamed to `calculate_diff_abundance()`. 
* `peptide_type()` has been renamed to `assign_peptide_type()`. 
* `median_normalisation()` has been renamed to `normalise()`. The normalisation method is now defined through an argument called `method`.
* `sequence_coverage()` has been renamed to `calculate_sequence_coverage()`.
* `network_analysis()` has been renamed to `analyse_functional_network()`. 
* `treatment_enrichment()` has been renamed to `calculate_treatment_enrichment()`.
* `go_enrichment()` has been renamed to `calculate_go_enrichment()`. 
* `kegg_enrichment()` has been renamed to `calculate_kegg_enrichment()`.
* `volcano_protti()` has been renamed to `volcano_plot()`. 
* `plot_drc_4p()` has been renamed to `drc_4p_plot()`. 
* `plot_peptide_profiles()` has been renamed to `peptide_profile_plot()`.
* `plot_pval_distribution()` has been renamed to `pval_distribution_plot()`. 


## Bug fixes and documentation updates

* `calculate_diff_abundance()` had a bug that could not correctly deal with condition names containing spaces when a moderated t-test or the proDA algorithm was used. This has been fixed.
* `scale_protti()` can now deal with the case that a vector of only equal numbers (e.g. `c(1, 1, 1)`) is provided. In this case it always scales all values to 1 for `method = "01"` and to 0 for `method = "center"`.
* `pval_distribution_plot()` had a bug centered each bin around values, which meant that the first and last bin were smaller. This has been fixed.
* The `header` argument was removed from `try_query()`. It can now be directly supplied to the elipsis (...) of this function depending on which read method is used. For `type = "text/tab-separated-values"` the argument `col_names` of the corresponding `readr` function `read_tsv()` can be supplied.
* `split_metal_name()` removes more non-metal specific words that lead to false identifications. Also many metal names containing a "-" are now considered. A bug was resolved that did not deal with iron-sulfur cluster information correctly due to a wrong variable name.
* `extract_metal_binders()` now also correctly identifies monovalent inorganic cations as metals.
* Improve documentation of "Dose-Response Data Analysis Workflow" vignette. Added code for the creation of the `passed_filter` example column that is used in functions afterwards. 
* Examples for most functions have been added.
* Function documentation and code has been formatted to fit mostly within 100 characters line width.

# protti 0.1.1

* Fixed adjustment for multiple testing in the `diff_abundance()` function. Previously, the p-values in the whole result table were adjusted together. Now adjustments are made by comparisons. This change only affects analyses that have more than one comparison (two conditions). 
* Fixed `fetch_uniprot()`, `fetch_uniprot_proteome()` and `fetch_kegg()`. They now fail gracefully with informative messages if there is no internet connection, the database times out or the data resource has changed. Specifically the underlying `try_query()` function was changed.
* Small fixes in some examples.
* Small changes in Vignettes. 

# protti 0.1.0

* First release version.
