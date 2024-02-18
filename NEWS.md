# protti 0.7.0.9000

## New features

* `calculate_go_enrichment()` got additional arguments.
  * `facet_n_col`: determines the number of columns the faceted plot should have if a group column is provided.
  * `plot_title`: specifies the title of the plot.
  * `min_n_detected_proteins_in_process`: argument for plotting that specifies the minimum number of proteins a GO term needs to be detected for.

# protti 0.7.0

## New features

* `correct_lip_for_abundance()` was added. It corrects LiP-peptides for changes in protein abundance and calculates their significance using a t-test. The function is based on the [MSstatsLiP](https://www.bioconductor.org/packages/release/bioc/html/MSstatsLiP.html) package developed by the Vitek Lab. Big thanks to @FehrAaron for implementing it! 
* `qc_cvs()` received a new argument called `max_cv` that specifies the maximum CV that should be included in the plot.
* `peptide_profile_plot()` received a new argument called `complete_sample`. If set to `TRUE`, each protein gets assigned all sample names that are found in the input data. This ensures that the plot always contains all samples on the x-axis even if there are no measured intensities for a specific sample. The default is `FALSE`, which is the original behaviour of the function.
* `volcano_plot()` received the `colour` argument that allows the user to provide custom colours for points.
* Increased the speed of `find_peptide()` and `assign_peptide_type()` by only computing on the smallest possible subset of data before joining back to the original data frame.
* `calculate_treatment_enrichment()` can now be applied on data frames with multiple different groups. The enrichment will be calculated for each group separately. If the data is plotted, each group is displayed in a separate facet. The group is provided to the new `group` argument.
* `qc_pca()`: If the condition argument is numeric a colour gradient is used instead.

## Bug fixes

* `volcano_plot()` now also works interactively if there are no significant hits.
* `fetch_chebi()`: fixed an issue caused by `na_if()` that changed its behaviour after the recent `dplyr` update.
* `qc_proteome_coverage()`: fixed the label order of fractions of proteins detected and not detected in the proteome. Fixes issue #194.
* `calculate_protein_abundance()` now correctly retains columns if `for_plot = TRUE`. Previously the columns to retain were not joined considering the precursor column, which lead to duplications of information where it did not belong. Fixes issue #197.
* `fetch_kegg()` now returns the pathway name correctly again.
* `qc_intensity_distribution()`, `qc_median_intensities()`, `qc_charge_states()`, `qc_contaminants()`, `qc_missed_cleavages()`, `qc_peptide_type()`, `qc_ids()`: If the provided sample column is of type factor, the level order won't be overwritten anymore.
*`fit_drc_4p()`: If there are no correlations an empty data frame is returned to prevent errors in `parallel_fit_drc_4p()`.
* `calculate_sequence_coverage()` does not fail anymore if a protein only contains `NA` peptide sequences.
* `qc_sequence_coverage()` does not return a plot anymore if `plot = FALSE`. This fixes issue #207.
* `qc_data_completeness()` if sample was of type `factor` the function did not properly facet the data when the `digestion` argument was provided. Now we filter out all 0% completeness values that come from factor levels that are not present in subsetted data. 

# protti 0.6.0

## New features

* `calculate_go_enrichment()` can now be applied on data frames with multiple different groups. The enrichment will be calculated for each group separately. If the data is plotted, each group is displayed in a separate facet. The group is provided to the new `group` argument. The `y_axis_free` argument determines if the y-axis of the faceted plot is "free" or "fixed".
* Added a `version` argument to `fetch_alphafold_prediction()` that specifies which verison of the database should be retrieved. The default is currently the newest version `"v4"`.
* `qc_ranked_intensities()` was added. It ranks protein, peptide or precursor intensities from highest to lowest. Ranked intensities can also be plotted using the `plot` argument.
* `fetch_chebi()` recieved a `timeout` argument that specifies after how many seconds the connection to the database should timeout. The default is 60 seconds as previously used.

## Bug fixes

* `pval_distribution_plot()` facets now have the correct style.
* `calculate_protein_abundance()` requires at least three distinct peptides for quantification. The function now applies this rule for each sample independently except for checking the whole dataset to contain at least three distinct peptides.

# protti 0.5.0

## New features

* `fetch_alphafold_aligned_error()` was added. It fetches the aligned error matrix for structure predictions from the AlphaFold EBI database.
* `predict_alphafold_domain()` was added. It uses a graph-based community clustering algorithm of AlphaFold predicted aligned errors in order to infer protein domains in AlphaFold predictions. The code is based on [python code](https://github.com/tristanic/pae_to_domains) by Tristan Croll.

## Bug fixes

* `assign_missingness()` now correctly deals with unequal replicate numbers of comparisons. In addition there is a message returned if an unequal number of replicates is detected for a comparison.
* `fetch_chebi()` fixed a bug that prevented the function from failing gracefully if there is a connection problem to the server.
* `extract_metal_binders()` now checks if the provided data frames are `NULL`. If yes, a message and `NULL` is returned.
* `fetch_mobidb()` was updated after the API changed.

## Additional changes

* Updated the "Protein Structure Analysis Workflow" vignette to include the `fetch_alphafold_aligned_error()` and `predict_alphafold_domain()` functions.

# protti 0.4.0

## New features

* Reintroduced the functionalities relying on the `iq` package to `protti`. `calculate_protein_abundance()` now has the method `"iq"` again as an option.
* `fetch_pdb()` now also retrieves information on engineered mutations, non-standard monomers, secondary structure and binding interfaces of ligands.
* `extract_metal_binders()` was completely redone. This was in response to the UniProt update and rework of the binding column provided by UniProt. This function extracts and concatenates all metal binding information available for a protein based on the UniProt and QuickGO databases. Therefore, this function now also takes gene ontology (GO) information from QuickGO as input. Instead of being able to provide column names to specific argument of the function you now only provide the data frames. This makes the function less flexible but reduces the amount of arguments required to achieve the same result. You just need to make sure that the input data frames contain columns with the correct names as stated by the function documentation.
* `fetch_quickgo()` was added. It fetches gene ontology (GO) information from the QuickGO EBI database. The retrieved information can either be GO annotations for provided UniProt IDs or Taxon identifiers, a list of all GO terms or a "slims" subset of GO IDs that can be generated based on provided GO IDs.
* `fetch_chebi()` now has the `stars` argument with which one can select the evidence levels for which entries should be retrieved.

## Bug fixes

* Fixed the `auth_seq_id` column that is part of the output of the `fetch_pdb()` function. Previously, the column could contain duplicated or missing positions. This was formerly identified by comparing the number of positions within the `auth_seq_id` column and the number of residues in the deposited `pdb_sequence`. Positions are now correct. The original output can be found in the `auth_seq_id_original` column.
* In the `calculate_diff_abundance()` function the intensity column can now be retained with the `retain_columns` argument. This was previously not possible until now since this column was used to reduce the annotation dataset. However, after reassessing the benefit of this filter step, it seemed not necessary. 
* We assumed that users would only retain columns in `calculate_diff_abundance()` that would not duplicate the data. However, this seems not to be the case, which can lead to wrong p-value adjustment. p-value adjustment was originally performed after the columns indicated in `retain_columns` are joined back to the data. Now p-value adjustment is performed prior to retaining columns as well as only on the subset of data that actually contains p-values. Previously we (by default `filter_NA_missingness = TRUE`) only filtered out `NA`s in the `missingness` column prior to p-value adjustment. However, it is possible that `missingness` is not `NA` but the p-value is `NA`. Now for all methods except for `"proDA"` we remove `NA` p-values before p-value adjustment. For `"proDA"` data is handled as previously since p-values are never `NA`. 

## Additional changes

* The default batchsize of `fetch_pdb()` was changed to 100 (from 200). This was done since more information is retrieved now, which slows to function down and is slightly improved when batch sizes are smaller.
* `try_query()` now only retries to retrieve information once if the returned message was "Timeout was reached". In addition, a `timeout` and `accept` argument have been added.
* The UniProt database has changed its API, therefore column names have changed as well as the format of data. We adjusted the `fetch_uniprot()` and `fetch_uniprot_proteome()` function accordingly. Please be aware that some columns names might have changed and your code might throw error messages if you did not adjust it accordingly. 

# protti 0.3.1

## Bug fixes

* Corrected the "Protein Structure Analysis Workflow" vignette. The example for `map_peptides_on_structure()` was still using "residue" as its input. We now use `find_peptide_in_structure()` to generate the correct input column.
* Fixed a bug in `fetch_uniprot()` and `fetch_uniprot_proteome()`. As UniProt has updated their website and their programmatic access, we now download the information from the legacy version temporarily. A real fix will follow.
* Fixed a bug in `fetch_kegg()`. The function did not retrieve any data after the API URL had changed.

# protti 0.3.0

## New features

* The "Protein Structure Analysis Workflow" vignette was added. It contains an example workflow for the analysis of structural proteomics data.
* `fetch_eco()` was added. It fetches evidence & conclusion ontology information from the EBI database.
* `qc_proteome_coverage()` now has the `reviewed` argument that specifies if only reviewed entries in UniProt should be considered as the proteome. The default is `TRUE` and stays the same as previously.
* `volcano_plot()` now has the `facet_scales` argument that specifies if the scales should be "free" or "fixed" when a faceted plot is created. The arguments that can be provided are the same that can be provided to the `scales` argument of `ggplot2::facet_wrap()`. The new default is now `"fixed"`.
* `pval_distribution_plot()` now has the optional `facet_by` variable that allows faceting of the plot.

## Bug fixes

* Fixed a bug in `map_peptides_on_structure()` that caused an error if the column provided to the `auth_seq_id` argument was called "residue".
* Fixed a bug in `volcano_plot()` that did not calculate the horizontal cutoff line correctly if there were multiple significance values that have the same adjusted significance value. Now it correctly uses the two p-values closest to the cutoff for the line position calculation. In addition, points were not correctly displayed if no horizontal cutoff line was created due to no significant values. Now all values are displayed correctly.
* Fixed a bug related to fetch functions not failing gracefully. The problem was that the internal `try_query()` function now returns errors as a character string if it encounters one. Functions using `try_query()` however, still expected `NULL` if there was an error. Also adjusted additional fetch functions that do not use `try_query()` to fail gracefully and to return informative messages upon encountering errors.

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
