# protti 0.1.1

* Fixed adjustment for multiple testing in the `diff_abundance()` function. Previously, the p-values in the whole result table were adjusted together. Now adjustments are made by comparisons. This change only affects analyses that have more than one comparison (two conditions). 
* Fixed `fetch_uniprot()`, `fetch_uniprot_proteome()` and `fetch_kegg()`. They now fail gracefully with informative messages if there is no internet connection, the database times out or the data resource has changed. Specifically the underlying `try_query()` function was changed.
* Small fixes in some examples.
* Small changes in Vignettes. 

# protti 0.1.0

* First release version.
