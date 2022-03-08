## Submission 

The issues for which the package was archived have been fixed: We tested the package using the suggested environmental variable (_R_CHECK_DEPENDS_ONLY_=true in addition to --run-donttest) and fixed all examples 
and tests that caused problems.

We also addressed the comments from Gregor Seyer:

"Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 12 .Rd files, e.g.:
    diff_abundance.Rd: \value
    go_enrichment.Rd: \value
    kegg_enrichment.Rd: \value
    median_normalisation.Rd: \value
    network_analysis.Rd: \value
    peptide_type.Rd: \value
    ..."
    
* The value field has been added to these .Rd files.

"\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}."

* We removed the \dontrun{} wrapper from as much functions as possible and wrote working examples. 

There is one Note regarding misspelled words in the DESCRIPTION file. Those words are not misspelled.

## Test environments
* macOS-latest (on GitHub actions), R 4.1.2
* windows-latest (on GitHub actions), R 4.1.2
* ubuntu-20.04 (on GitHub actions), R 4.1.2
* ubuntu-20.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* Ubuntu Linux 20.04.1 LTS (R-hub), r-release

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 1 notes ✓

