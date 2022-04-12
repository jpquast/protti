## Submission 

* We fixed the following issue raised by Prof. Brian Ripley. It was a simple bug in some of the fetch functions that has been fixed.

"Dear maintainer,

The donntest error showing is

checking examples with --run-donttest ... [101s/574s] ERROR
Running examples in ‘protti-Ex.R’ failed
The error most likely occurred in:

> ### Name: fetch_uniprot_proteome
> ### Title: Fetch proteome data from UniProt
> ### Aliases: fetch_uniprot_proteome
>
> ### ** Examples
>
> ## No test:
> head(fetch_uniprot_proteome(9606))
Error in `colnames<-`(`*tmp*`, value = column_names) :
attempt to set 'colnames' on an object with less than two dimensions
Calls: head -> fetch_uniprot_proteome -> colnames<-
Execution halted

Please correct before 2022-04-24 to safely retain your package on CRAN.

It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

This needs correction whether or not the resource recovers.

The CRAN Team"

## Test environments
* macOS-latest (on GitHub actions), R 4.1.3
* windows-latest (on GitHub actions), R 4.1.3
* ubuntu-20.04 (on GitHub actions), R 4.1.3
* ubuntu-20.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* Ubuntu Linux 20.04.1 LTS (R-hub), r-release

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

