## Submission 

This submission fixes the issue raised by Prof. Brian Ripley (below). In addition small bug fixes to a few functions were made. 

* A macOS run overnight gave 
Running the tests in ‘tests/testthat.R’ failed.
where you are blaming the user for your package's error!

It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

Please correct before 2021-05-18 to safely retain your package on CRAN.  This needs correction whether or not the resource recovers.

	* This is now fixed. Functions fetch_uniprot, fetch_uniprot_proteome and fetch_kegg now fail with an informative message if there is no internet connection, the connection timed out or the data resource has changed. 

## Test environments
* macOS-latest (on GitHub actions), R 4.0.5
* windows-latest (on GitHub actions), R 4.0.5
* ubuntu-16.04 (on GitHub actions), R 4.0.5
* ubuntu-16.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* ubuntu-gcc-release (R-hub), r-release

## R CMD check results

> checking CRAN incoming feasibility ... NOTE  Maintainer: 'Jan-Philipp Quast <quast@imsb.biol.ethz.ch>'    Possibly mis-spelled words in DESCRIPTION:
  Feng (16:158)

0 errors | 0 warnings | 1 notes

This is a name not a mis-spelled word.

