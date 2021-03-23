## Resubmission 2
This is the second resubmission responding to the feedback of Uwe Ligges on first resubmission below:

* Overall checktime 24 min > 10 min. Pls reduce to less than 10 min.
	* Fixed the test time as requested. Non-essential functionalities are now only checked locally. New check time in seconds: 426

## Resubmission 1
This is a resubmission. In this version I have:

* Fixed the (possibly) invalid URLs.
* Included the checks from R-hub and win-builder that actually show me the notes I get from CRAN.
* Not fixed the note of possibly mis-spelled words in the DESCRIPTION file. They are not mis-spelled. 

## Test environments
* macOS-latest (on GitHub actions), R 4.0.4
* windows-latest (on GitHub actions), R 4.0.4
* ubuntu-16.04 (on GitHub actions), R 4.0.4
* ubuntu-16.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* ubuntu-gcc-release (R-hub), r-release

## R CMD check results

> checking CRAN incoming feasibility ... NOTE  Maintainer: 'Jan-Philipp Quast <quast@imsb.biol.ethz.ch>'    Possibly mis-spelled words in DESCRIPTION:    LiP (2:33, 16:141)    Proteomics (2:18)    proteolysis (16:110)    proteomics (16:49, 16:171)    workflows (16:35)

0 errors | 0 warnings | 1 notes

* This is a new release.
