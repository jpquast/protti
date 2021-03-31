## Resubmission 4
This is the fourth resubmission responding to the feedback of Gregor Seyer on third resubmission below:

* You have examples for unexported functions.
  calculate_imputation() in:
     calculate_imputation.Rd
  drc_4p() in:
     drc_4p.Rd
  find_all_subs() in:
     find_all_subs.Rd
  find_chebis() in:
     find_chebis.Rd
  replace_identified_by_x() in:
     replace_identified_by_x.Rd
  split_metal_name() in:
     split_metal_name.Rd
  try_query() in:
     try_query.Rd
Please either omit these examples or export the functions.

	* Done. Removed examples.

* \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
Does not seem necessary.

Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.

	* Unwrapped functions with examples executable by the user. Many functions require an extensive data frame that cannot be easily provided so the example is not executable. These examples should mainly exemplify how to provide function arguments. The following functions now run their example code: create_queue, create_synthetic_data, find_peptide, median_normalisation, peptide_type, scale_protti, sequence_coverage, ttest_protti. The following functions were moved from \dontrun{} to \donttest{}: fetch_chebi, fetch_go, fetch_kegg, fetch_mobidb, fetch_uniprot_proteome, fetch_uniprot.

* Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies.
In your examples/vignettes/tests you can write to tempdir().

	* Functions that were writing by default to the working directory (create_queue) do now not export anything by default anymore. The user needs to specify that they want to export data by setting export to TRUE. Vignettes, examples and tests should not export anything. 

* Please ensure that you do not use more than 2 cores in your examples, vignettes, etc.

	* Done. Should not use multiple cores anymore. Removed read_protti function from CRAN tests.

I hope these changes fix the issues. If there is anything else related to your suggestions that I missed please let me know, I am happy to implement any further suggestions. Thanks a lot for your time and reviews!

## Resubmission 3
This is the third resubmission responding to the feedback of Uwe Ligges on second resubmission below:

* Found the following (possibly) invalid URLs:
   URL: https://www.thermofisher.com/de/de/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html (moved to https://www.thermofisher.com/at/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html)
     From: inst/doc/input_preparation_workflow.html
     Status: 200
     Message: OK

Please change http --> https, add trailing slashes, or follow moved content as appropriate.

	* The URL changes depending on the country it is accessed from. We changed it to "de/de" since this prevented the CRAN note at first due to the test servers being in Germany. However, it seems that there are also severs in Austria as it says that the URL has moved to "at/en". So if this is tested on servers in two different countries one of them will always give this note. I would suggest leaving this URL as it is. If you have any other suggestions of how to deal with this differently I am happy to implement them. 

* Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?

	* Done

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
