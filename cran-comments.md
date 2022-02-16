## Submission 

This submission fixes some issues that were raised by Prof. Brian Ripley:

"Packages in Suggests should be used conditionally: see 'Writing R Extensions'.
This needs to be corrected even if the missing package(s) become available.
It can be tested by checking with _R_CHECK_DEPENDS_ONLY_=true."

We tested the package using the suggested environmental variable and fixed all examples 
and tests that caused problems.

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

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

