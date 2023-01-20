## Submission 

* We fixed a bug that caused an example to fail if there 
was a connectivity issue to a database. This should now be resolved.
* The bug caused an error to appear when \donttest examples were tested.
* This release therefore fixes the problem reported by Prof. Brian Ripley.

## Test environments
* macOS-latest (on GitHub actions), R 4.2.2
* windows-latest (on GitHub actions), R 4.2.2
* ubuntu-20.04 (on GitHub actions), R 4.2.2
* ubuntu-20.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* Ubuntu Linux 20.04.1 LTS (R-hub), r-release

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

