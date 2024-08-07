## Submission 

* We specifically addressed and fixed the issue raised by Prof. Brian Ripley:
  * The `analyse_functional_network()` function did not fail gracefully.
  * We implemented a `try_catch()` that specifically rescues the cases in which the `STRINGdb` package does not fail gracefully. This fixes the issue.
* Additionally we added new features and fixed bugs.

## Test environments
* macOS-latest (on GitHub actions), R 4.4.1
* windows-latest (on GitHub actions), R 4.4.1
* ubuntu-20.04 (on GitHub actions), R 4.4.1
* ubuntu-20.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* Ubuntu Linux 20.04.1 LTS (R-hub), r-release

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

