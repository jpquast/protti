## Submission 

This submission introduces new features and fixes a few bugs.

## Test environments
* macOS-latest (on GitHub actions), R 4.1.1
* windows-latest (on GitHub actions), R 4.1.1
* ubuntu-16.04 (on GitHub actions), R 4.1.1
* ubuntu-16.04 (on GitHub actions), r-devel
* windows-ix86+x86_64 (win-builder), r-devel
* fedora-clang-devel (R-hub), r-devel
* windows-x86_64-devel (R-hub), r-devel
* Ubuntu Linux 20.04.1 LTS (R-hub), r-release

## R CMD check results

> checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Jan-Philipp Quast <quast@imsb.biol.ethz.ch>'
  
Found the following (possibly) invalid URLs:
  URL: https://metalpdb.cerm.unifi.it
    From: man/fetch_metal_pdb.Rd
    Status: Error
    Message: SSL certificate problem: unable to get local issuer certificate
  URL: https://metalpdb.cerm.unifi.it/api_help
    From: man/fetch_metal_pdb.Rd
    Status: Error
    Message: SSL certificate problem: unable to get local issuer certificate
  URL: https://metalpdb.cerm.unifi.it/perGeometry
    From: man/fetch_metal_pdb.Rd
    Status: Error
    Message: SSL certificate problem: unable to get local issuer certificate
    
I get these notes. I can access the URL through the browser. I assume this is an issue specific to this database. If this is a problem with CRAN I can always remove the links, less convenient for the user but possible.
    
  URL: https://www.maxquant.org/
    From: inst/doc/data_analysis_dose_response_workflow.html
          inst/doc/input_preparation_workflow.html
          inst/doc/quality_control_workflow.html
    Status: Error
    Message: SSL certificate problem: certificate has expired

I would assume they will renew their certificate soon since it is the main website this widely used software is distributed through. Again if it is a problem with CRAN I can remove the link for now.


