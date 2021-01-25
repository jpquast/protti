.onAttach <- function(libname, pkgname) {
  if(.Platform$OS.type == "unix"){
  packageStartupMessage("\U1F469\U1F3FD\U200D\U1F52C Welcome to protti version ", utils::packageVersion("protti"), "! \U1F468\U1F3FC\U200D\U1F4BB 
                            \n\U1F52C Have fun analysing your data! \U1F4BB")
  }
  if(.Platform$OS.type == "windows"){
    packageStartupMessage("Welcome to protti version ", utils::packageVersion("protti"), "! 
                            \nHave fun analysing your data!")
  }
}