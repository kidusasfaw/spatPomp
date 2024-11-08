## set up default verbosity option

.onAttach <- function(libname,pkgname){
  options(spatPomp_verbose=FALSE)
  packageStartupMessage("kidusasfaw/spatPomp is no longer under development.\nUse spatPomp-org/spatPomp for the latest GitHub updates of spatPomp.")
}

