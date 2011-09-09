###
###
###

## Ensure that the global package options are set to their default
.onAttach <- function(libpath, pkg) {
  options(quasr.cores=1)
  options(quasr.maxmem=2048)
}


##.onLoad <- function(libname, pkgname)

