###
###
###

## Ensure that the global package options are set to their default
.onAttach <- function(libpath, pkg) {
    options(quasr.countApp="R")
    if(suppressWarnings(require(parallel, quietly=TRUE)))
        options(quasr.cores=max(detectCores()-2,1))
    else
        options(quasr.cores=1)
  options(quasr.maxmem=2048)
}


##.onLoad <- function(libname, pkgname)

## TODO check if detach of dyn library is done

