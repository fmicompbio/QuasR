###
###
###

## Ensure that the global package options are set to their default
.onAttach <- function(libpath, pkg) {
    options(quasr.countApp="R")
    options(quasr.cores=max(detectCores()-2,1))
    options(quasr.clusterSize=2)
    options(quasr.maxmem=2048)
    options(quasr.quiet=FALSE)
}
