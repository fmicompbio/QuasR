
countAlignments <- function(file, destination, maxHits=getOption("quasr.maxhits"), ...,
             overwrite=FALSE, indexDestination=TRUE)
{
    destination <- unlist(strsplit(file, "\\.bam$"))
    d0 <- paste(destination, "bam", sep=".")

    allowedMaxHits <- .Call(.getAllowedMaxHits)
    if(maxHits > allowedMaxHits) ## TODO error handling
      warning("Maximal allowed hits is '", allowedMaxHits, "'.")

    cntFile <- tempfile()
    sortFile <- tempfile()
    on.exit(unlink(c(cntFile, sortFile)))
    
    tryCatch({
        if (!overwrite && file.exists(d0)) {
            msg <-
                sprintf("'destination' exists, 'overwrite' is FALSE\n  destination.bam: %s",
                        "destination", "overwrite", d0)
            stop(msg)
        }
        sortFile <- sortBam(file, sortFile, byQname=TRUE)        
        cntFile <- .Call(.countAlignments, sortFile, cntFile, maxHits)
        if (!file.exists(cntFile))
            stop("failed to create 'BAM' file")
        if(indexDestination) {
            destination <- sortBam(cntFile, destination)
            indexBam(destination)
        } else {
            destination <- d0
            .file.rename(cntFile, destination)
        }
    }, error=function(err) {
        msg <- sprintf("'countAlignments' %s\n  SAM file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })
    destination
}
