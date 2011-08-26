
.weightAlignments <- function(file, destination, maxHits=NULL, ..., overwrite=FALSE, indexDestination=TRUE)
{
    # TODO check length(file) == 1L
    destination <- unlist(strsplit(file, "\\.bam$"))
    d0 <- paste(destination, "bam", sep=".")

    allowedMaxHits <- .Call(.getAllowedMaxHits)
    if(is.null(maxHits) || maxHits > allowedMaxHits){
        maxHits <- allowedMaxHits
        warning("Maximal hits is set to '", allowedMaxHits, " maximal allowed hits'.")
    }

    cntFile <- tempfile()
    on.exit(unlink(cntFile))

    tryCatch({
        if (!overwrite && file.exists(d0)) {
            msg <-
                sprintf("'destination' exists, 'overwrite' is FALSE\n  destination.bam: %s", "destination", "overwrite", d0)
            stop(msg)
        }
        ## sort bamfile by shortread name      
        if(any(scanBamHeader(file)[[1]]$text$`@HD` == "SO:queryname")){
            sortFile <- file
        } else {
            sortFile <- tempfile()
            on.exit(unlink(sortFile))
            sortFile <- sortBam(file, sortFile, byQname=TRUE)
        }
             
        cntFile <- .Call(.weight_alignments, sortFile, cntFile, maxHits)
        if (!file.exists(cntFile))
            stop("failed to create 'BAM' file")
        if(indexDestination) {
            destination <- sortBam(cntFile, destination)
            indexBam(destination)
        } else {
            destination <- d0
            file.rename(cntFile, destination)
        }
    }, error=function(err) {
        msg <- sprintf("'countAlignments' %s\n  SAM file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })
    return(destination)
}
