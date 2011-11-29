
.weightAlignments <- function(file, destination, index, qproject, overwrite=FALSE, indexDestination=TRUE)
{
    # TODO check length(file) == 1L
    destination <- unlist(strsplit(file, "\\.bam$"))
    d0 <- paste(destination, "bam", sep=".")

    maxHits <- qproject@env$maxHits
    allowedMaxHits <- .Call(.getAllowedMaxHits)
    if(is.null(qproject@env$maxHits) || qproject@env$maxHits > allowedMaxHits){
        maxHits <- allowedMaxHits
        warning("Maximal hits is set to '", allowedMaxHits, " maximal allowed hits'.")
    }

    cntFile <- tempfile()
    on.exit(unlink(cntFile))

    if (!overwrite && file.exists(d0)) {                
        stop(sprintf("Destination '%s' exists and parameter 'overwrite' is set to %s", d0, overwrite))
    }
    
    tryCatch({
        ## sort bamfile by shortread name
        if(any(scanBamHeader(file)[[1]]$text["@HD"] == "SO:queryname")){ ## TODO should be all not any
            sortFile <- file
        } else {
            sortFile <- tempfile()
            on.exit(unlink(sortFile), add=T)
            sortFile <- sortBam(file, sortFile, byQname=TRUE)
        }
        
        ## create sam header file
        headerFile <- .createSamHeaderFile(sortFile, index=index, qproject=qproject)
        on.exit(unlink(headerFile), add=T)

        hitcount <- .Call(.weight_alignments, sortFile, cntFile, headerFile, maxHits)
        names(hitcount) <- c(0:maxHits, paste(">", maxHits,sep=""))

        if(!file.exists(cntFile))
            stop("failed to create 'BAM' file")
        
        if(indexDestination) {
            destination <- sortBam(cntFile, destination)
            indexBam(destination)
        } else {
            destination <- d0
            file.rename(cntFile, destination)
        }
    }, error=function(err) {
        msg <- sprintf("'weightAlignments' %s\n  SAM file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })

    return(hitcount)
}

.createSamHeaderFile <- function(bamfile, headerFile, index, qproject, hitcount){
    ## read current bam_header
    h <- scanBamHeader(bamfile)[[1]]$text 
    pgLineIdx <- names(h) == "@PG"
    ppValues <- unlist(strsplit(h[pgLineIdx][[1]][grep("ID:", h[pgLineIdx])],":"))[2]
    h <- paste(names(h), unlist(lapply(h, paste, collapse="\t")), sep="\t")
    ## create new pg tag for quasr
    quasrVersion <- installed.packages()['QuasR', 'Version']
    # TODO alignmentParameter is never null
    alignmentParameter <- ifelse(is.null(qproject@env$alignmentParameter), # TODO alignmentParameter is never null
                                 "", 
                                 sprintf("\tap:%s", paste(qproject@env$alignmentParameter, collapse=",")))
    alignmentTarget <- sprintf("\tat:%s", index$name)
    #alignmentTarget <- sprintf("\tat:%s-%s", index$name, index$md5sum)
    hitcount <- ifelse(!missing(hitcount),
                        sprintf("\tml:%s\tmc:%s", 
                                paste(names(hitcount), collapse=";"), 
                                paste(hitcount, collapse=";"))
                        , "")
    pgLine <- paste("@PG\tID:QuasR\tPN:QuasR\tVN:",
                    quasrVersion,
                    alignmentTarget,
                    alignmentParameter,
                    hitcount,
                    "\tCL:\"Add or modifiy IH Tag\"\tPP:",
                    ppValues, sep="")

    ## check if pg tag for quasr exists already
    quasrLineIdx <- grep("ID:QuasR", h[pgLineIdx])
    ## add new pg tag to bam_header
    if(length(quasrLineIdx) && quasrLineIdx == 1)
        h <- c(h[!pgLineIdx], pgLine, h[pgLineIdx][-1]) ## replase the quasr pg tag
    else 
        h <- c(h[!pgLineIdx], pgLine, h[pgLineIdx])
    if(missing(headerFile))
        headerFile <- tempfile(pattern="header_", fileext=".sam")
    ## create sam_header file
    cat(paste(h, collapse="\n"), file=headerFile)
    return(headerFile)
}