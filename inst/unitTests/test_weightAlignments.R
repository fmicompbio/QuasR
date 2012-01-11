test_weightAlignments <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    on.exit(unlink(system.file(package="QuasR", "extdata", "phage_genomes", "RbowtieIndex_phage_genomes"),
                   recursive=TRUE), add=TRUE)
        
    samFilename <- system.file(package="QuasR", "unitTests", "case", "test.sam")
    bamFilename <- asBam(samFilename, tempfile(), indexDestination=FALSE, overwrite=TRUE)
    on.exit(unlink(bamFilename), add=TRUE)
    
    ## get mapped Alignments
    aln1 <- readBamGappedAlignments(bamFilename, param=ScanBamParam(what="qname"))
    count <- tapply(values(aln1)$qname, values(aln1)$qname, length)
    
    project <- qProject(sampleFile, genomeName, bamfileDir=td, maxHits=90L)
    
    allowedMaxHits <- .Call(QuasR:::.getAllowedMaxHits)
    checkEqualsNumeric(32766, allowedMaxHits)
    
    hitcountDist <- QuasR:::.weightAlignments(bamFilename, bamFilename, 
                                              index=project@env$index$genome,
                                              qproject=project, 
                                              overwrite=TRUE, indexDestination=FALSE)
    cat(hitcountDist)
    
    ## check if same filename
    ##checkIdentical(bamFilename, newbamFilename)
    
    ## check hitcount distribution
    # TODO
    
    ## check if there is an error if bamFile exists already
    checkException(QuasR:::.weightAlignments(bamFilename, bamFilename,
                                             index=project@env$index$genome,
                                             qproject=project,
                                             overwrite=FALSE))
    
    ## check if there is a warning for to big maxhits
    #QuasR:::.weightAlignments(bamFilename, tempfile(), maxHits=allowedMaxHits+1, overwrite=TRUE)
    
    ## get Alignments with weights
    param <- ScanBamParam(what=c("qname","rname"), tag="IH") 
    list <- as.data.frame(scanBam(bamFilename, param=param)[[1]], stringsAsFactors=FALSE)
    ## check weights
    mapped <- list[list$IH > 0,]
    lapply(1:nrow(mapped), function(row){ 
        checkEqualsNumeric(count[mapped[row,]$qname], mapped[row,]$IH)
        })       
    notmapped <- list[list$IH <= 0,]
    lapply(1:nrow(notmapped), function(row){
        checkEqualsNumeric(0, notmapped[row,]$IH)
        checkTrue(is.na(count[notmapped[row,]$qname]))
        })
    ## check max hits
    maxHits <- max(count)-1 # TODO
    project <- qProject(sampleFile, genomeName, bamfileDir=td, maxHits=maxHits)
    # checkException() TODO checkWarning()
    hitcountDist2 <- suppressWarnings(QuasR:::.weightAlignments(bamFilename, bamFilename, 
                                               index=project@env$index$genome,
                                               qproject=project,
                                               overwrite=TRUE, indexDestination=FALSE))
    cat(hitcountDist2)
    ## check if equel hitcounts
    list <- as.data.frame(scanBam(bamFilename, param=param)[[1]], stringsAsFactors=FALSE)

    mapped <- list[list$IH > 0 & list$IH <= maxHits,]
    lapply(1:nrow(mapped), function(row){ 
        checkEqualsNumeric(count[mapped[row,]$qname], mapped[row,]$IH, 
                           msg=sprintf("\nTest mapped qname '%s': Not equal numeric  target '%i' current '%i'",
                                       mapped[row,]$qname, count[mapped[row,]$qname], mapped[row,]$IH))
        })
    
    notmapped <- list[list$IH <= 0,]
    lapply(1:nrow(notmapped), function(row){
        checkEqualsNumeric(0, notmapped[row,]$IH, 
                           msg=sprintf("\nTest notmapped qname '%s': Not equal numeric  target '%i' current '%i'",
                                       notmapped[row,]$qname, 0, notmapped[row,]$IH))
        checkTrue(is.na(count[notmapped[row,]$qname]))
        })    
    
    overmapped <- list[list$IH > maxHits,]
    lapply(1:nrow(overmapped), function(row){ 
        checkEqualsNumeric(allowedMaxHits, overmapped[row,]$IH, 
                           msg=sprintf("\nTest overmapped qname '%s': Not equal numeric  target '%i' current '%i'", 
                                       overmapped[row,]$qname, allowedMaxHits, overmapped[row,]$IH))
        })
}

# test_createSamHeaderFile <- function()
# {
#     hf <- tempfile(tmpdir=td, fileext=".sam")
#     cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
#         sprintf("at:%s", genomeName),
#         "CL:\"Add or modifiy IH Tag\"", sep="\t", file=hf)
#     bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
#     listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
#     checkException(QuasR:::.getBamFile(listBamFilenames, genomeName, ""))
# 
#     ## missing genomeName in bamheader
#     ## uncertain situation the process should not continue
#     hf <- tempfile(tmpdir=td, fileext=".sam")
#     cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
#         sprintf("at:%s", genomeName),
#         "CL:\"Add or modifiy IH Tag\"", sep="\t", file=hf)
#     bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
#     listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
# }