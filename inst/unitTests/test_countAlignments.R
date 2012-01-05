test_countAlignments <- function()
{    
    samFilename <- system.file(package="QuasR", "unitTests", "case", "test.sam")
    bamFilename <- asBam(samFilename, tempfile(), indexDestination=FALSE, overwrite=TRUE)
    on.exit(unlink(bamFilename))
    
    ## get mapped Alignments
    aln <- readBamGappedAlignments(bamFilename, param=ScanBamParam(what="qname"))
    gr <- as(aln, "GRanges")
        
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", shift=0L, minoverlap=1L, maxHits=NULL)
    
    DEACTIVATED("Not implemented yet")
#     QuasR:::.countAlignments(bamFilename, gr, stranded=TRUE,
#                              overlap=c("any", "within", "startwithin", "endwithin"),
#                              shift=0L, minoverlap=1L, maxHits=NULL)
}