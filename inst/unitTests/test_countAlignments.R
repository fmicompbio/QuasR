test_countAlignments <- function()
{    
    samFilename <- system.file(package="QuasR", "unitTests", "case", "test.sam")
    bamFilename <- asBam(samFilename, tempfile(), overwrite=TRUE)
    on.exit(unlink(c(bamFilename, paste(bamFilename, ".bai", sep=""))))

    ## get mapped Alignments
    aln <- readBamGappedAlignments(bamFilename, param=ScanBamParam(what="qname"))
    readTable <- table(values(aln)$qname)
    gr <- sort(as(aln, "GRanges"))
    checkEquals(gr, reduce(gr))

    ## shift 0
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", shift=0L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="within", shift=0L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="startwithin", shift=0L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="endwithin", shift=0L)
    checkEquals(length(readTable), sum(ans))

    ## shift +1
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", shift=1L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="within", shift=1L)
    checkEquals(0, sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="startwithin", shift=1L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="endwithin", shift=1L)
    checkEquals(0, sum(ans))
    
    ## shift -1
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", shift=-1L)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="within", shift=-1L)
    checkEquals(0, sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="startwithin", shift=-1L)
    checkEquals(0, sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="endwithin", shift=-1L)
    checkEquals(length(readTable), sum(ans))
    
    ## stranded
    grAntiSens <- gr
    strand(grAntiSens) <- ifelse(strand(gr) == "-", "+", "-")
    ans <- QuasR:::.countAlignmentsR(bamFilename, grAntiSens, stranded=FALSE, overlap="any")
    checkEquals(length(readTable), sum(ans))
    
    ans <- QuasR:::.countAlignmentsR(bamFilename, grAntiSens, stranded=TRUE, overlap="any")
    checkEquals(0, sum(ans))
    
    ## minoverlap
    mo <- max(width(gr))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", minoverlap=mo)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", minoverlap=mo+1)
    checkEquals(0, sum(ans))

    ## maxHits TODO
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", maxHits=100)
    checkEquals(length(readTable), sum(ans))
    ans <- QuasR:::.countAlignmentsR(bamFilename, gr, stranded=TRUE, overlap="any", maxHits=1)
    checkEquals(sum(readTable == 1), sum(ans))
}