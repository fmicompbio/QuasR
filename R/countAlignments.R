
.countAlignments <- function(file, grange, stranded=FALSE, overlap=c("any", "within"))
{
    ## check arguments, TODO more
    overlap <- match.arg(overlap)

    ## translate seqname to tid
    seqnamemap <- .Call(.seqname, file[1])
    regions <- data.frame(tid=seqnamemap$tid[ IRanges::as.vector(IRanges::match(seqnames(grange), seqnamemap$seqnames)) ],
                          start=as.integer(start(grange)),
                          end=as.integer(end(grange)),
                          strand=as.character(strand(grange)),
                          stringsAsFactors=FALSE
                          )
    ## get counts
    count <- .Call(.count_alignments, file, file, regions, stranded, overlap)

    return(count)
}
