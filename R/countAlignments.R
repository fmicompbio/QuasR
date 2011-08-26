
.countAlignments <- function(file, grange, type=c("any"))
{
    ## TODO check if single filename
    
    seqnamemap <- .Call(.seqname, file)
    ## translate seqname to tid.
    regions <- data.frame(tid=seqnamemap$tid[ IRanges::as.vector(IRanges::match(seqnames(grange), seqnamemap$seqnames)) ],
                          start=as.integer(start(ranges(grange))),
                          end=as.integer(end(ranges(grange)))
                          )
    
    count <- .Call(.count_alignments, file, file, regions, type)
    
    return(count)
}
