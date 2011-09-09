
.countAlignments <- function(file, grange, type=c("any"))
{
    ## TODO check if single filename
#     if(length(file) != 1L)
#        stop("collapse TODO")

    ## translate seqname to tid
    seqnamemap <- .Call(.seqname, file[1])
    regions <- data.frame(tid=seqnamemap$tid[ IRanges::as.vector(IRanges::match(seqnames(grange), seqnamemap$seqnames)) ],
                          start=as.integer(start(ranges(grange))),
                          end=as.integer(end(ranges(grange)))
                          )
    ## get counts
    #cat(file, class(file), mode(file), sep="-")
    count <- .Call(.count_alignments, file, file, regions, type)

    return(count)
}
