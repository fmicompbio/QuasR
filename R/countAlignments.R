
.countAlignments <- function(file, grange)
{
    ## TODO check if single filename
    names <- .Call(.seqname, file)
    refseq <- as.factor(as.character(seqnames(grange)))
    ## HACK translate seqname to tid. Is there a cleaner version?
    ##names$tid[names$seqnames=="chr4"]
    tidlevels <- unlist(lapply(levels(refseq), function(level){
            as.character(names$tid[names$seqnames==level])
        }))
    levels(refseq) <- tidlevels

    ## TODO tid, start, end as a single IRanges object instead of three integer object
    tid <- as.integer(as.character(refseq))
    start <- as.integer(start(ranges(grange)))
    end <- as.integer(end(ranges(grange)))

    ## TODO cleanup the index file handling
    count <- unlist(lapply(1:(length(tid)), function(i, tid, start, end){
            .Call(.count_alignments, file, file, tid[i], start[i], end[i])
        }, tid, start, end))

    return(count)
}