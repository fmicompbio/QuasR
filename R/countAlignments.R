
.countAlignments <- function(file, grange)
{
    names <- .Call(.seqname, file)
    refseq <- as.factor(as.character(seqnames(grange)))
    ##names$tid[names$seqnames=="chr4"]
    tidlevels <- unlist(lapply(levels(refseq), function(level){
            as.character(names$tid[names$seqnames==level])
        }))
    levels(refseq) <- tidlevels

    tid <- as.integer(as.character(refseq))
    start <- as.integer(start(ranges(grange)))
    end <- as.integer(end(ranges(grange)))

    count <- unlist(lapply(1:(length(tid)), function(i, tid, start, end){
            .Call(.count_alignments, file, file, tid[i], start[i], end[i])
        }, tid, start, end))

    return(count)
}