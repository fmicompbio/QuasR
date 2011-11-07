.countAlignments <- function(bamfile, grange, stranded=FALSE,
                             overlap=c("any", "within", "startwithin", "endwithin"),
                             shift=0L, minoverlap=1L, maxHits=NULL)
{

    ## check arguments
    if(!is.character(bamfile))
        stop("The parameter 'bamfile' must be of type character.")
    if(!is(grange, "GRanges"))
        stop("The parameter 'grange' must be of type GRanges.")
    if(!is.logical(stranded))
        stop("The parameter 'stranded' must be of type logical.")
    overlap <- match.arg(overlap)
    if(length(shift) != 1L && length(shift) != length(bamfile))
        stop("The parameter 'shift' must be of length one or of equal length as the parameter 'file'.")
    shift <- as.integer(shift)
    minoverlap <- as.integer(minoverlap)
    if(!is.null(maxHits))
        maxHits <- as.integer(maxHits)

    bamview <- BamViews(bamPaths=bamfile)

    count <- unlist(lapply(split(grange), function(range){
            
        ## TODO reduce resize(range,)
        nameRange <- names(range)
        if(shift > 0){
            expandedRange <- GRanges(ranges=IRanges(start=start(range)-max(shift), end=end(range)+max(shift)),
                                seqnames=seqnames(range), strand=strand(range))
            bamRanges(bamview) <- expandedRange
            names(range) <- sprintf("%s:%s-%s", seqnames(expandedRange), start(expandedRange), end(expandedRange))
        } else {
            bamRanges(bamview) <- range
            names(range) <- sprintf("%s:%s-%s", seqnames(range), start(range), end(range))
        }
        
        param <- ScanBamParam(what=c("qname","pos","cigar", "strand", "rname"),
                              tag="IH",
                              flag=scanBamFlag(isUnmappedQuery=FALSE))
        scanBamRes <- suppressWarnings(scanBam(bamview, param=param))
        reads <- GRangesList(lapply(names(scanBamRes[[1]]), function(region){
            reads <- GRangesList(lapply(scanBamRes, function(bf){
                if(length(bf[[region]]$pos)){
                    GRanges(ranges=IRanges(start=bf[[region]]$pos, width=cigarToWidth(bf[[region]]$cigar)),
                            strand=bf[[region]]$strand,
                            seqnames=bf[[region]]$rname,
                            names=bf[[region]]$qname,
                            IH=if(is.null(bf[[region]]$tag$IH))
                                      1
                                else
                                      bf[[region]]$tag$IH      
                            )
                } else
                    GRanges()
                }))
            reads <- unlist(reads, use.names=F)
        }))
        names(reads) <- names(scanBamRes[[1]])
           
        cnt <- unlist(lapply(names(range), function(region){
            if(shift > 0 && length(reads[[region]])){
                reads[[region]] <- shift(reads[[region]],
                                         ifelse(as.vector(strand(reads[[region]]) == "+"), shift, -shift))
            }

            idx <- switch(overlap,
                   startwithin = {
                       start(reads[[region]]) <- ifelse(as.vector(strand(reads[[region]]) == "+"), 
                                                   start(reads[[region]]), end(reads[[region]]))
                       width(reads[[region]]) <- 1
                       .isWithin(range[region], reads[[region]], 1L)
                   }, 
                   endwithin = {
                       start(reads[[region]]) <- ifelse(as.vector(strand(reads[[region]]) == "-"), 
                                                   start(reads[[region]]), end(reads[[region]]))
                       width(reads[[region]]) <- 1
                       .isWithin(range[region], reads[[region]], 1L)
                   },
                   within = .isWithin(range[region], reads[[region]], minoverlap),
                   any = .isAny(range[region], reads[[region]], minoverlap)
                   )
            if(!is.null(maxHits))
                idx <- idx & values(reads[[region]])$IH <= maxHits
            if(stranded && as.vector(strand(range[region]) != "*"))
                idx <- idx & as.vector(strand(range[region]) == strand(reads[[region]]))                          
            sum(1/values(reads[[region]][idx,])$IH)
        }))
        return(cnt)
                    
    }))
    names(count) <- names(grange)

    if("collapseQuery" %in% names(elementMetadata(grange))){
        count <- unlist(lapply(split(count, elementMetadata(grange)[,"collapseQuery"]), sum))
    }
    return(count)
}
            
.isWithin <- function(grange, reads, minoverlap){
    start(grange) <= start(reads) &
        end(reads) <= end(grange) &
        minoverlap <= width(reads)
}

.isAny <- function(grange, reads, minoverlap){
    mo <- minoverlap - 1L
    end(reads) - start(grange) >= mo &
        end(grange) - start(reads) >= mo &
        minoverlap <= width(reads)
}


.countAlignmentsC <- function(bamfile, grange, stranded=FALSE, overlap=c("any", "within", "startwithin", "endwithin"), shift=0L, minoverlap=1L, maxHits=NULL)
{
    ## check arguments
    if(!is.character(bamfile))
        stop("The parameter 'bamfile' must be of type character.")
    if(!is(grange, "GRanges"))
        stop("The parameter 'grange' must be of type GRanges.")
    if(!is.logical(stranded))
        stop("The parameter 'stranded' must be of type logical.")
    overlap <- match.arg(overlap)
    if(length(shift) != 1L && length(shift) != length(bamfile))
        stop("The parameter 'shift' must be of length one or of equal length as the parameter 'file'.")
    if(length(shift) == 1L)
        shift <- rep(shift, length(bamfile))
    shift <- as.integer(shift)
    minoverlap <- as.integer(minoverlap)
    if(is.null(maxHits))
        maxHits <- 0L
    maxHits <- as.integer(maxHits)
    # TODO check parameters

    ## translate seqname to tid
    seqnamemap <- .Call(.seqname, bamfile[1])
    regions <- data.frame(tid=seqnamemap$tid[ IRanges::as.vector(IRanges::match(seqnames(grange), seqnamemap$seqnames)) ],
                          start=as.integer(start(grange)),
                          end=as.integer(end(grange)),
                          strand=as.character(strand(grange)),
                          stringsAsFactors=FALSE
                          )
    ## get counts
    ## TODO shift(grange, shift)
    count <- .Call(.count_alignments, bamfile, bamfile, regions, stranded, overlap, minoverlap, shift, maxHits)
    if("collapseQuery" %in% names(elementMetadata(grange))){
        count <- unlist(lapply(split(count, elementMetadata(grange)[,"collapseQuery"]), sum))
    }
    #rownames(counts) <- names(query)

    return(count)
}
