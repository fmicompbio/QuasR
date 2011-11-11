.countAlignments <- function(bamfile, grange, stranded=FALSE,
                             overlap=c("any", "within", "startwithin", "endwithin"),
                             shift=0L, minoverlap=1L, maxHits=NULL){
    if(getOption("quasr.countApp") == "R")
        return(.countAlignmentsR(bamfile=bamfile, grange=grange, stranded=stranded,
                                 overlap=overlap,  shift=shift, minoverlap=minoverlap, maxHits=maxHits))
    else
        return(.countAlignmentsC(bamfile=bamfile, grange=grange, stranded=stranded,
                                 overlap=overlap,  shift=shift, minoverlap=minoverlap, maxHits=maxHits))
}

.countAlignmentsR <- function(bamfile, grange, stranded=FALSE,
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

    count <- unlist(lapply(split(grange,1), function(range){
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

        count <- unlist(lapply(seq(scanBamRes[[1]]), function(iRegion){
            cnt <- unlist(lapply(seq(scanBamRes), function(iBamfile){
                ## name of current region
                regionName <- names(scanBamRes[[iBamfile]])[iRegion]
                ## calculate width
                readWidth <- cigarToWidth(scanBamRes[[iBamfile]][[iRegion]]$cigar)
                ## shift
                readIsPlusStrand <- scanBamRes[[iBamfile]][[iRegion]]$strand == "+"
                if(shift > 0){ # TODO check if length(shift) > 1
                    scanBamRes[[iBamfile]][[iRegion]]$pos <-
                        scanBamRes[[iBamfile]][[iRegion]]$pos +
                            ifelse(readIsPlusStrand, shift, -shift)
                }
                ## in region
                idx <- switch(overlap,
                              startwithin = {
                                  scanBamRes[[iBamfile]][[iRegion]]$pos <-
                                      ifelse(readIsPlusStrand,
                                             scanBamRes[[iBamfile]][[iRegion]]$pos,
                                             scanBamRes[[iBamfile]][[iRegion]]$pos + readWidth)
                                  start(range[regionName]) <= scanBamRes[[iBamfile]][[iRegion]]$pos &
                                  scanBamRes[[iBamfile]][[iRegion]]$pos <= end(range[regionName])
                              },
                              endwithin = {
                                  scanBamRes[[iBamfile]][[iRegion]]$pos <-
                                      ifelse(readIsPlusStrand,
                                             scanBamRes[[iBamfile]][[iRegion]]$pos + readWidth,
                                             scanBamRes[[iBamfile]][[iRegion]]$pos)
                                  start(range[regionName]) <= scanBamRes[[iBamfile]][[iRegion]]$pos &
                                  scanBamRes[[iBamfile]][[iRegion]]$pos <= end(range[regionName])
                              },
                              within = {
                                  start(range[regionName]) <=
                                      scanBamRes[[iBamfile]][[iRegion]]$pos &
                                  scanBamRes[[iBamfile]][[iRegion]]$pos + readWidth <= end(range[regionName]) &
                                  minoverlap <= readWidth
                              },
                              any = {
                                  mo <- minoverlap - 1L
                                  (scanBamRes[[iBamfile]][[iRegion]]$pos + readWidth) - start(range[regionName]) >= mo &
                                  end(range[regionName]) - scanBamRes[[iBamfile]][[iRegion]]$pos >= mo &
                                  minoverlap <= readWidth
                              }
                              )
                ## check maxhits
                if(!is.null(maxHits) && !is.null(scanBamRes[[iBamfile]][[iRegion]]$tag$IH))
                    idx <- idx & scanBamRes[[iBamfile]][[iRegion]]$tag$IH <= maxHits
                ## check stranded
                if(stranded && as.vector(strand(range[regionName]) != "*"))
                    idx <- idx & as.vector(strand(range[regionName]) == scanBamRes[[iBamfile]][[iRegion]]$strand)
                ## sum weigths
                if(is.null(scanBamRes[[iBamfile]][[iRegion]]$tag$IH)) ## set IH tag to 1 if missing
                    sum(idx)
                else
                    sum(1/scanBamRes[[iBamfile]][[iRegion]]$tag$IH[idx])
            }))
            sum(cnt)
        }))
    }))
    names(count) <- names(grange)

    if("collapseQuery" %in% names(elementMetadata(grange))){
        count <- unlist(lapply(split(count, elementMetadata(grange)[,"collapseQuery"]), sum))
    }
    return(count)
}

.countAlignmentsC <- function(bamfile, grange, stranded=FALSE,
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
                          start=as.integer(start(grange)-1), ## samtool sw has 0-based start
                          end=as.integer(end(grange)),
                          strand=as.character(strand(grange)),
                          stringsAsFactors=FALSE
                          )
    ## get counts
    count <- .Call(.count_alignments, bamfile, bamfile, regions, stranded, overlap, minoverlap, shift, maxHits)
#     system.time(count <- .Call(.count_alignments, bamfile, bamfile, regions, stranded, overlap, minoverlap, shift, maxHits))
#     browser()
#     system.time(countt <- unlist(mclapply(1:dim(regions)[1], function(idx){
#         .Call(.count_alignments, bamfile, bamfile, regions[idx,], stranded, overlap, minoverlap, shift, maxHits)
#         }, mc.cores=getOption("quasr.cores"))))
    names(count) <- names(grange)

    ## collapse counts
    if("collapseQuery" %in% names(elementMetadata(grange))){
        count <- unlist(lapply(split(count, elementMetadata(grange)[,"collapseQuery"]), sum))
    }
    return(count)
}
