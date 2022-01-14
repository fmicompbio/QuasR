## query : GRanges : value is list of matrices:
##                   one list element per sample with matrix:
##                            one(normal) or three(allelic) rows per unique query name and max(width(query)) columns
##                   e.g.   res[['Sample1']]['HCP_TSS','-100']

#' Quantify alignments by relative position
#' 
#' Quantify alignments from sequencing data, relative to their position in 
#' query regions.
#' 
#' \code{qProfile} is used to count alignments in each sample from a 
#' \code{qProject} object, relative to their position in query regions.
#' 
#' Most arguments are identical to the ones of \code{\link[QuasR]{qCount}}.
#' 
#' The \code{query} argument is a \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' object that defines the regions for the profile. All regions in
#' \code{query} will be aligned to one another at their anchor position,
#' which corresponds to their biological start position (\code{start(query)}
#' for regions on strand \dQuote{+} or \dQuote{*}, \code{end(query)} for 
#' regions on strand \dQuote{-}).
#' 
#' This anchor position will be extended (with regard to strand) by
#' the number of bases specified by \code{upstream} and \code{downstream}.
#' In the return value, the anchor position will be at position zero.
#' 
#' If \code{binSize} is greater than one, \code{upstream} and \code{downstream}
#' will be slightly increased in order to include the complete first and last
#' bins of \code{binSize} bases.
#' 
#' Regions with identical names in \code{names{query}} will be summed, and
#' profiles will be padded with zeros to accomodate the length of all profiles.
#' 
#' @param proj A \code{\linkS4class{qProject}} object representing a 
#'   sequencing experiment as returned by \code{\link[QuasR]{qAlign}}
#' @param query An object of type \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'   with the regions to be profiled. All regions in \code{query} will be
#'   anchored at their biological start position (\code{start(query)} for
#'   regions on strand \dQuote{+} or \dQuote{*}, \code{end(query)} for
#'   regions on strand \dQuote{-}). This position will become position zero 
#'   in the return value.
#' @param upstream An \dQuote{integer} vector of length one or the same 
#'   length as \code{query} indicating the number of bases upstream of the 
#'   anchor position to include in the profile.
#' @param downstream An \dQuote{integer} vector of length one or the same 
#'   length as \code{query} indicating the number of bases downstream of the
#'   anchor position to include in the profile.
#' @param selectReadPosition defines the part of the alignment that has to be 
#'   contained within a query region to produce an overlap (see Details), and 
#'   that is used to calculate the relative position within the query region. 
#'   Possible values are:
#'   \itemize{
#'     \item \code{start} (default): start of the alignment
#'     \item \code{end}: end of the alignment
#'   }
#' @param shift controls the shifting alignments towards their 3'-end before 
#'   quantification. \code{shift} can be one of:
#'   \itemize{
#'     \item an \dQuote{integer} vector of the same length as the
#'     number of alignment files
#'     \item a single \dQuote{integer} value
#'     \item the character string \code{"halfInsert"} (only available for
#'     paired-end experiments)
#'   }
#'   The default of \code{0} will not shift any alignments.
#' @param orientation sets the required orientation of the alignments relative 
#'   to the query region in order to be counted, one of:
#'   \itemize{
#'     \item \code{any} (default): count alignment on the same and opposite strand
#'     \item \code{same} : count only alignment on the same strand
#'     \item \code{opposite} : count only alignment on the opposite strand
#'   }
#' @param useRead For paired-end experiments, selects the read mate whose 
#'   alignments should be counted, one of:
#'   \itemize{
#'     \item \code{any} (default): count all alignments
#'     \item \code{first} : count only alignments from the first read
#'     \item \code{last} : count only alignments from the last read
#'   }
#' @param auxiliaryName Which bam files to use in an experiments with 
#'   auxiliary alignments (see Details).
#' @param mask If not \code{NULL}, a \code{\link[GenomicRanges:GRanges-class]{GRanges}} 
#'   object with reference regions to be masked, i.e. excluded from the 
#'   quantification, such as unmappable or highly repetitive regions (see 
#'   Details).
#' @param collapseBySample If \code{TRUE} (the default), sum alignment 
#'   counts from bam files with the same sample name.
#' @param includeSpliced If \code{TRUE} (the default), include spliced 
#'   alignments when counting. A spliced alignment is defined as an 
#'   alignment with a gap in the read of at least 60 bases.
#' @param includeSecondary If \code{TRUE} (the default), include alignments 
#'   with the secondary bit (0x0100) set in the \code{FLAG} when counting.
#' @param mapqMin Minimal mapping quality of alignments to be included when 
#'   counting (mapping quality must be greater than or equal to \code{mapqMin}). 
#'   Valid values are between 0 and 255. The default (0) will include all 
#'   alignments.
#' @param mapqMax Maximal mapping quality of alignments to be included when 
#'   counting (mapping quality must be less than or equal to \code{mapqMax}).
#'   Valid values are between 0 and 255. The default (255) will include all 
#'   alignments.
#' @param absIsizeMin For paired-end experiments, minimal absolute insert 
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included when 
#'   counting. Valid values are greater than 0 or \code{NULL} (default), 
#'   which will not apply any minimum insert size filtering.
#' @param absIsizeMax For paired-end experiments, maximal absolute insert 
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included when 
#'   counting. Valid values are greater than 0 or \code{NULL} (default), 
#'   which will not apply any maximum insert size filtering.
#' @param maxInsertSize Maximal fragment size of the paired-end experiment.  
#'   This parameter is used if \code{shift="halfInsert"} and will ensure that 
#'   query regions are made wide enough to emcompass all alignment pairs whose 
#'   mid falls into the query region. The default value is \code{500} bases.
#' @param binSize Numeric scalar giving the size of bins (must be an odd number).
#'   The default value (\code{1}) gives back counts for single bases. Otherwise,
#'   alignments are counted in adjacent, non-overlapping windows of size 
#'   \code{binSize} that tile the interval defined by \code{upstream} and 
#'   \code{downstream}.
#' @param clObj A cluster object to be used for parallel processing (see 
#'   \sQuote{Details}).
#' 
#' @name qProfile
#' @aliases qProfile
#' 
#' @return 
#' A \code{list} of matrices with \code{length(unique(names(query)))} rows
#' with profile names, and \code{max(upstream)+max(downstream)+1} columns
#' indicating relative position (for \code{binsize=1}).
#' 
#' For \code{binSize} values greater than 1, the number of columns corresponds to
#' the number of bins (tiles), namely 
#' \code{ceiling(max(upstream)/binSize)+ceiling(max(downstream)/binSize)}.
#' A middle bin of size \code{binSize} is always positioned centered at the anchor
#' of each region. Additional bins are positioned upstream and downstream, adjacent
#' to that middle bin, in order to include at least \code{upstream} and
#' \code{downstream} bases, respectively (potentially more in order to fill the 
#' first and last bins).
#' 
#' The relative positions are given as column names (for \code{binSize > 1}
#' they refer to the bin mid). In that case, the bins are "right-open". For
#' example, if \code{binSize = 10}, the bin with the midpoint "-50" contains
#' counts for the alignments in [-55,-45).
#' 
#' The first list element is called \dQuote{coverage} and contains, for each
#' profile and relative position, the number of overlapping regions
#' that contributed to the profile.
#' 
#' Subsequent list elements contain the alignment counts for individual
#' sequence files (\code{collapseBySample=FALSE}) or samples 
#' (\code{collapseBySample=TRUE}) in \code{proj}.
#' 
#' For projects with allele-specific quantification, i.e. if a file with
#' single nucleotide polymorphisms was supplied to the \code{snpFile}
#' argument of \code{\link[QuasR]{qAlign}}, there will be three rows
#' instead of one row with counts per unique region name, with numbers
#' of alignments for Reference, Unknown and Alternative genotypes
#' (suffixed _R, _U and _A).
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' 
#' @keywords utilities misc
#' 
#' @export
#' 
#' @seealso 
#' \code{\link[QuasR]{qCount}},
#' \code{\link[QuasR]{qAlign}},
#' \code{\linkS4class{qProject}},
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
#' 
#' @examples 
#' # copy example data to current working directory
#' file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
#' 
#' # create alignments (single-end experiment)
#' genomeFile <- "extdata/hg19sub.fa"
#' sampleFile <- "extdata/samples_chip_single.txt"
#' proj <- qAlign(sampleFile, genomeFile)
#' 
#' # load transcript start site coordinates
#' library(rtracklayer)
#' annotationFile <- "extdata/hg19sub_annotation.gtf"
#' tssRegions <- import.gff(annotationFile, format="gtf",
#'                          feature.type="start_codon")
#' 
#' # obtain a combined TSS profile
#' pr1 <- qProfile(proj, tssRegions)
#' lapply(pr1, dim)
#' lapply(pr1, "[", , 1:5)
#' 
#' prComb <- do.call("+", lapply(pr1[-1], function(x) x/pr1[[1]]))
#' barplot(prComb, xlab="Position", ylab="Mean no. of alignments")
#' 
#' # obtain TSS profiles for individual regions
#' names(tssRegions) <- mcols(tssRegions)$transcript_id
#' pr2 <- qProfile(proj, tssRegions)
#' lapply(pr2, dim)
#' lapply(pr2, "[", 1:3, 1:5)
#' 
#' @importFrom Rsamtools scanBamHeader
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end strand
#' @importFrom GenomeInfoDb seqnames seqlevels
#' @importFrom parallel clusterEvalQ clusterMap
#' 
qProfile <- function(proj,
                     query,
                     upstream = 1000,
                     downstream = upstream,
                     selectReadPosition = c("start", "end"),
                     shift = 0L,
                     orientation = c("any", "same", "opposite"),
                     useRead = c("any", "first", "last"),
                     auxiliaryName = NULL,
                     mask = NULL,
                     collapseBySample = TRUE,
                     includeSpliced = TRUE,
                     includeSecondary = TRUE,
                     mapqMin = 0L,
                     mapqMax = 255L,
                     absIsizeMin = NULL,
                     absIsizeMax = NULL,
                     maxInsertSize = 500L,
                     binSize = 1L,
                     clObj = NULL) {
    ## setup variables from 'proj' ---------------------------------------------
    ## 'proj' is correct type?
    if (!inherits(proj, "qProject", which = FALSE))
        stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")
    
    samples <- proj@alignments$SampleName
    nsamples <- length(samples)
    bamfiles <-
        if (is.null(auxiliaryName))
            proj@alignments$FileName
    else if (!is.na(i <- match(auxiliaryName, proj@aux$AuxName)))
        unlist(proj@auxAlignments[i, ], use.names = FALSE)
    else
        stop("unknown 'auxiliaryName', should be one of: NULL, ",
             paste(sprintf("'%s'", proj@aux$AuxName), collapse = ", "))

        
    ## validate parameters -----------------------------------------------------
    if (!is.numeric(upstream))
        stop("'upstream' must be of type 'numeric'")
    if (!is.numeric(downstream))
        stop("'downstream' must be of type 'numeric'")
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    selectReadPosition <- match.arg(selectReadPosition)
    orientation <- match.arg(orientation)
    useRead <- match.arg(useRead)
    if (!is.logical(collapseBySample) || length(collapseBySample) != 1L)
        stop("'collapseBySample' must be either TRUE or FALSE")
    if (!is.logical(includeSpliced) || length(includeSpliced) != 1L)
        stop("'includeSpliced' must be either TRUE or FALSE")
    if (!is.logical(includeSecondary) || length(includeSecondary) != 1L)
        stop("'includeSecondary' must be either TRUE or FALSE")
    if ((!is.null(absIsizeMin) || !is.null(absIsizeMax)) && proj@paired == "no")
        stop("'absIsizeMin' and 'absIsizeMax' can only be used for paired-end experiments")
    if (is.null(absIsizeMin)) # -1L -> do not apply TLEN filtering
        absIsizeMin <- -1L
    if (is.null(absIsizeMax))
        absIsizeMax <- -1L
    if (!is.numeric(maxInsertSize) || length(maxInsertSize) != 1L)
        stop("'maxInsertSize' must be a numerical scalar")
    if (!is.numeric(binSize) || length(binSize) != 1L || binSize <= 0)
        stop("'binSize' must be a single numerical value greater than zero")
    if (binSize %% 2 != 1)
        stop("'binSize' must be an odd number")
    
    ## check shift
    if (length(shift) == 1 && shift == "halfInsert") {
        if (proj@paired == "no") {
            stop("'shift=\"halfInsert\"' can only be used for paired-end experiments")
        } else {
            shifts <- rep(-1000000L, nsamples)
            broaden <- as.integer(ceiling(maxInsertSize/2))
        }
    } else {
        if (!is.numeric(shift) || (length(shift) > 1 && length(shift) != nsamples))
            stop(sprintf("'shift' must be 'halfInsert', a single integer or an integer vector with %d values", nsamples))
        else if (length(shift) == 1)
            shifts <- rep(as.integer(shift), nsamples)
        else
            shifts <- as.integer(shift)
        broaden <- 0L
    }
    
    ## all query chromosomes present in all bamfiles?
    trTab <- table(unlist(lapply(Rsamtools::scanBamHeader(bamfiles), 
                                 function(bh) names(bh$targets))))
    trCommon <- names(trTab)[trTab == length(bamfiles)]
    if (any(f <- !(GenomeInfoDb::seqlevels(query) %in% trCommon)))
        stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                     paste(GenomeInfoDb::seqlevels(query)[f], collapse = ", ")))
    
    ## 'useRead' set but not a paired-end experiment?
    if (useRead != "any" && proj@paired == "no")
        warning("ignoring 'useRead' for single read experiments")
    
    
    ## preprocess query --------------------------------------------------------
    ##    --> extract 'querynames', 'refpos', 'queryWin', 'maxUp', 'maxDown', 'maxUpBin', 'maxDownBin'
    ##        split into 'queryWinL', 'refposL', 'querynamesIntL' by name
    ##    GRanges query --------------------------------------------------------
    if (inherits(query,"GRanges")) {
        # define 'querynames'
        if (!is.null(names(query)))
            querynames <- names(query)
        else
            querynames <- rep("query", length(query)) # combine all regions
        
        # define 'refpos' and 'queryWin'
        if (length(upstream) == 1)
            upstream <- rep(upstream, length(query))
        else if (length(upstream) != length(query))
            stop(sprintf("the length of 'upstream' (%d) must be one or the same as 'query' (%d)",
                         length(upstream), length(query)))
        if (length(downstream) == 1)
            downstream <- rep(downstream, length(query))
        else if (length(downstream) != length(query))
            stop(sprintf("the length of 'downstream' (%d) must be one or the same as 'query' (%d)",
                         length(downstream), length(query)))
        
        # adjust 'upstream', 'downstream' to include full bins
        upstream <- ceiling((upstream - (binSize - 1)/2) / binSize) * binSize + 
            (binSize - 1)/2
        downstream <- ceiling((downstream - (binSize - 1)/2) / binSize) * binSize + 
            (binSize - 1)/2
        
        # define 'maxUp', 'maxDown', 'maxUpBin' and 'maxDownBin'
        # (have one bin that is centered at anchor point, relative position of bin mid at zero)
        #   maxUp, maxDown: number of bases upstream/downstream of start(query) --> includes half of middle bin, but not the anchor posittion
        #                   (total number of covered bases: maxUp + maxDown + 1)
        #   maxUpBin, maxDownBin: number of bins upstream/downstream of middle bin
        #                         (total number of bins: maxUpBin + maxDownBin + 1)
        maxUpBin   = as.integer(ceiling((max(upstream)   + 1 - (binSize + 1)/2)/binSize))
        maxDownBin = as.integer(ceiling((max(downstream) + 1 - (binSize + 1)/2)/binSize))
        maxUp   = as.integer(maxUpBin   * binSize + (binSize - 1)/2)
        maxDown = as.integer(maxDownBin * binSize + (binSize - 1)/2)
        # err = (maxUp + maxDown + 1) %% binSize # has to be zero
        binNames <- as.character(seq(-maxUpBin * binSize, maxDownBin * binSize, by = binSize))
        if ((maxUpBin + maxDownBin + 1) > 20000L)
            warning(sprintf("profiling over large region (%d basepairs, %d bins) - may be slow and require a lot of memory",
                            maxUp + maxDown + 1, maxUpBin + maxDownBin))
        
        plusStrand <- as.character(BiocGenerics::strand(query)) != "-"
        refpos <- ifelse(plusStrand, BiocGenerics::start(query), 
                         BiocGenerics::end(query))
        queryWin <- GenomicRanges::GRanges(
            GenomeInfoDb::seqnames(query),
            IRanges::IRanges(
              start = pmax(1, ifelse(plusStrand, 
                                     refpos - upstream, refpos - downstream)),
              end = ifelse(plusStrand, refpos + downstream, refpos + upstream)),
            strand = BiocGenerics::strand(query)
        )

        # split 'queryWin' --> 'queryWinL' and 'refpos' --> 'refposL' by 'querynames'
        if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
            profileChunkL <- split(seq_len(length(queryWin)), 
                                   factor(querynames, levels = unique(querynames)))
            approxNumRegionsPerChunk <- length(queryWin) / (length(clObj) / length(bamfiles))
            profileChunkToTask <- round(cumsum(vapply(profileChunkL, length, 1))
                                        / approxNumRegionsPerChunk)
            taskL <- lapply(split(seq_along(profileChunkL), profileChunkToTask), function(i) do.call(c, profileChunkL[i]))
        } else {
            taskL <- list(seq_len(length(queryWin))) # single task per bam file
        }
        queryWinL <- lapply(taskL, function(i) queryWin[i])
        refposL <- lapply(taskL, function(i) refpos[i])
        querynamesInt <- as.integer(factor(querynames, levels = unique(querynames)))
        querynamesIntL <- lapply(taskL, function(i) querynamesInt[i])
        
        # calculate coverage
        cvgBaseL <- lapply(seq_along(queryWin), function(i) 
            (maxUp - upstream[i] + 1):(maxUp + downstream[i] + 1))
        cvgBase <- do.call(
          rbind, lapply(split(seq_along(queryWin),
                              factor(querynames, levels = unique(querynames))),
                        function(i) tabulate(unlist(cvgBaseL[i]), 
                                             nbins = maxUp + maxDown + 1))
        )
        cvg <- do.call(
          cbind, lapply(split(seq.int(maxUp + maxDown + 1), 
                              rep(seq(-maxUpBin, maxDownBin), each = binSize)),
                        function(i) rowSums(cvgBase[, i, drop = FALSE]))
        )
        colnames(cvg) <- binNames
        
    } else {
        stop("'query' must be an object of type 'GRanges'")
    }
    ## from now on, only use 'queryWinL' (named list of GRanges objects) with reference positions in 'refposL'
    
    
    ## apply 'mask' to queryWinL -----------------------------------------------
    if (!is.null(mask)) {
        if (!inherits(mask, "GRanges"))
            stop("'mask' must be an object of type 'GRanges'")
        stop("'mask' is not yet supported for qProfile")
    }        
    
    
    ## setup tasks for parallelization -----------------------------------------
    if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
        loadQuasR(clObj)
        myapply <- function(...) parallel::clusterMap(clObj, ..., SIMPLIFY = FALSE, 
                                                      .scheduling = "dynamic")
    } else {
        myapply <- function(...) mapply(..., SIMPLIFY = FALSE)
    }
    
    ## count alignments --------------------------------------------------------
    message("profiling alignments...", appendLF = FALSE)
    res <- myapply(profileAlignments,
                   bamfile = rep(bamfiles, each = length(queryWinL)),
                   queryids = querynamesIntL,
                   regions = queryWinL,
                   refpos = refposL,
                   shift = rep(shifts, each = length(queryWinL)),
                   MoreArgs = list(
                       selectReadPosition = selectReadPosition,
                       orientation = orientation,
                       useRead = useRead,
                       broaden = broaden,
                       allelic = !is.na(proj@snpFile),
                       maxUp = maxUp,
                       maxDown = maxDown,
                       maxUpBin = maxUpBin,
                       maxDownBin = maxDownBin,
                       includeSpliced = includeSpliced,
                       includeSecondary = includeSecondary,
                       mapqmin = as.integer(mapqMin)[1],
                       mapqmax = as.integer(mapqMax)[1],
                       absisizemin = as.integer(absIsizeMin)[1],
                       absisizemax = as.integer(absIsizeMax)[1],
                       binsize = as.integer(binSize),
                       binNames = binNames))
    message("done")
    
    ## fuse by input file, rename and collapse by sample
    if (!is.na(proj@snpFile)) {
        res <- lapply(split(seq_along(res), rep(displayNames(proj), each = length(queryWinL))),
                      function(i) {
                          tmpR <- do.call(rbind, lapply(res[i], "[[", "R"))
                          tmpU <- do.call(rbind, lapply(res[i], "[[", "U"))
                          tmpA <- do.call(rbind, lapply(res[i], "[[", "A"))
                          rownames(tmpR) <- rownames(tmpU) <- rownames(tmpA) <-
                              unique(querynames)
                          list(R = tmpR, U = tmpU, A = tmpA)
                      })
        if (collapseBySample) {
            res <- lapply(split(seq_along(res), factor(samples, levels = unique(samples))),
                          function(i) list(R = Reduce("+", lapply(res[i], "[[", "R")),
                                           U = Reduce("+", lapply(res[i], "[[", "U")),
                                           A = Reduce("+", lapply(res[i], "[[", "A"))))
        }
        nms <- paste(rep(names(res), each = 3), c("R", "U", "A"), sep = "_")
        res <- do.call(c, res)
        names(res) <- nms
    } else {
        res <- lapply(split(seq_along(res), rep(displayNames(proj), each = length(queryWinL))),
                      function(i) {
                          tmp <- do.call(rbind, res[i])
                          rownames(tmp) <- unique(querynames)
                          tmp
                      })
        if (collapseBySample) {
            res <- lapply(split(seq_along(res), factor(samples, levels = unique(samples))),
                          function(i) Reduce("+", res[i]))
        }
    }
    
    ## add the region coverage as first elemnt
    res <- c(list(coverage = cvg), res)
    
    ## return results
    return(res)
}

## profile alignments (with the C-function) for single bamfile, multiple regions, single shift
## return a numeric vector with maxWidth elements corresponding to the positions within the query regions
#' @keywords internal
#' @importFrom Rsamtools scanBamHeader
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand match as.vector
profileAlignments <- function(bamfile, queryids, regions, refpos, shift, 
                              selectReadPosition, orientation, useRead, 
                              broaden, allelic, maxUp, maxDown, maxUpBin, 
                              maxDownBin, includeSpliced, includeSecondary,
                              mapqmin, mapqmax, absisizemin, absisizemax, 
                              binsize, binNames) {
    tryCatch({ # try catch block contains whole function
        
        # translate seqnames to tid and create region data.frame
        seqnamesBamHeader <- names(Rsamtools::scanBamHeader(bamfile)[[1]]$targets)
        
        # prepare region vectors
        tid <- BiocGenerics::as.vector(BiocGenerics::match(GenomeInfoDb::seqnames(regions),
                                                           seqnamesBamHeader) - 1L) 
        s <- BiocGenerics::start(regions) - 1L # Samtools library has 0-based inclusive start
        e <- BiocGenerics::end(regions) # Samtools library has 0-based exclusive end
        rp <- refpos - 1L # Samtools library has 0-based inclusive start
        
        ## swap selstrand for 'orientation="opposite"'
        regstrand <- as.character(BiocGenerics::strand(regions))
        if (orientation == "any")
            selstrand <- rep("*", length(regions))
        else if (orientation == "opposite")
            selstrand <- c("+"="-", "-"="+", "*"="*")[regstrand]
        else # orientation == "same"
            selstrand <- regstrand
        
        ## translate useRead parameter
        BAM_FREAD1 <- 64L
        BAM_FREAD2 <- 128L
        if (useRead == "any")
            readBitMask <- BAM_FREAD1 + BAM_FREAD2
        else if (useRead == "first")
            readBitMask <- BAM_FREAD1
        else if (useRead == "last")
            readBitMask <- BAM_FREAD2

        ## translate includeSecondary parameter
        BAM_FSECONDARY <- 256L
        if (includeSecondary)
            readBitMask <- readBitMask + BAM_FSECONDARY
  
        ## count alignments by position
        if (!allelic) {
            count <- t(.Call(profileAlignmentsNonAllelic, bamfile, queryids, tid, 
                             s, e, rp, selstrand, regstrand,
                             selectReadPosition, readBitMask, shift, broaden, 
                             maxUp, maxDown, maxUpBin, maxDownBin,
                             includeSpliced, mapqmin, mapqmax, absisizemin, 
                             absisizemax, binsize, binNames))
        } else {
            count <- lapply(.Call(profileAlignmentsAllelic, bamfile, queryids, tid, 
                                  s, e, rp, selstrand, regstrand,
                                  selectReadPosition, readBitMask, shift, broaden, 
                                  maxUp, maxDown, maxUpBin, maxDownBin,
                                  includeSpliced, mapqmin, mapqmax, absisizemin, 
                                  absisizemax, binsize, binNames), t)
        }
        
        return(count)

    }, error = function(ex) {
        reg <- regions[c(1, length(regions))]
        emsg <- paste("Internal error on ", Sys.info()['nodename'], 
                      ", bamfile ", bamfile," with regions\n\t", 
                      paste(GenomeInfoDb::seqnames(reg), ":", 
                            BiocGenerics::start(reg), "-" , BiocGenerics::end(reg),
                            ":", BiocGenerics::strand(reg), sep = "", 
                            collapse = "\n\t...\n\t"), 
                      "\n Error message is: ", ex$message, sep = "")
        stop(emsg)
    })
}
