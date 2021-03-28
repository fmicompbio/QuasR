# export alignment to wig file
# - each read (pair) "lives" on a single base:
#     pairs on the middle of the fragment
#     singles on the 5'-end base shifted by "shift" towards their 3'-end
# - multiple samples can be automatically normalized to one another
# - multiple bam files with identical sample name can be combined into the same file (track)
# - wig files can be compressed
#
# proj        : qProject object
# file        : wig file name(s) (will be compressed if file QuasR:::compressedFileFormat(file) != "none"
# collapseBySample : create one track per unique sample name
# binsize     : stepInterval and windowSize for the fixedStep wig file
# shift       : only for single read projects; shift read
# strand      : only include alignments on '+', '-' or '*' (any) strand
# scaling     : scale multiple tracks to one another?
# tracknames  : names for display in track header
# log2p1      : transform alignment count by log2(x+1)
# colors      : colors for tracks
# mapqMin     : minimum mapping quality (MAPQ >= mapqMin)
# mapqMax     : maximum mapping quality (MAPQ <= mapqMax)
# absIsizeMin : minimum absolute insert size (TLEN >= absIsizeMin)
# absIsizeMax : maximum absolute insert size (TLEN <= absIsizeMax)
# useRead     : for paired-end data, what read to use
# pairedAsSingle : for paired-end data, treat as single read data (do not calculate fragment mid-points)
#' QuasR wig file export
#' 
#' Create a fixed-step wig file from the alignments in the genomic bam files 
#' of the \sQuote{QuasR} project.
#' 
#' \code{qExportWig()} uses the genome bam files in \code{proj} as input
#' to create wig or bigWig files with the number of alignments (pairs)
#' per window of \code{binsize} nucleotides. By default
#' (\code{collapseBySample=TRUE}), one file per unique sample will be
#' created. If \code{collapseBySample=FALSE}, one file per genomic bam
#' file will be created. See \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html} 
#' for the definition of the wig format, and 
#' \url{http://genome.ucsc.edu/goldenPath/help/bigWig.html} for the definition 
#' of the bigWig format.
#' 
#' The genome is tiled with sequential windows of length \code{binsize},
#' and alignments in the bam file are assigned to these windows: Single
#' read alignments are assigned according to their 5'-end coordinate
#' shifted by \code{shift} towards the 3'-end (assuming that the 5'-end
#' is the leftmost coordinate for plus-strand alignments, and the rightmost
#' coordinate for minus-strand alignments). Paired-end alignments are
#' assigned according to the base in the middle between the leftmost and
#' rightmost coordinates of the aligned pair of reads. Each pair of reads
#' is only counted once, and not properly paired alignments are
#' ignored. If \code{useRead} is set to select only the first or last
#' read in a paired-end experiment, the selected read will be treated as
#' reads from a single read experiment. Secondary alignments can be
#' excluded by setting \code{includeSecondary=FALSE}. In paired-end
#' experiments, \code{absIsizeMin} and \code{absIsizeMax} can be used to select
#' alignments based on their insert size (TLEN field in SAM Spec v1.4).
#' 
#' For \code{scaling=TRUE}, the number of alignments per bin \eqn{n}
#' for the sample \eqn{i} are linearly scaled to the mean total
#' number of alignments over all samples in \code{proj} according to:
#' \eqn{n_s = n /N[i] *mean(N)} where \eqn{n_s} is the scaled number
#' of alignments in the bin and \eqn{N} is a vector with the total
#' number of alignments for each sample. Alternatively, if scaling is set
#' to a positive numerical value \eqn{s}, this value is used instead of
#' \eqn{\textnormal{mean}(N)}{mean(N)}, and values are scaled according
#' to: \eqn{n_s = n /N[i] *s}.
#' 
#' \code{mapqMin} and \code{mapqMax} allow to select alignments
#' based on their mapping qualities. \code{mapqMin} and \code{mapqMax} can
#' take integer values between 0 and 255 and equal to
#' \eqn{-10 log_{10} Pr(\textnormal{mapping position is wrong})}{-10
#' log10 Pr(mapping position is wrong)}, rounded to the nearest
#' integer. A value 255 indicates that the mapping quality is not available.
#' 
#' If \code{createBigWig=FALSE} and \code{file} ends with \sQuote{.gz},
#' the resulting wig file will be compressed using gzip and is suitable
#' for uploading as a custom track to your favorite genome browser
#' (e.g. UCSC or Ensembl).
#' 
#' @param proj A \code{qProject} object as returned by \code{qAlign}.
#' @param file A character vector with the name(s) for the wig or bigWig
#'   file(s) to be generated. Either \code{NULL} or a vector of the same
#'   length as the number of bam files (for \code{collapseBySample=FALSE}) 
#'   or the number of unique sample names (for \code{collapseBySample=TRUE}) 
#'   in \code{proj}. If \code{NULL}, the wig or bigWig file names are generated 
#'   from the names of the genomic bam files or unique sample names with an 
#'   added \dQuote{.wig.gz} or \dQuote{.bw} extension.
#' @param collapseBySample If \code{TRUE}, genomic bam files with identical
#'   sample name will be combined (summed) into a single track.
#' @param binsize A numerical value defining the bin and step size for the 
#'   wig or bigWig file(s). \code{binsize} will be coerced to \code{integer()}.
#' @param shift Either a vector or a scalar value defining the read shift (e.g.
#'   half of fragment length, see \sQuote{Details}). If \code{length(shift)>1},
#'   the length must match the number of bam files in \sQuote{proj}, and
#'   the i-th sample will be converted to wig or bigWig using the value in
#'   \code{shift[i]}. \code{shift} will be coerced to \code{integer()}. For
#'   paired-end alignments, \code{shift} will be ignored, and a warning
#'   will be issued if it is set to a non-zero value (see \sQuote{Details}).
#' @param strand Only count alignments of \code{strand}. The default 
#'   (\dQuote{*}) will count all alignments.
#' @param scaling If TRUE or a numerical value, the output values in the wig 
#'   or bigWig file(s) will be linearly scaled by the total number of aligned 
#'   reads per sample to improve comparability (see \sQuote{Details}).
#' @param tracknames A character vector with the names of the tracks to appear 
#'   in the track header. If \code{NULL}, the sample names in \code{proj} 
#'   will be used.
#' @param log2p1 If \code{TRUE}, the number of alignments \code{x} per bin will
#'   be transformed using the formula \code{log2(x+1)}.
#' @param colors A character vector with R color names to be used for the tracks.
#' @param includeSecondary If \code{TRUE} (the default), include alignments 
#'   with the secondary bit (0x0100) set in the \code{FLAG}.
#' @param mapqMin Minimal mapping quality of alignments to be included
#'   (mapping quality must be greater than or equal to \code{mapqMin}). 
#'   Valid values are between 0 and 255. The default (0) will include all 
#'   alignments.
#' @param mapqMax Maximal mapping quality of alignments to be included
#'   (mapping quality must be less than or equal to \code{mapqMax}). 
#'   Valid values are between 0 and 255. The default (255) will include all
#'   alignments.
#' @param absIsizeMin For paired-end experiments, minimal absolute insert
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included. Valid 
#'   values are greater than 0 or \code{NULL} (default), which will not 
#'   apply any minimum insert size filtering.
#' @param absIsizeMax For paired-end experiments, maximal absolute insert 
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included. Valid 
#'   values are greater than 0 or \code{NULL} (default), which will not apply 
#'   any maximum insert size filtering.
#' @param createBigWig If \code{TRUE}, first a temporary wig file will be 
#'   created and then converted to BigWig format (file extension \dQuote{.bw}) 
#'   using the \code{\link[rtracklayer]{wigToBigWig}} function from 
#'   package \pkg{rtracklayer}.
#' @param useRead For paired-end experiments, selects the read mate whose 
#'   alignments should be counted, one of:
#'   \itemize{
#'     \item \code{any} (default): count all alignments
#'     \item \code{first} : count only alignments from the first read
#'     \item \code{last} : count only alignments from the last read
#'   }
#'   For single-read alignments, this argument will be ignored. For
#'   paired-end alignments, setting this argument to a value different
#'   from the default (\code{any}) will cause \code{qExportWig} not to
#'   automatically use the mid of fragments, but to treat the selected
#'   read as if it would come from a single-read experiment (see 
#'   \sQuote{Details}).
#' @param pairedAsSingle If \code{TRUE}, treat paired-end data single read 
#'   data, which means that instead of calculating fragment mid-points for 
#'   each read pair, the 5-prime ends of the reads is used. This is for example 
#'   useful when analyzing paired-end DNAse-seq or ATAC-seq data, in which 
#'   the read starts are informative for chromatin accessibility.
#' 
#' @return (invisible) The file name of the generated wig or bigWig file(s).
#' 
#' @export
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' 
#' @seealso 
#' \code{\linkS4class{qProject}}, \code{\link{qAlign}},
#' \code{\link[rtracklayer]{wigToBigWig}}
#' 
#' @keywords utilities
#' 
#' @name qExportWig
#' @aliases qExportWig
#' 
#' @importFrom Rsamtools scanBamHeader
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom grDevices col2rgb colorRampPalette
#' @importFrom rtracklayer wigToBigWig
#' 
#' @examples 
#' # copy example data to current working directory
#' file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
#' 
#' # create alignments
#' sampleFile <- "extdata/samples_chip_single.txt"
#' genomeFile <- "extdata/hg19sub.fa"
#' proj <- qAlign(sampleFile, genomeFile)
#' 
#' # export wiggle file
#' qExportWig(proj, binsize=100L, shift=0L, scaling=TRUE)
#' 
qExportWig <- function(proj,
                       file = NULL,
                       collapseBySample = TRUE,
                       binsize = 100L,
                       shift = 0L,
                       strand = c("*", "+", "-"),
                       scaling = TRUE,
                       tracknames = NULL,
                       log2p1 = FALSE,
                       colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                  "#66A61E", "#E6AB02", "#A6761D", "#666666"),
                       includeSecondary = TRUE,
                       mapqMin = 0L,
                       mapqMax = 255L,
                       absIsizeMin = NULL,
                       absIsizeMax = NULL,
                       createBigWig = FALSE,
                       useRead = c("any", "first", "last"),
                       pairedAsSingle = FALSE) {
    # validate parameters
    # ...proj
    if (!is(proj, "qProject"))
        stop("'proj' must be a 'qProject' object")
    
    if (collapseBySample) {
        bamfiles <- split(
            proj@alignments$FileName, 
            as.factor(proj@alignments$SampleName))[unique(proj@alignments$SampleName)]
        samplenames <- names(bamfiles)
    } else {
        bamfiles <- as.list(proj@alignments$FileName)
        samplenames <- proj@alignments$SampleName
    }
    n <- length(bamfiles)
    paired <- proj@paired != "no"
    
    # ...strand
    strand <- match.arg(strand) # checks if strand has length 1
    
    # ...tracknames
    if (is.null(tracknames)) {
        tracknames <- if (collapseBySample) samplenames else displayNames(proj)
        if (strand[1] != "*")
            tracknames <- sprintf("%s (%s)", tracknames, strand)
    }
    
    # ...file
    if (is.null(file)) {
        fileExt <- if (createBigWig) ".bw" else ".wig.gz"
        if (collapseBySample)
            file <- paste0(samplenames, fileExt)
        else
            file <- paste0(displayNames(proj), fileExt)
    } else if (length(file) != n) {
        stop(sprintf("the length of 'file' (%d) does not match the number of wig files to be generated (%d)", length(file), n))
    }
    if (createBigWig && any(!grepl(".bw$", file)))
        stop("file names have to end with '.bw' for createBigWig=TRUE")
    compressFormat <- compressedFileFormat(file)
    if (!all(compressFormat %in% c("none", "gzip")))
        stop("only gzip compressed wig files (extension '.gz') are supported")
    compress <- compressFormat == "gzip"
    if (length(compress) == 1)
        compress <- rep(compress, n)

    # ...binsize
    if (length(binsize) != 1)
        stop("'binsize' must be a single integer value")
    binsize <- as.integer(binsize)
    if (is.na(binsize) || binsize < 1)
        stop("'binsize' must be a positive integer value")

    # ...shift
    shift <- as.integer(shift)
    if (any(is.na(shift)))
        stop("'shift' has to be a vector of integer values")
    if (length(shift) != 1 && length(shift) != n)
        stop(sprintf("'shift' has to contain either a single value or one value per output wig file (%d)", n))
    if (paired && shift != 0L) {
        warning("ignoring 'shift' value for paired-end alignments (will calculate alignment-specific values)")
        shift <- 0L
    }
    if (length(shift) == 1)
        shift <- rep(shift, n)

    # ...scaling
    if (!(length(scaling) == 1L && (is.logical(scaling) || is.numeric(scaling))))
        stop("'scaling' needs to be a logical(1) or a numeric(1)")
    fact <- rep(1, n)
    if (is.numeric(scaling) || (is.logical(scaling) && scaling)) {
        message("collecting mapping statistics for scaling...", appendLF = FALSE)
        tmp <- alignmentStats(proj, collapseBySample = collapseBySample)
        N <- tmp[grepl(":genome$", rownames(tmp)), 'mapped']
        names(N) <- sub(":genome$", "", names(N))
        if (is.logical(scaling)) {
            #fact <- min(N) / N
            fact <- mean(N) / N
        } else if(is.numeric(scaling) && length(scaling) == 1L && scaling > 0) {
            fact <- scaling / N
        } else {
            stop("'scaling' must be greater than zero")
        }
        message("done")
    }

    # ...log2p1
    if (length(log2p1) != 1 || !is.logical(log2p1))
        stop("'log2p1' has to be either 'TRUE' or 'FALSE'")
    log2p1 <- rep(log2p1,n)
    
    # ...colors
    if (length(colors) < n)
        colors <- grDevices::colorRampPalette(colors)(n)
    colors <- apply(grDevices::col2rgb(colors), 2, paste, collapse = ",")

    # ...includeSecondary
    if (length(includeSecondary) != 1 || !is.logical(includeSecondary))
        stop("'includeSecondary' must be of type logical(1)")

    # ...mapping qualities
    if (length(mapqMin) != 1 || !is.integer(mapqMin) || 
        any(is.na(mapqMin)) || min(mapqMin) < 0L || max(mapqMax) > 255L)
        stop("'mapqMin' must be of type integer(1) and have a values between 0 and 255")
    mapqMin <- rep(mapqMin, n)
    if (length(mapqMax) != 1 || !is.integer(mapqMax) || 
        any(is.na(mapqMax)) || min(mapqMax) < 0L || max(mapqMax) > 255L)
        stop("'mapqMax' must be of type integer(1) and have a values between 0 and 255")
    mapqMax <- rep(mapqMax, n)

    # ...absolute insert size
    if ((!is.null(absIsizeMin) || !is.null(absIsizeMax)) && !paired)
        stop("'absIsizeMin' and 'absIsizeMax' can only be used for paired-end experiments")
    if (is.null(absIsizeMin)) # -1L -> do not apply TLEN filtering
        absIsizeMin <- -1L
    if (is.null(absIsizeMax))
        absIsizeMax <- -1L

    # ...useRead
    useRead <- match.arg(useRead)
    if (useRead != "any" && !paired) {
        warning("ignoring 'useRead' for single read experiments")
        useRead <- "any"
    } else if (paired && useRead != "any") {
        message("'useRead' is set - will treat alignments as single reads (no calculation of fragment midpoints)")
        paired <- FALSE
    } else if (pairedAsSingle) {
        if (paired) {
            message("'pairedAsSingle' is set - will treat alignments as single reads (no calculation of fragment midpoints)")
            paired <- FALSE
        } else {
            warning("ignoring 'pairedAsSingle' for single read experiments")
        }
    }

    # translate useRead parameter
    BAM_FREAD1 <- 64L
    BAM_FREAD2 <- 128L
    if (useRead == "any") {
        readBitMask <- BAM_FREAD1 + BAM_FREAD2
    } else if (useRead == "first") {
        readBitMask <- BAM_FREAD1
    } else if (useRead == "last") {
        readBitMask <- BAM_FREAD2
    }
    
    # generate the wig file(s)
    message("start creating ", if (createBigWig) "bigWig" else "wig"," file", 
            if (n > 1) "s" else "", "...")
    tempwigfile <- if (createBigWig) sapply(1:n, function(i) tempfile(fileext = ".wig")) else file
    if (createBigWig) {
        tmp <- Rsamtools::scanBamHeader(bamfiles[[1]][1])[[1]]$targets
        si <- GenomeInfoDb::Seqinfo(names(tmp), tmp)
    }
    lapply(1:n, function(i) {
        message("  ",file[i]," (",tracknames[i],")")
        .Call(bamfileToWig, as.character(bamfiles[[i]]), 
              as.character(tempwigfile[i]), as.logical(paired[1]),
              as.integer(binsize[1]), as.integer(shift[i]), 
              as.character(strand[1]), as.numeric(fact[i]),
              as.character(tracknames[i]), as.logical(log2p1[i]),
              as.character(colors[i]), as.logical(compress[i]), 
              as.logical(includeSecondary[1]),
              mapqMin[i], mapqMax[i], as.integer(absIsizeMin), 
              as.integer(absIsizeMax), readBitMask)
        if (createBigWig) {
            rtracklayer::wigToBigWig(tempwigfile[i], si, file[i])
            unlink(tempwigfile[i])
        }
    })
    message("done")
    
    return(invisible(file))
}
