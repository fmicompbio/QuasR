# proj       : qProject object
# query      : NULL (whole genome)
#              GRanges object (only quantify C's in these regions)
# reportLevel: "C" or "alignment"
# mode       : "allC"            : all C's (+/- strands separate)
#              "CpG"             : only C's in CpG context (+/- strands separate)
#              "CpGcomb"(default): only C's in CpG context (+/- strands collapsed)
#              "var"             : variant detection (all C's, +/- strands separate)
# collapseBySample : combine (sum) counts from bamfiles with the same sample name
# collapseByQueryRegion : combine (sum) counts for C's per query region
# asGRanges  : return value as GRanges object or data.frame
# mask       : mask genomic regions (e.g. unmappable regions)
# reference  : source of bam files ("genome" or AuxName)
# keepZero   : return C's with total==0 in results
# mapqMin    : minimal mapping quality
# mapqMax    : maximal mapping quality
# clObj      : cluster object for parallelization
#
# value      : if asGRanges==TRUE, GRanges object with one region per quantified C and two metadata columns per sample (_T, _M)
#              else, data.frame with one row per quantified C and in the columns the coordinates of the C, as well as two count (_T, _M) per sample

# TODO:
# - support for gapped aligments (parsing of cigar strings at C level)
# - ignore presumable PCR duplicates (multiple alignment pairs with same external coordinates)
# - ignore overlapping read pairs (insert size smaller than twice the read length)

#' Quantify DNA methylation
#' 
#' Quantify methylation of cytosines from bisulfite sequencing data.
#' 
#' \code{qMeth} can be used on a \code{qProject} object from a bisulfite
#' sequencing experiment (sequencing of bisulfite-converted DNA), such as
#' the one returned by \code{\link[QuasR]{qAlign}} when its parameter
#' \code{bisulfite} is set to a different value than \dQuote{no}.
#' 
#' \code{qMeth} quantifies DNA methylation by counting total and
#' methylated events for individual cytosines, using the alignments that
#' have been generated in converted (three-letter) sequence space for
#' example by \code{\link[QuasR]{qAlign}}. A methylated event corresponds 
#' to a C/C match in the alignment, an unmethylated event to a T/C mismatch
#' (or G/G matches and A/G mismatches on the opposite strand). For paired-end
#' samples, the part of the left fragment alignment that overlaps
#' with the right fragment alignment is ignored, preventing the
#' use of redundant information coming from the same molecule.
#' 
#' Both directed (\code{bisulfite}=\dQuote{dir}) and undirected
#' (\code{bisulfite}=\dQuote{undir}) experimental protocols are supported
#' by \code{\link[QuasR]{qAlign}} and \code{qMeth}.
#' 
#' By default, results are returned per C nucleotide. If
#' \code{reportLevel}=\dQuote{alignment}, results are reported separately
#' for individual alignments. In that case, \code{query} has to be a
#' \code{GRanges} object with exactly one region, \code{mode} has to be
#' either \dQuote{CpG} or \dQuote{allC}, the arguments
#' \code{collapseByQueryRegion}, \code{asGRanges}, \code{mask} and
#' \code{keepZero} have no effect and allele-specific projects are
#' treated in the same way as normal (non-allele specific) projects.
#' 
#' Using the parameter \code{mode}, quantification can be limited to
#' cytosines in CpG context, and counts obtained for the two cytosines on
#' opposite strands within a single CpG can be combined (summed).
#' 
#' The quantification of methylation for all cytosines in the query region(s)
#' (\code{mode}=\dQuote{allC}) should be done with care, especially for
#' large query regions, as the return value may require a large amount of
#' memory.
#' 
#' If \code{mode} is set to \dQuote{var}, \code{qMeth} only counts reads
#' from the strand opposite of the cytosine and reports total and
#' matching alignments. For a position identical to the reference
#' sequence, only matches (and very few sequencing errors) are
#' expected, independent on the methylation state of the cytosine. A
#' reduced fraction of alignments matching the reference are indicative
#' of sequence variations in the sequenced sample.
#' 
#' \code{mapqMin} and \code{mapqMax} allow to select alignments
#' based on their mapping qualities. \code{mapqMin} and \code{mapqMax} can
#' take integer values between 0 and 255 and equal to
#' \eqn{-10 log_{10} Pr(\textnormal{mapping position is wrong})}{-10
#' log10 Pr(mapping position is wrong)}, rounded to the nearest
#' integer. A value 255 indicates that the mapping quality is not available.
#' 
#' If an object that inherits from class \code{cluster} is provided to
#' the \code{clObj} argument, for example an object returned by
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel},
#' the quantification task is split into multiple chunks and processed in
#' parallel using \code{\link[parallel:clusterApply]{clusterApplyLB}} from package
#' \pkg{parallel}. Not all tasks will be efficiently parallelized: For
#' example, a single query region and a single (group of) bam files will
#' not be split into multiple chunks.
#' 
#' @param proj A \code{qProject} object from a bisulfite sequencing experiment
#' @param query A \code{GRanges} object with the regions to be
#'   quantified. If \code{NULL}, all available target sequences (e.g. the
#'   whole genome) will be analyzed. Available target sequences are
#'   extracted from the header of the first bam file.
#' @param reportLevel Report results combined for C's
#'   (\code{reportLevel}=\dQuote{C}), the default) or individually for single
#'   alignments (\code{reportLevel}=\dQuote{alignment}). The latter imposes
#'   further restrictions on some arguments (see \sQuote{Details}).
#' @param mode Cytosine quantification mode, one of:
#'   \itemize{
#'     \item \code{CpGcomb} : only C's in CpG context (strands combined)
#'     \item \code{CpG} : only C's in CpG context (strands separate)
#'     \item \code{allC} : all C's (strands separate)
#'     \item \code{var} : variant detection (all C's, strands separate)
#'   }
#'   \code{CpGcomb} is the default.
#' @param collapseBySample If \code{TRUE}, combine (sum) counts from
#'   bamfiles with the same sample name.
#' @param collapseByQueryRegion If \code{TRUE}, combine (sum) counts for
#'   all cytosines contained in the same query region.
#' @param asGRanges If \code{TRUE}, return results as a \code{GRanges} object;
#'   if \code{FALSE}, the results are returned as a \code{data.frame}.
#' @param mask An optional \code{GRanges} object with genomic regions to
#'   be masked, i.e. excluded from the analysis (e.g. unmappable regions).
#' @param reference Source of bam files; can be either \dQuote{genome}
#'   (then the alignments against the genome are used) or the name of an
#'   auxiliary target sequence (then alignments against this target
#'   sequence will be used). The auxiliary name must correspond to the
#'   name contained in the auxiliary file refered by the
#'   \code{auxiliaryFile} argument of \code{\link[QuasR]{qAlign}}.
#' @param keepZero If \code{FALSE}, only cytosines covered by at least
#'   one alignment will be returned; \code{keepZero} must be \code{TRUE}
#'   if multiple samples have the same sample name and
#'   \code{collapseBySample} is \code{TRUE}.
#' @param mapqMin Minimal mapping quality of alignments to be included when
#'   counting (mapping quality must be greater than or equal to
#'   \code{mapqMin}). Valid values are between 0 and 255. The default (0)
#'   will include all alignments.
#' @param mapqMax Maximal mapping quality of alignments to be included when
#'   counting (mapping quality must be less than or equal to \code{mapqMax}).
#'   Valid values are between 0 and 255. The default (255) will include
#'   all alignments.
#' @param clObj A cluster object to be used for parallel processing of
#'   multiple files (see \sQuote{Details}).
#' 
#' @export
#' 
#' @return 
#' For \code{reportLevel}=\dQuote{C}, a \code{GRanges} object if
#' \code{asGRanges}=\code{TRUE}, otherwise a \code{data.frame}.
#' 
#' Each row contains the coordinates of individual cytosines for
#' \code{collapseByQueryRegion}=\code{FALSE} or query regions
#' for \code{collapseByQueryRegion}=\code{TRUE}.
#' 
#' In addition to the coordinates columns (or \code{seqnames},
#' \code{ranges} and \code{strand} slots for \code{GRanges} objects),
#' each row contains per bam file:
#' 
#' Two values (total and methylated events, with suffixes _T and _M), or
#' if the \code{qProject} object was created including a SNP table,
#' six values (total and methylated events for Reference, Unknown and
#' Alternative genotypes, with suffixed _TR, _TU, _TA, _MR, _MU and _MA).
#' In the latter case, C's or CpG's that overlap with SNPs in the table
#' are removed.
#' 
#' If \code{collapseBySample}=\code{TRUE}, groups of bam files with
#' identical sample name are combined (summed) and will be represented by
#' a single set of total and methylated count columns.
#' 
#' If \code{mode}=\dQuote{var}, the _T and _M columns correspond to total
#' and matching alignments overlapping the guanine paired to the cytosine.
#' 
#' For \code{reportLevel}=\dQuote{alignment}, a \code{list} with one
#' element per bam file or sample (depending on \code{collapseBySample}). 
#' Each list element is another list with the elements:
#' \itemize{
#'   \item \code{aid}: character vector with unique alignment identifiers
#'   \item \code{Cid}: integer vector with genomic coordinate of C base
#'   \item \code{strand}: character vector with the strand of the C base
#'   \item \code{meth}: integer vector with methylation state for
#'   alignment and C defined by \code{aid} and \code{Cid}. The values are
#'   1 for methylated or 0 for unmethylated states.
#' }
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' @keywords utilities misc
#' 
#' @seealso 
#' \code{\link[QuasR]{qAlign}},
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
#' 
#' @name qMeth
#' @aliases qMeth
#' 
#' @examples 
#' # copy example data to current working directory
#' file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
#' 
#' # create alignments
#' sampleFile <- "extdata/samples_bis_single.txt"
#' genomeFile <- "extdata/hg19sub.fa"
#' proj <- qAlign(sampleFile, genomeFile, bisulfite="dir")
#' proj
#' 
#' # calculate methylation states
#' meth <- qMeth(proj, mode="CpGcomb")
#' meth
#' 
#' @importFrom Rsamtools scanBamHeader scanFaIndex
#' @importFrom parallel clusterEvalQ clusterApplyLB
#' @importFrom GenomeInfoDb seqlevels seqinfo
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
qMeth <- function(proj,
                  query = NULL,
                  reportLevel = c("C", "alignment"),
                  mode = c("CpGcomb", "CpG", "allC", "var"),
                  collapseBySample = TRUE,
                  collapseByQueryRegion = FALSE,
                  asGRanges = TRUE,
                  mask = NULL,
                  reference = "genome",
                  keepZero = TRUE,
                  mapqMin = 0L,
                  mapqMax = 255L,
                  clObj = NULL) {
    ## setup variables from 'proj' ---------------------------------------------
    # 'proj' is correct type?
    if (!inherits(proj, "qProject", which = FALSE))
        stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")
    if (proj@bisulfite == "no")
        stop("'proj' is not a bisufite-seq project")
    if (proj@splicedAlignment)
        stop("'spliceAlignment==TRUE' is not supported by qMeth")
    if (proj@aligner == "Rhisat2")
        stop("Rhisat2 is not supported by qMeth")
    
    samples <- proj@alignments$SampleName
    nsamples <- length(samples)
    if (reference == "genome") {
        bamfiles <- proj@alignments$FileName
        referenceFormat <- proj@genomeFormat
        referenceSource <- proj@genome
    } else if (!is.na(i <- match(reference, rownames(proj@auxAlignments)))) {
        bamfiles <- unlist(proj@auxAlignments[i, ], use.names = FALSE)
        referenceFormat <- "file"
        referenceSource <- proj@aux[i, 'FileName']
    } else {
        stop("unknown 'reference', should be one of: ",
             paste(sprintf("'%s'", c("genome", rownames(proj@auxAlignments))), collapse = ", "))
    }
    reportLevel <- match.arg(reportLevel)
    mode <- match.arg(mode)

    
    ## validate parameters -----------------------------------------------------
    # 'query' is correct type?
    if (is.null(query)) {
        if (reportLevel == "alignment")
            stop("'query' must be an object of type 'GRanges' with length 1 for reportLevel='alignment'")
        tr <- Rsamtools::scanBamHeader(bamfiles[1])[[1]]$targets
        query <- GenomicRanges::GRanges(names(tr), 
                                        IRanges::IRanges(start = 1, end = tr))
    } else if(!inherits(query, "GRanges")) {
        stop("'query' must be either NULL or an object of type 'GRanges'")
    }

    if (!is.logical(keepZero) || length(keepZero) != 1)
        stop("'keepZero' must be either TRUE or FALSE")
    if (!is.logical(asGRanges) || length(asGRanges) != 1)
        stop("'asGRanges' must be either TRUE or FALSE")
    
    if (!keepZero && ((collapseBySample && length(unique(samples)) > 1) || 
                      (!collapseBySample && length(samples) > 1)))
        stop("'keepZero' must be TRUE if there are multiple non-collapsable samples")
    
    if (mode == "var" && collapseByQueryRegion)
        stop("'collapseByQueryRegion' must be FALSE for variant detection mode")
    if (mode == "var" && !is.na(proj@snpFile))
        stop("allele-specific mode cannot be combined with variant detection mode")
    
    if (reportLevel == "alignment") {
        if (!(mode %in% c("CpG", "allC")))
            stop("'mode' must be 'CpG' or 'allC' for reportLevel='alignment'")
        if (length(query) != 1)
            stop("'query' must be an object of type 'GRanges' with length 1 for reportLevel='alignment'")
    }
    
    # all query chromosomes present in all bamfiles?
    trTab <- table(unlist(lapply(Rsamtools::scanBamHeader(bamfiles), 
                                 function(bh) names(bh$targets))))
    trCommon <- names(trTab)[trTab == length(bamfiles)]
    if (any(f <- !(GenomeInfoDb::seqlevels(query) %in% trCommon)))
        stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                     paste(GenomeInfoDb::seqlevels(query)[f], collapse = ", ")))
    
    
    ## apply 'mask' to query ---------------------------------------------------
    if (!is.null(mask)) {
        if (reportLevel == "alignment")
            warning("ignoring 'mask' for reportLevel='alignment'")
        stop("'mask' (masking of query regions, e.g. unmappable genomic ", 
             "regions) is not implemented yet")
    }
    
    ## setup tasks for parallelization -----------------------------------------
    ## TODO: create several chunks per chromosome?
    taskIByQuery <- split(seq.int(length(query)), 
                          as.factor(GenomeInfoDb::seqnames(query)))
    taskIByQuery <- taskIByQuery[sapply(taskIByQuery,length) > 0]
    nChunkQuery <- length(taskIByQuery)
    
    if (collapseBySample) {
        taskBamfiles <- split(bamfiles, samples)
        sampleNames <- names(taskBamfiles)
    } else {
        taskBamfiles <- as.list(bamfiles)
        sampleNames <- displayNames(proj)
    }
    nChunkBamfile <- length(taskBamfiles)
    
    nChunk <- nChunkQuery * nChunkBamfile
    taskIByQuery <- rep(taskIByQuery, nChunkBamfile)
    taskBamfiles <- rep(taskBamfiles, each = nChunkQuery)

    if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
        message("preparing to run on ", min(nChunk, length(clObj)), " nodes...", 
                appendLF = FALSE)
        ret <- parallel::clusterEvalQ(clObj, library("QuasR")) # load package on nodes
        if (!all(sapply(ret, function(x) "QuasR" %in% x)))
            stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
        myapply <- function(...)
            parallel::clusterApplyLB(clObj, ...)
        message("done")
    } else {
        myapply <- function(...)
            lapply(...)
    }
    
    ## quantify methylation  ---------------------------------------------------
    if (reportLevel == "alignment") {
        # ...per alignment reporting mode: list(nSamples) of list(3) with 
        # "aid","Cid","meth" elements
        resL <- myapply(
            seq_len(nChunk), #nChunk is always equal to nChunkBamfile (length(taskBamfiles))
            function(i) quantifyMethylationBamfilesRegionsSingleChromosomeSingleAlignments(
                taskBamfiles[[i]],
                query[taskIByQuery[[i]]],
                collapseByQueryRegion,
                mode,
                referenceFormat,
                referenceSource,
                as.integer(mapqMin)[1],
                as.integer(mapqMax)[1])
        )
        names(resL) <- sampleNames
        res <- resL
        
    } else if (!is.na(proj@snpFile)) {
        # allele-bis-mode: 6 columns per sample in output (TR, MR, TU, MU, TA, MA) -------------
        resL <- myapply(
            seq_len(nChunk),
            function(i) quantifyMethylationBamfilesRegionsSingleChromosomeAllele(
                taskBamfiles[[i]],
                query[taskIByQuery[[i]]],
                collapseByQueryRegion,
                mode,
                referenceFormat,
                referenceSource,
                proj@snpFile,
                keepZero,
                as.integer(mapqMin)[1],
                as.integer(mapqMax)[1])
        )
        # combine and reorder chunks
        # ...rbind chunks for the same sample
        resL <- lapply(split(seq_len(nChunk), rep(seq_len(nChunkBamfile), each = nChunkQuery)),
                       function(i) do.call(rbind, resL[i]))
        names(resL) <- sampleNames
        
        if (any(unlist(lapply(resL, nrow), use.names = FALSE) != nrow(resL[[1]])))
            stop("error while combining partial results (chunks are incompatable)")
        
        # ...cbind TR/MR/TU/MU/TA/MA columns for different samples
        res <- cbind(resL[[1]][, c("chr", "start", "end", "strand")],
                     do.call(cbind, lapply(resL, "[", c("TR", "MR", "TU", "MU", "TA", "MA"))),
                     stringsAsFactors = FALSE)
        colnames(res)[5:ncol(res)] <- sprintf("%s_%s", rep(sampleNames, each = 6),
                                              c("TR", "MR", "TU", "MU", "TA", "MA"))
    } else if (mode == "var") {
        # variant detection mode: 2 columns per sample in output (total and match) -------------
        resL <- myapply(
            seq_len(nChunk),
            function(i) detectVariantsBamfilesRegionsSingleChromosome(
                taskBamfiles[[i]],
                query[taskIByQuery[[i]]],
                referenceFormat,
                referenceSource,
                keepZero,
                as.integer(mapqMin)[1],
                as.integer(mapqMax)[1])
        )
        # combine and reorder chunks
        # ...rbind chunks for the same sample
        resL <- lapply(split(seq_len(nChunk), rep(seq_len(nChunkBamfile), each = nChunkQuery)),
                       function(i) do.call(rbind, resL[i]))
        names(resL) <- sampleNames
        
        if (any(unlist(lapply(resL, nrow), use.names = FALSE) != nrow(resL[[1]])))
            stop("error while combining partial results (chunks are incompatable)")
        
        # ...cbind T/M columns for different samples
        res <- cbind(resL[[1]][, c("chr", "start", "end", "strand")],
                     do.call(cbind, lapply(resL, "[", c("T", "M"))),
                     stringsAsFactors = FALSE)
        colnames(res)[5:ncol(res)] <- sprintf("%s_%s", rep(sampleNames, each = 2), c("T", "M"))
    } else {
        # normal mode: 2 columns per sample in output (T and M) -------------
        resL <- myapply(
            seq_len(nChunk),
            function(i) quantifyMethylationBamfilesRegionsSingleChromosome(
                taskBamfiles[[i]],
                query[taskIByQuery[[i]]],
                collapseByQueryRegion,
                mode,
                referenceFormat,
                referenceSource,
                keepZero,
                as.integer(mapqMin)[1],
                as.integer(mapqMax)[1])
        )
        # combine and reorder chunks
        # ...rbind chunks for the same sample
        resL <- lapply(split(seq_len(nChunk), rep(seq_len(nChunkBamfile), each = nChunkQuery)),
                       function(i) do.call(rbind, resL[i]))
        names(resL) <- sampleNames
        
        if (any(unlist(lapply(resL, nrow), use.names = FALSE) != nrow(resL[[1]])))
            stop("error while combining partial results (chunks are incompatable)")
        
        # ...cbind T/M columns for different samples
        res <- cbind(resL[[1]][, c("chr", "start", "end", "strand")],
                     do.call(cbind, lapply(resL, "[", c("T", "M"))),
                     stringsAsFactors = FALSE)
        colnames(res)[5:ncol(res)] <- sprintf("%s_%s", rep(sampleNames, each = 2), c("T", "M"))
    }
    
    if (asGRanges && reportLevel != "alignment") {
        if (referenceFormat == "file") {
            si <- GenomeInfoDb::seqinfo(Rsamtools::scanFaIndex(referenceSource))
        } else {
            library(referenceSource, character.only = TRUE)
            gnmObj <- get(referenceSource)
            si <- GenomeInfoDb::seqinfo(gnmObj)
        }
        res <- GenomicRanges::GRanges(seqnames = res$chr, 
                                      IRanges::IRanges(start = res$start, 
                                                       end = res$end),
                       strand = res$strand, seqinfo = si, res[, 5:ncol(res)])
    }
    
    ## return results
    return(res)
}


# detect variants for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, never collapsed)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+2*nSamples vectors: chr, start, end, strand of C, counts of T (total) and M (match) reads
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom Rsamtools scanFaIndex scanFa
#' @importFrom BSgenome getSeq
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start end
#' 
detectVariantsBamfilesRegionsSingleChromosome <- function(bamfiles, regions, 
                                                          referenceFormat, 
                                                          reference, keepZero, 
                                                          mapqmin, mapqmax) {
    ## verify parameters
    if (length(chr <- as.character(unique(GenomeInfoDb::seqnames(regions)))) != 1)
        stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")
    
    ## collapse regions
    regionsStart <- as.integer(min(BiocGenerics::start(regions)))
    regionsEnd   <- as.integer(max(BiocGenerics::end(regions)))
    regionsGr    <- GenomicRanges::GRanges(
        chr, IRanges::IRanges(start = regionsStart, end = regionsEnd)
    )
    
    ## get sequence string from...
    #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
    if (referenceFormat == "file") { # genome file
        chrLen <- as.integer(GenomeInfoDb::seqlengths(Rsamtools::scanFaIndex(reference))[chr])
        seqstr <- as.character(Rsamtools::scanFa(reference, regionsGr)[[1]])
    } else {                        # BSgenome object
        library(reference, character.only = TRUE)
        referenceObj <- get(reference) # access the BSgenome
        chrLen <- as.integer(length(referenceObj[[chr]]))
        seqstr <- BSgenome::getSeq(referenceObj, regionsGr, as.character = TRUE)
    }
    #message("done")
    
    ## call CPP function (multiple bam files, single region)
    #message("detecting single nucleotide variations...", appendLF=FALSE)
    resL <- .Call(detectSNVs, bamfiles, chr, chrLen, regionsStart, seqstr,
                  keepZero, mapqmin, mapqmax)
    #message("done")
    
    ## filter out C's that do not fall into 'regions'
    #message("processing results...", appendLF=FALSE)
    ov <- GenomicRanges::findOverlaps(
        query = GenomicRanges::GRanges(
            chr, IRanges::IRanges(start = resL$position, width = 1
            )
        ), 
        subject = regions)
    resL <- lapply(resL, "[", unique(S4Vectors::queryHits(ov)))
    
    res <- data.frame(chr = resL$chr,
                      start = resL$position,
                      end = resL$position,
                      strand = rep("*", length(resL$chr)),
                      T = resL$nTotal,
                      M = resL$nMatch,
                      stringsAsFactors = FALSE)
    #message("done")
    
    res
}

# quantify methylation for (report for INDIVIDUAL ALIGMENTS):
#  - multiple bamfiles (will allways be collapsed)
#  - a single region
#  - mode (defines which and how C's are quantified)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a list(4) with elements "aid","Cid","strand","meth"
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom Rsamtools scanFaIndex scanFa
#' @importFrom BSgenome getSeq
#' @importFrom BiocGenerics start end
#' 
quantifyMethylationBamfilesRegionsSingleChromosomeSingleAlignments <-
    function(bamfiles, regions, collapseByQueryRegion, mode = c("CpG","allC"), 
             referenceFormat, reference, mapqmin, mapqmax) {
        ## verify parameters
        if (length(regions) != 1)
            stop("'regions' must be of length 1 for 'quantifyMethylationBamfilesRegionsSingleChromosomeSingleAlignments'")
        mode <- c("CpG" = 1L, "allC" = 2L)[match.arg(mode)]
        chr <- as.character(unique(GenomeInfoDb::seqnames(regions)))
        regionsStart <- as.integer(BiocGenerics::start(regions))
        regionsEnd   <- as.integer(BiocGenerics::end(regions))
        regionsGr    <- GenomicRanges::GRanges(
            chr, IRanges::IRanges(start = regionsStart, end = regionsEnd)
        )
        
        ## get sequence string from...
        if (referenceFormat == "file") { # genome file
            chrLen <- as.integer(GenomeInfoDb::seqlengths(
                Rsamtools::scanFaIndex(reference))[chr]
            )
            seqstr <- as.character(Rsamtools::scanFa(reference, regionsGr)[[1]])
        } else {                         # BSgenome object
            library(reference, character.only = TRUE)
            referenceObj <- get(reference) # access the BSgenome
            chrLen <- as.integer(length(referenceObj[[chr]]))
            seqstr <- BSgenome::getSeq(referenceObj, regionsGr, as.character = TRUE)
        }
        
        ## call CPP function (multiple bam files, single region)
        resL <- .Call(quantifyMethylationSingleAlignments, bamfiles, chr, chrLen,
                      regionsStart, seqstr, mode, mapqmin, mapqmax)
        
        return(resL)
    }

# quantify methylation for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, may be collapsed if collapseByRegion==TRUE)
#  - mode (defines which and how C's are quantified)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+2*nSamples vectors: chr, start, end, strand of C, counts of T (total) and M (methylated) reads
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom Rsamtools scanFaIndex scanFa
#' @importFrom BSgenome getSeq
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics start end
quantifyMethylationBamfilesRegionsSingleChromosome <- function(bamfiles, regions, 
                                                               collapseByRegion, 
                                                               mode = c("CpGcomb", "CpG", "allC"),
                                                               referenceFormat, reference,
                                                               keepZero, mapqmin, mapqmax) {
    ## verify parameters
    if (length(chr <- as.character(unique(GenomeInfoDb::seqnames(regions)))) != 1)
        stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")
    mode <- c("CpGcomb" = 0L, "CpG" = 1L, "allC" = 2L)[match.arg(mode)]
    Cwidth <- ifelse(mode == 0, 2L, 1L)
    
    ## collapse regions
    regionsStart <- as.integer(min(BiocGenerics::start(regions)))
    regionsEnd   <- as.integer(max(BiocGenerics::end(regions)))
    regionsGr    <- GenomicRanges::GRanges(
        chr, IRanges::IRanges(start = regionsStart, end = regionsEnd)
    )
    
    ## get sequence string from...
    #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
    if (referenceFormat == "file") { # genome file
        chrLen <- as.integer(GenomeInfoDb::seqlengths(
            Rsamtools::scanFaIndex(reference))[chr]
        )
        seqstr <- as.character(Rsamtools::scanFa(reference, regionsGr)[[1]])
    } else {                         # BSgenome object
        library(reference, character.only = TRUE)
        referenceObj <- get(reference) # access the BSgenome
        chrLen <- as.integer(length(referenceObj[[chr]]))
        seqstr <- BSgenome::getSeq(referenceObj, regionsGr, as.character = TRUE)
    }
    #message("done")
    
    ## call CPP function (multiple bam files, single region)
    #message("quantifying methylation...", appendLF=FALSE)
    resL <- .Call(quantifyMethylation, bamfiles, chr, chrLen, regionsStart,
                  seqstr, mode, keepZero, mapqmin, mapqmax)
    #message("done")
    
    ## collapse by region
    #message("processing results...", appendLF=FALSE)
    ov <- GenomicRanges::findOverlaps(
        query = GenomicRanges::GRanges(chr, IRanges::IRanges(start = resL$position,
                                                             width = Cwidth)),
        subject = regions)
    
    if (collapseByRegion) {
        tmpT <- tmpM <- numeric(length(regions))
        tmp <- tapply(resL$T[S4Vectors::queryHits(ov)], S4Vectors::subjectHits(ov), sum)
        tmpT[as.integer(names(tmp))] <- tmp
        tmp <- tapply(resL$M[S4Vectors::queryHits(ov)], S4Vectors::subjectHits(ov), sum)
        tmpM[as.integer(names(tmp))] <- tmp
        res <- data.frame(chr = rep(chr, length(regions)),
                          start = start(regions),
                          end = end(regions),
                          strand = as.character(strand(regions)),
                          T = tmpT,
                          M = tmpM,
                          stringsAsFactors = FALSE)
    } else {
        ## filter out C's that do not fall into 'regions'
        resL <- lapply(resL, "[", unique(S4Vectors::queryHits(ov)))
        
        res <- data.frame(chr = resL$chr,
                          start = resL$position,
                          end = resL$position+Cwidth-1L,
                          strand = resL$strand,
                          T = resL$T,
                          M = resL$M,
                          stringsAsFactors = FALSE)
    }
    #message("done")
    
    res
}


# quantify ALLELE-SPECIFIC methylation for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, may be collapsed if collapseByRegion==TRUE)
#  - mode (defines which and how C's are quantified)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+6*nSamples vectors: chr, start, end, strand of C, counts of TR, TU, TA (total) and MR, MU, MA (methylated) reads
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom Rsamtools scanFaIndex scanFa
#' @importFrom BSgenome getSeq
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start end strand
quantifyMethylationBamfilesRegionsSingleChromosomeAllele <-
    function(bamfiles, regions, collapseByRegion, mode = c("CpGcomb", "CpG", "allC"),
             referenceFormat, reference, snpFile, keepZero, mapqmin, mapqmax) {
        ## verify parameters
        if (length(chr <- as.character(unique(GenomeInfoDb::seqnames(regions)))) != 1)
            stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")
        mode <- c("CpGcomb" = 0L, "CpG" = 1L, "allC" = 2L)[match.arg(mode)]
        Cwidth <- ifelse(mode == 0, 2L, 1L)
        
        ## collapse regions
        regionsStart <- as.integer(min(BiocGenerics::start(regions)))
        regionsEnd   <- as.integer(max(BiocGenerics::end(regions)))
        regionsGr    <- GenomicRanges::GRanges(
            chr, IRanges::IRanges(start = regionsStart, end = regionsEnd)
        )
        
        ## get sequence string from...
        #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
        if (referenceFormat == "file") { # genome file
            chrLen <- as.integer(GenomeInfoDb::seqlengths(
                Rsamtools::scanFaIndex(reference))[chr]
            )
            seqstr <- as.character(Rsamtools::scanFa(reference, regionsGr)[[1]])
        } else {                         # BSgenome object
            library(reference, character.only = TRUE)
            referenceObj <- get(reference) # access the BSgenome
            chrLen <- as.integer(length(referenceObj[[chr]]))
            seqstr <- BSgenome::getSeq(referenceObj, regionsGr, as.character = TRUE)
        }
        #message("done")
        
        ## call CPP function (multiple bam files, single region)
        #message("quantifying methylation...", appendLF=FALSE)
        resL <- .Call(quantifyMethylationAllele, bamfiles, chr, chrLen, regionsStart,
                      seqstr, mode, keepZero, mapqmin, mapqmax)
        #message("done")

        ## filter out CpGs that overlap SNPs (may not be possible to discriminate allele from methylation status)
        #message("removing C's overlapping SNPs...", appendLF=FALSE)
        snpL <- scan(snpFile, what = list(chr = "", pos = 1L, R = "", A = ""), quiet = TRUE)
        snp <- GenomicRanges::GRanges(
            snpL$chr, IRanges::IRanges(start = snpL$pos, width = nchar(snpL$R))
        )
        ikeep <- !IRanges::overlapsAny(
            GenomicRanges::GRanges(
                chr, IRanges::IRanges(start = resL$position, width = Cwidth)
            ), snp
        )
        resL <- lapply(resL, "[", ikeep)
        #message(sprintf("removed %d C's, done",sum(!ikeep)))

        ## collapse by region
        #message("processing results...", appendLF=FALSE)
        ov <- GenomicRanges::findOverlaps(
            query = GenomicRanges::GRanges(chr, IRanges::IRanges(start = resL$position,
                                                                 width = Cwidth)), 
            subject = regions
        )

        if (collapseByRegion) {
            tmpTR <- tmpTU <- tmpTA <- tmpMR <- tmpMU <- tmpMA <- numeric(length(regions))
            tmp <- tapply(resL$TR[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpTR[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$TU[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpTU[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$TA[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpTA[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MR[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpMR[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MU[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpMU[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MA[S4Vectors::queryHits(ov)], 
                          S4Vectors::subjectHits(ov), sum)
            tmpMA[as.integer(names(tmp))] <- tmp
            res <- data.frame(chr = rep(chr,length(regions)),
                              start = BiocGenerics::start(regions),
                              end = BiocGenerics::end(regions),
                              strand = as.character(BiocGenerics::strand(regions)),
                              TR = tmpTR,
                              MR = tmpMR,
                              TU = tmpTU,
                              MU = tmpMU,
                              TA = tmpTA,
                              MA = tmpMA,
                              stringsAsFactors = FALSE)
        } else {
            ## filter out C's that do not fall into 'regions'
            resL <- lapply(resL, "[", unique(S4Vectors::queryHits(ov)))

            res <- data.frame(chr = resL$chr,
                              start = resL$position,
                              end = resL$position + Cwidth - 1L,
                              strand = resL$strand,
                              TR = resL$TR,
                              MR = resL$MR,
                              TU = resL$TU,
                              MU = resL$MU,
                              TA = resL$TA,
                              MA = resL$MA,
                              stringsAsFactors = FALSE)
        }
        #message("done")
        
        res
    }

