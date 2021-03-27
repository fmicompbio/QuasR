## query : GRanges : value has one row per unique name in GRanges
## query : GRangesList : value has one row per name of GRangesList
## query : TxDb : value has has one row per 'reportLevel'
## query : NULL & reportLevel=="junction": value has has one row per 'junction'

#' Quantify alignments
#' 
#' Quantify alignments from sequencing data.
#' 
#' \code{qCount} is used to count alignments in each sample from a
#' \code{qProject} object. The features to be quantified, together with
#' the mode of quantification, are specified by the \code{query}
#' argument, which is one of:
#' \itemize{
#'   \item \code{\link[GenomicRanges:GRanges-class]{GRanges}}: Overlapping alignments
#'   are counted separately for each coordinate region. If multiple
#'   regions have identical names, their counts will be summed, counting
#'   each alignment only once even if it overlaps more than one of these
#'   regions. Alignments may be counted more than once if they overlap
#'   multiple regions that have different names.
#'   This mode is for example used to quantify ChIP-seq alignments in
#'   promoter regions, or gene expression levels in an RNA-seq experiment
#'   (using a \code{query} with exon regions named by gene).
#'   \item \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}: Alignments are
#'   counted and summed for each list element in \code{query} if they
#'   overlap with any of the regions contained in the list element. The
#'   order of the list elements defines a hierarchy for quantification:
#'   Alignment will only be counted for the first element (the one with
#'   the lowest index in \code{query}) that they overlap, but not for any
#'   potential further list elements containing overlapping regions.
#'   This mode can be used to hierarchically and uniquely count (assign)
#'   each alignment to a one of several groups of regions (the elements
#'   in \code{query}), for example to estimate the fractions of different
#'   classes of RNA in an RNA-seq experiment (rRNA, tRNA, snRNA, snoRNA,
#'   mRNA, etc.)
#'   \item \code{\link[GenomicFeatures:TxDb-class]{TxDb}}: Used to extract
#'   regions from annotation and report alignment counts depending on the
#'   value of \code{reportLevel}. If \code{reportLevel="exon"},
#'   alignments overlapping each exon in \code{query} are counted.
#'   If \code{reportLevel="gene"}, alignment counts for all exons of a
#'   gene will be summed, counting each alignment only once even if it
#'   overlaps multiple annotated exons of a gene. These are useful to
#'   calculate exon or gene expression levels in RNA-seq experiments
#'   based on the annotation in a \code{TxDb} object. If
#'   \code{reportLevel="promoter"}, the \code{promoters} function from package
#'   \pkg{GenomicFeatures} is used with default arguments to extract
#'   promoter regions around transcript start sites, e.g. to quantify
#'   alignments inf a ChIP-seq experiment.
#'   \item any of the above or \code{NULL} for
#'   \code{reportLevel="junction"}: The \code{query} argument is ignored
#'   if \code{reportLevel} is set to \code{"junction"}, and \code{qCount}
#'   will count the number of alignments supporting each exon-exon
#'   junction detected in any of the samples in \code{proj}. The
#'   arguments \code{selectReadPosition}, \code{shift},
#'   \code{orientation}, \code{useRead} and \code{mask} will have no
#'   effect in this quantification mode.
#' }
#' 
#' The additional arguments allow to fine-tune the quantification:
#' 
#' \code{selectReadPosition} defines the part of the alignment that has
#' to be contained within a query region for an overlap. The values
#' \code{start} (default) and \code{end} refer to the biological start
#' (5'-end) and end (3'-end) of the alignment. For example, the
#' \code{start} of an alignment on the plus strand is its leftmost
#' (lowest) base, and the \code{end} of an alignment on the minus strand
#' is also the leftmost base.
#' 
#' \code{shift} allows on-the-fly shifting of alignments towards their
#' 3'-end prior to overlap determination and counting. This can be
#' helpful to increase resolution of ChIP-seq experiments by moving
#' alignments by half the immuno-precipitated fragment size towards the
#' middle of fragments. \code{shift} is either an \dQuote{integer} vector
#' with one value per alignment file in \code{proj}, or a single
#' \dQuote{integer} value, in which case all alignment files will be
#' shifted by the same value. For paired-end experiments, it can be
#' alternatively set to "halfInsert", which will estimate the true
#' fragment size from the distance between aligned read pairs and shift
#' the alignments accordingly.
#' 
#' \code{orientation} controls the interpretation of alignment strand
#' when counting, relative to the strand of the query region. \code{any}
#' will count all overlapping alignments, irrespective of the alignment
#' strand (e.g. used in an unstranded RNA-seq experiment). \code{same}
#' will only count the alignments on the same strand as the query region
#' (e.g. in a stranded RNA-seq experiment), and \code{opposite} will only
#' count the alignments on the opposite strand from the query region
#' (e.g. to quantify anti-sense transcription in a stranded RNA-seq
#' experiment).
#' 
#' \code{includeSpliced} and \code{includeSecondary} can be used to
#' include or exclude spliced or secondary alignments,
#' respectively. \code{mapqMin} and \code{mapqMax} allow to select alignments
#' based on their mapping qualities. \code{mapqMin} and \code{mapqMax} can
#' take integer values between 0 and 255 and equal to
#' \eqn{-10 log_{10} Pr(\textnormal{mapping position is wrong})}{-10
#' log10 Pr(mapping position is wrong)}, rounded to the nearest
#' integer. A value 255 indicates that the mapping quality is not available.
#' 
#' In paired-end experiments, \code{useRead} allows to quantify either
#' all alignments (\code{useRead="any"}), or only the first
#' (\code{useRead="first"}) or last (\code{useRead="last"}) read from a
#' read pair or read group. Note that for \code{useRead="any"} (the
#' default), an alignment pair that is fully contained within a query
#' region will contribute two counts to the value of that
#' region. \code{absIsizeMin} and \code{absIsizeMax} can be used to
#' select alignments based on their insert size (TLEN field in SAM Spec
#' v1.4).
#' 
#' \code{auxiliaryName} selects the reference sequence for which
#' alignments should be quantified. \code{NULL} (the default) will
#' select alignments against the genome. If set to a character string
#' that matches one of the auxiliary target names (as specified in
#' the \code{auxiliaryFile} argument of \code{\link[QuasR]{qAlign}}), the
#' corresponding alignments will be counted.
#' 
#' \code{mask} can be used to specify a
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with regions in the
#' reference sequence to be excluded from quantification. The regions
#' will be considered unstranded (\code{strand="*"}). Alignments that
#' overlap with a region in \code{mask} will not be counted. Masking may
#' reduce the effective width of query regions reported by \code{qCount},
#' even down to zero for regions that are fully contained in \code{mask}.
#' 
#' If \code{clObj} is set to an object that inherits from class
#' \code{cluster}, for example an object returned by
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel}, the
#' quantification task is split into multiple chunks and processed in
#' parallel using \code{\link[parallel:clusterApply]{clusterMap}}. Currently, not all
#' tasks will be efficiently parallelized: For example, a single query
#' region and a single (group of) bam files will not be split into
#' multiple chunks.
#' 
#' @param proj A \code{\linkS4class{qProject}} object representing a
#'   sequencing experiment as returned by \code{\link[QuasR]{qAlign}}
#' @param query An object of type \code{\link[GenomicRanges:GRanges-class]{GRanges}}, 
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} or 
#'   \code{\link[GenomicFeatures:TxDb-class]{TxDb}} with the regions to be
#'   quantified. The type of \code{query} will determine the mode of
#'   quantification (see \sQuote{Details}). For
#'   \code{reportLevel="junction"}, \code{query} is ignored and can also
#'   be \code{NULL}.
#' @param reportLevel Level of quantification (\code{query} of type
#'   \code{TxDb} or \code{NULL}), one of
#'   \itemize{
#'     \item \code{gene} (default): one value per gene
#'     \item \code{exon}: one value per exon
#'     \item \code{promoter}: one value per promoter
#'     \item \code{junction}: one count per detected exon-exon junction
#"     (\code{query} will be ignored in this case)
#'   }
#' @param selectReadPosition defines the part of the alignment that has
#'   to be contained within a query region to produce an overlap (see
#'   \sQuote{Details}). Possible values are:
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
#'   auxiliary alignments (see \sQuote{Details}).
#' @param mask If not \code{NULL}, a \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#'   object with reference regions to be masked, i.e. excluded from the
#'   quantification, such as unmappable or highly repetitive regions (see
#'   \sQuote{Details}).
#' @param collapseBySample If \code{TRUE} (the default), sum alignment
#'   counts from bam files with the same sample name.
#' @param includeSpliced If \code{TRUE} (the default), include spliced
#'   alignments when counting. A spliced alignment is defined as an
#'   alignment with a gap in the read of at least 60 bases.
#' @param includeSecondary If \code{TRUE} (the default), include alignments
#'   with the secondary bit (0x0100) set in the \code{FLAG} when counting.
#' @param mapqMin Minimal mapping quality of alignments to be included when
#'   counting (mapping quality must be greater than or equal to
#'   \code{mapqMin}). Valid values are between 0 and 255. The default (0)
#'   will include all alignments.
#' @param mapqMax Maximal mapping quality of alignments to be included when
#'   counting (mapping quality must be less than or equal to \code{mapqMax}).
#'   Valid values are between 0 and 255. The default (255) will include
#'   all alignments.
#' @param absIsizeMin For paired-end experiments, minimal absolute insert
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included when
#'   counting. Valid values are greater than 0 or \code{NULL} (default),
#'   which will not apply any minimum insert size filtering.
#' @param absIsizeMax For paired-end experiments, maximal absolute insert
#'   size (TLEN field in SAM Spec v1.4) of alignments to be included when
#'   counting. Valid values are greater than 0 or \code{NULL} (default),
#'   which will not apply any maximum insert size filtering.
#' @param maxInsertSize Maximal fragment size of the paired-end experiment. 
#'   This parameter is used if \code{shift="halfInsert"} and will
#'   ensure that query regions are made wide enough to emcompass all
#'   alignment pairs whose mid falls into the query region. The default
#'   value is \code{500} bases.
#' @param clObj A cluster object to be used for parallel processing (see
#'   \sQuote{Details}).
#' 
#' @name qCount
#' @aliases qCount
#' 
#' @return 
#' A \code{matrix} with effective query regions width in the first
#' column, and alignment counts in subsequent columns, or a
#' \code{GRanges} object if \code{reportLevel="junction"}.
#' 
#' The effective query region width returned as first column in the
#' matrix is calculated by the number of unique, non-masked bases in the
#' reference sequence that contributed to the count of this query
#' name (irrespective if the bases were covered by alignments or not).
#' An effective width of zero indicates that the region was fully
#' masked and will have zero counts in all samples.
#' 
#' The alignment counts in the matrix are contained from column two
#' onwards. For projects with allele-specific quantification, i.e. if a
#' file with single nucleotide polymorphisms was supplied to the
#' \code{snpFile} argument of \code{\link[QuasR]{qAlign}}, there will be
#' three columns per bam file (number of alignments for Reference, 
#' Unknown and Alternative genotypes, with suffixed _R, _U and
#' _A). Otherwise there is a single columns per bam file.
#' 
#' If \code{collapseBySample}=\code{TRUE}, groups of bam files with identical 
#' sample name are combined by summing their alignment counts.
#' 
#' For \code{reportLevel="junction"}, the return value is a
#' \code{GRanges} object. The start and end coordinates correspond to the
#' first and last base in each detected intron. Plus- and minus-strand
#' alignments are quantified separately, so that in an unstranded RNA-seq
#' experiment, the same intron may be represented twice; once for each
#' strand. The counts for each sample are contained in the \code{mcols}
#' of the \code{GRanges} object.
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' @keywords utilities misc
#' 
#' @export
#' 
#' @seealso 
#' \code{\link[QuasR]{qAlign}},
#' \code{\linkS4class{qProject}},
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
#' 
#' @examples 
#' library(GenomicRanges)
#' library(Biostrings)
#' library(Rsamtools)
#' 
#' # copy example data to current working directory
#' file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
#' 
#' # load genome sequence
#' genomeFile <- "extdata/hg19sub.fa"
#' gseq <- readDNAStringSet(genomeFile)
#' chrRegions <- GRanges(names(gseq), IRanges(start=1,width=width(gseq),names=names(gseq)))
#' 
#' # create alignments (paired-end experiment)
#' sampleFile <- "extdata/samples_rna_paired.txt"
#' proj <- qAlign(sampleFile, genomeFile, splicedAlignment=TRUE)
#' 
#' # count reads using a "GRanges" query
#' qCount(proj, query=chrRegions)
#' qCount(proj, query=chrRegions, useRead="first")
#' 
#' # hierarchical counting using a "GRangesList" query
#' library(rtracklayer)
#' annotationFile <- "extdata/hg19sub_annotation.gtf"
#' gtfRegions <- import.gff(annotationFile, format="gtf", feature.type="exon")
#' names(gtfRegions) <- mcols(gtfRegions)$source
#' gtfRegionList <- split(gtfRegions, names(gtfRegions))
#' names(gtfRegionList)
#' 
#' res3 <- qCount(proj, gtfRegionList)
#' res3
#' 
#' # gene expression levels using a "TxDb" query
#' library("GenomicFeatures")
#' genomeRegion <- scanFaIndex(genomeFile)
#' chrominfo <- data.frame(chrom=as.character(seqnames(genomeRegion)),
#'                         length=end(genomeRegion),
#'                         is_circular=rep(FALSE, length(genomeRegion)))
#' txdb <- makeTxDbFromGFF(annotationFile, 
#'                         format="gtf", 
#'                         chrominfo=chrominfo,
#'                         dataSource="Ensembl modified",
#'                         organism="Homo sapiens")
#'
#' res4 <- qCount(proj, txdb, reportLevel="gene")
#' res4
#' 
#' # exon-exon junctions
#' res5 <- qCount(proj, NULL, reportLevel="junction")
#' res5
#' 
qCount <- function(proj,
                   query,
                   reportLevel = c(NULL, "gene", "exon", "promoter", "junction"),
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
                   clObj = NULL) {
    ## setup variables from 'proj' ---------------------------------------------
    ## 'proj' is correct type?
    if (!inherits(proj, "qProject", which = FALSE))
        stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")
    
    samples <- proj@alignments$SampleName
    nsamples <- length(samples)
    bamfiles <-
        if(is.null(auxiliaryName))
            proj@alignments$FileName
    else if(!is.na(i <- match(auxiliaryName, proj@aux$AuxName)))
        unlist(proj@auxAlignments[i, ], use.names = FALSE)
    else
        stop("unknown 'auxiliaryName', should be one of: NULL, ",
             paste(sprintf("'%s'", proj@aux$AuxName), collapse = ", "))
    
    
    ## validate parameters -----------------------------------------------------
    reportLevel <- match.arg(reportLevel)
    selectReadPosition <- match.arg(selectReadPosition)
    orientation <- match.arg(orientation)
    useRead <- match.arg(useRead)
    if ((!is.null(absIsizeMin) || !is.null(absIsizeMax)) && proj@paired == "no")
        stop("'absIsizeMin' and 'absIsizeMax' can only be used for paired-end experiments")
    if (is.null(absIsizeMin)) # -1L -> do not apply TLEN filtering
        absIsizeMin <- -1L
    if (is.null(absIsizeMax))
        absIsizeMax <- -1L
    
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
        else if(length(shift) == 1)
            shifts <- rep(as.integer(shift), nsamples)
        else
            shifts <- as.integer(shift)
        broaden <- 0L
    }
    
    ## all query chromosomes present in all bamfiles?
    trTab <- table(unlist(lapply(scanBamHeader(bamfiles), function(bh) names(bh$targets))))
    trCommon <- names(trTab)[trTab == length(bamfiles)]
    queryseqs <- NULL
    if (inherits(query, c("GRanges","GRangesList"))) {
        queryseqs <- seqlevelsInUse(query)
    } else if (inherits(query,"TxDb")) { # only use active sequences
        queryseqs <- seqlevels(query) # isActiveSeq is depricated; no seqlevelsInUse for TxDb yet
    }
    if (!is.null(query) && any(f <- !(queryseqs %in% trCommon)))
        stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                     paste(seqlevels(query)[f], collapse = ", ")))
    
    ## 'query' is correct type?
    if (reportLevel == "junction") {
        if (proj@splicedAlignment != TRUE && proj@samplesFormat != "bam")
            stop("reportLevel=\"junction\" cannot be used for non-spliced alignments")
        if (absIsizeMin != -1L || absIsizeMax != -1L) {
            warning("ignoring 'absIsizeMin' and 'absIsizeMax' for reportLevel=\"junction\"")
            absIsizeMin <- absIsizeMax <- -1L
        }
        if (!is.null(query))
            warning("ignoring 'query' for reportLevel=\"junction\"")
        
        ### reportLevel == "junction" ------------------------------------------
        ## setup tasks for parallelization -------------------------------------
        if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
            message("preparing to run on ", length(clObj), " nodes...", appendLF = FALSE)
            ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
            if (!all(sapply(ret, function(x) "QuasR" %in% x)))
                stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
            taskTargets <- rep(trCommon, nsamples)
            bamfiles <- rep(bamfiles, each = length(trCommon))
            iByBamfile <- split(seq_along(bamfiles), bamfiles)[unique(bamfiles)]
            myapply <- function(...) {
                ret <- clusterMap(clObj, ..., SIMPLIFY = FALSE, .scheduling = "dynamic")
                # fuse
                if(!is.na(proj@snpFile)){ # ret is a list of list(id,R,U,A)
                    ret <- lapply(iByBamfile, function(i) 
                        list(id = do.call(c, lapply(ret[i], "[[", "id")),
                             R = do.call(c, lapply(ret[i], "[[", "R")),
                             U = do.call(c, lapply(ret[i], "[[", "U")),
                             A = do.call(c, lapply(ret[i], "[[", "A"))))
                } else {                  # ret is a list of named vectors
                    ret <- lapply(iByBamfile, function(i) do.call(c, unname(ret[i])))
                }
                ret
            }
            message("done")
        } else {
            taskTargets <- rep(list(NULL), length(bamfiles))
            myapply <- function(...) {
                ret <- mapply(..., SIMPLIFY = FALSE)
                ret
            }
        }
        
        ## count junctions -----------------------------------------------------
        message("counting junctions...", appendLF = FALSE)
        resL <- myapply(countJunctionsOneBamfile,
                        bamfile = bamfiles,
                        targets = taskTargets,
                        MoreArgs = list(allelic = !is.na(proj@snpFile),
                                        includeSecondary = includeSecondary,
                                        mapqmin = as.integer(mapqMin)[1],
                                        mapqmax = as.integer(mapqMax)[1]))
        message("done")
        
        ## make result rectangular and collapse (sum) counts by sample if necessary
        if (!is.na(proj@snpFile)) { # ret is a list of list(id,R,U,A)
            allJunctions <- unique(Reduce(c, lapply(resL, "[[", "id")))
            res <- matrix(0, nrow = length(allJunctions), 
                          ncol = 3*length(resL), dimnames = list(allJunctions, NULL))
            for (i in 1:length(resL)) # make res a matrix with 3 columns per sample
                res[resL[[i]]$id, ((i-1)*3+1):((i-1)*3+3)] <- 
                do.call(cbind, resL[[i]][c("R", "U", "A")])
            if (nsamples > length(unique(samples))) {
                if (collapseBySample) {
                    message("collapsing counts by sample...", appendLF = FALSE)
                    iBySample <- split(seq_len(nsamples), samples)[unique(samples)]
                    res <- do.call(cbind, lapply(iBySample, function(i)
                        cbind(R = rowSums(res[, (i-1)*3+1, drop = FALSE]),
                              U = rowSums(res[, (i-1)*3+2, drop = FALSE]),
                              A = rowSums(res[, (i-1)*3+3, drop = FALSE]))))
                    colnames(res) <- paste(rep(names(iBySample), each = 3), 
                                           c("R", "U", "A"), sep = "_")
                    message("done")
                    
                } else {
                    # unify non-collapsed identical sample names
                    colnames(res) <- paste(rep(displayNames(proj), each = 3),
                                           c("R", "U", "A"), sep = "_")
                }
            } else {
                colnames(res) <- paste(rep(samples, each = 3), c("R", "U", "A"), 
                                       sep = "_")
            }
            
        } else {                  # ret is a list of named vectors
            allJunctions <- unique(Reduce(c, lapply(resL, names)))
            res <- matrix(0, nrow = length(allJunctions), ncol = length(resL),
                          dimnames = list(allJunctions, NULL))
            for (i in 1:length(resL))
                res[names(resL[[i]]),i] <- resL[[i]]
            if (nsamples > length(unique(samples))) {
                if (collapseBySample) {
                    message("collapsing counts by sample...", appendLF = FALSE)
                    iBySample <- split(seq_len(nsamples), samples)[unique(samples)]
                    res <- do.call(cbind, lapply(iBySample, function(i) 
                        rowSums(res[, i, drop = FALSE])))
                    message("done")
                } else {
                    # unify non-collapsed identical sample names
                    colnames(res) <- displayNames(proj)
                }
            } else {
                colnames(res) <- samples
            }
        }

        ## make GRanges object
        res2 <- GRanges(seqnames = sub("^(.+):[0-9]+:[0-9]+:.$", "\\1", allJunctions),
                        ranges = IRanges(
                            start = as.integer(sub("^.+:([0-9]+):[0-9]+:.$", 
                                                   "\\1", allJunctions)),
                            end = as.integer(sub("^.+:[0-9]+:([0-9]+):.$",
                                                 "\\1", allJunctions))),
                        strand = sub("^.+:[0-9]+:[0-9]+:(.)$", "\\1", allJunctions))
        mcols(res2) <- res
        
        ## return results
        return(res2)
    } else {
        ### reportLevel != "junction" ------------------------------------------
        if (!inherits(query, c("GRanges", "GRangesList", "TxDb")))
            stop("'query' must be either an object of type 'GRanges', 'GRangesList' or 'TxDb', or NULL for reportLevel=\"junction\"")
        
        ## 'useRead' set but not a paired-end experiment?
        if (useRead != "any" && proj@paired == "no")
            warning("ignoring 'useRead' for single read experiments")
        
        
        ## preprocess query ----------------------------------------------------
        ##    --> create 'flatquery', 'querynames', 'querylengths' and 'zeroquerynames'
        ##    GRanges query ----------------------------------------------------
        if (inherits(query,"GRanges")) {
            if (!is.null(names(query)) && length(query) > length(unique(names(query)))) {
                # remove redundancy from 'query' by names
                tmpquery <- reduce(split(query, names(query))[unique(names(query))])
                flatquery <- unlist(tmpquery, use.names = FALSE)
                querynames <- rep(names(tmpquery), elementNROWS(tmpquery))
                rm(tmpquery)
            } else {
                flatquery <- query
                querynames <- if (is.null(names(query))) 
                    as.character(seq_len(length(query))) else names(query)
            }
            querylengths <- width(flatquery)
            zeroquerynames <- character(0)
            
            ##    GRangesList query --------------------------------------------
        } else if (inherits(query, "GRangesList")) {
            if (any(i <- elementNROWS(query) == 0)) {
                warning(sprintf("removing %d elements from 'query' with zero regions: %s",
                                sum(i), paste(names(query)[i], collapse = ", ")))
                query <- query[-i]
            }
            # hierarchically remove redundancy from 'query'
            message("hierarchically removing redundancy from 'query'...", 
                    appendLF = FALSE)
            if (orientation == "any")
                strand(query) <- endoapply(strand(query), function(x) 
                    Rle(factor("*", levels = c("+", "-", "*")), lengths = length(x)))
            query <- reduce(query)
            if (length(query) > 1) {
                cumquery <- query[[1]]
                for (i in 2:length(query)) {
                    query[[i]] <- setdiff(query[[i]], cumquery)
                    cumquery <- c(query[[i]], cumquery)
                }
            }
            message("done")
            flatquery <- unlist(query, use.names = FALSE)
            querynames <- rep(if(is.null(names(query))) as.character(seq_len(length(query))) else names(query),
                              elementNROWS(query))
            querylengths <- unlist(width(query), use.names = FALSE)
            zeroquerynames <- (if(is.null(names(query))) as.character(seq_len(length(query))) else names(query))[elementNROWS(query) == 0]
            
            ##    TxDb query -----------------------------------------------
        } else if (inherits(query, "TxDb")) {
            if (is.null(reportLevel))
                stop("'reportLevel' must be set to a non-NULL value for 'query' of type 'TxDb'")
            message(sprintf("extracting %s regions from TxDb...", reportLevel), 
                    appendLF = FALSE)
            if (reportLevel == "gene") {
                tmpquery <- reduce(exonsBy(query, by = "gene"))
                flatquery <- unlist(tmpquery, use.names = FALSE)
                querynames <- rep(names(tmpquery), elementNROWS(tmpquery))
                querylengths <- unlist(width(tmpquery), use.names = FALSE)
                rm(tmpquery)
                
            } else if (reportLevel == "exon") {
                flatquery <- exons(query, columns = "exon_id")
                querynames <- as.character(mcols(flatquery)$exon_id)
                querylengths <- width(flatquery)
                
            } else if (reportLevel == "promoter") {
                flatquery <- promoters(query, columns = c("tx_id", "tx_name"))
                querynames <- paste(as.character(mcols(flatquery)$tx_id),
                                    as.character(mcols(flatquery)$tx_name), sep = ";")
                querylengths <- width(flatquery)
            }
            zeroquerynames <- character(0)
            message("done")
        }
        if (length(flatquery) == 0)
            stop("'query' is empty - nothing to do")
        ## from now on, only use 'flatquery' (GRanges object) with names in 'querynames' and lengthes in 'querylengths'
        
        
        ## apply 'mask' to flatquery -------------------------------------------
        if (!is.null(mask)) {
            if (!inherits(mask,"GRanges"))
                stop("'mask' must be an object of type 'GRanges'")
            message("removing 'mask' ranges from 'query'...", appendLF = FALSE)
            strand(mask) <- "*"
            mask <- reduce(mask)
            ov <- findOverlaps(flatquery, mask)
            qOM <- unique(queryHits(ov))
            gr1 <- GRanges(seqnames = Rle(qOM), 
                           ranges = IRanges(start = start(flatquery)[qOM], 
                                            end = end(flatquery)[qOM]))
            gr2 <- GRanges(seqnames = Rle(queryHits(ov)), 
                           ranges = IRanges(start = start(mask)[subjectHits(ov)], 
                                            end = end(mask)[subjectHits(ov)]))
            SD <- setdiff(gr1, gr2)
            notOverlappingMaskInd <- which(!((1:length(flatquery)) %in% qOM))
            completelyMaskedInd <- 
                qOM[!(qOM %in% as.numeric(as.character(unique(seqnames(SD)))))]
            
            # 'zeroquery' contains regions that are completely masked (will get zero count and length)
            #zeroquery <- flatquery[completelyMaskedInd]
            tmpnames <- querynames[completelyMaskedInd]
            zeroquerynames <- c(zeroquerynames, tmpnames[!(tmpnames %in% querynames[-completelyMaskedInd])])
            #zeroquerylengths <- rep(0L, length(completelyMaskedInd))
            
            # masked 'flatquery' maybe split into several non-masked pieces
            flatquery <- c(flatquery[notOverlappingMaskInd],
                           GRanges(seqnames = seqnames(flatquery)[as.numeric(as.character(seqnames(SD)))],
                                   ranges = ranges(SD), 
                                   strand = strand(flatquery)[as.numeric(as.character(seqnames(SD)))],
                                   seqlengths = seqlengths(flatquery)),
                           ignore.mcols = TRUE)
            querynames <- querynames[c(notOverlappingMaskInd, as.numeric(as.character(seqnames(SD))))]
            querylengths <- width(flatquery)
            message("done")
        }
        
        ## setup tasks for parallelization -------------------------------------
        ## TODO: if sum(width(flatquery)) close to sum(seqlengths(genome)) -> select variant counting algorithm (sequential walk through bamfiles)
        if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
            message("preparing to run on ", length(clObj), " nodes...", appendLF = FALSE)
            ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
            if (!all(sapply(ret, function(x) "QuasR" %in% x)))
                stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
            taskIByFlatQuery <- splitIndices(nx = length(flatquery),
                                             ncl = ceiling(length(clObj)/nsamples*2))
            if (inherits(taskIByFlatQuery, "integer", which = FALSE))
                taskIByFlatQuery <- list(taskIByFlatQuery) # make sure taskIByFlatQuery is a list, even if ceiling(length(clObj) /nsamples *2)==1
            taskSamples <- rep(samples, each = length(taskIByFlatQuery))
            taskBamfiles <- rep(bamfiles, each = length(taskIByFlatQuery))
            flatquery <- lapply(taskIByFlatQuery, function(i) flatquery[i])
            shifts <- rep(shifts, each = length(taskIByFlatQuery))
            myapply <- function(...) {
                ret <- clusterMap(clObj, ..., SIMPLIFY = FALSE, .scheduling = "dynamic")
                ## fuse
                iBySample <- split(seq_along(ret),names(ret))[unique(names(ret))]
                names(ret) <- NULL
                if (!is.na(proj@snpFile)) {
                    ret <- do.call(cbind, lapply(iBySample, function(i) 
                        do.call(rbind, ret[i])))
                    postfix <- substring(colnames(ret), nchar(colnames(ret)))
                    ## rename
                    dimnames(ret) <- list(querynames, 
                                          paste(rep(samples, each = 3), 
                                                postfix, sep = "_"))
                } else {
                    ret <- do.call(cbind, lapply(iBySample, function(i) 
                        do.call(c, ret[i])))
                    ## rename
                    dimnames(ret) <- list(querynames, samples)
                }
                ret
            }
            message("done")
        } else {
            taskSamples <- samples
            taskBamfiles <- bamfiles
            flatquery <- list(flatquery)
            myapply <- function(...) {
                ret <- do.call(cbind, mapply(..., SIMPLIFY = FALSE))
                ## rename
                if (!is.na(proj@snpFile))
                    dimnames(ret) <- list(querynames, 
                                          paste(rep(samples, each = 3),
                                                substring(colnames(ret), 
                                                          nchar(colnames(ret))), 
                                                sep = "_"))
                else
                    dimnames(ret) <- list(querynames, samples)
                ret
            }
        }
        
        ## count alignments ----------------------------------------------------
        message("counting alignments...", appendLF = FALSE)
        res <- myapply(countAlignments,
                       bamfile = taskBamfiles,
                       regions = flatquery,
                       shift = shifts,
                       MoreArgs = list(
                           selectReadPosition = selectReadPosition,
                           orientation = orientation,
                           useRead = useRead,
                           broaden = broaden,
                           allelic = !is.na(proj@snpFile),
                           includeSpliced = includeSpliced,
                           includeSecondary = includeSecondary,
                           mapqmin = as.integer(mapqMin)[1],
                           mapqmax = as.integer(mapqMax)[1],
                           absisizemin = as.integer(absIsizeMin)[1],
                           absisizemax = as.integer(absIsizeMax)[1]))
        message("done")
        
        
        ## collapse (sum) counts by sample if necessary
        if (nsamples > length(unique(samples))) {
            if (collapseBySample) {
                message("collapsing counts by sample...", appendLF = FALSE)
                if (is.na(proj@snpFile))
                    iBySample <- split(seq_len(nsamples), samples)[unique(samples)]
                else
                    iBySample <- split(seq_len(ncol(res)), 
                                       colnames(res))[unique(colnames(res))]
                res <- do.call(cbind, lapply(iBySample, function(i) 
                    rowSums(res[, i, drop = FALSE])))
                message("done")
                
            } else {
                # unify non-collapsed identical sample names
                if (is.na(proj@snpFile))
                    colnames(res) <- displayNames(proj)
                else
                    colnames(res) <- paste(rep(displayNames(proj), each = 3),
                                           substring(colnames(res), 
                                                     nchar(colnames(res))), sep = "_")
            }
        }
        
        ## add the region width as first column
        res <- cbind(width = querylengths, res)
        rm(querylengths)
        
        ## collapse (sum) counts by 'querynames'
        if (length(querynames) > length(unique(querynames))) {
            message("collapsing counts by query name...", appendLF = FALSE)
            iByQuery <- split(seq_len(nrow(res)), querynames)[unique(querynames)]
            res <- do.call(rbind, lapply(iByQuery, function(i) 
                colSums(res[i, 1:ncol(res), drop = FALSE])))
            rownames(res) <- querynames <- names(iByQuery)
            message("done")
        }
        if (length(zeroquerynames) > length(unique(zeroquerynames)))
            zeroquerynames <- unique(zeroquerynames)
        
        ## combine with zeroquery and reorder according to 'query'
        res2 <- matrix(0, nrow = length(querynames) + length(zeroquerynames), 
                       ncol = ncol(res),
                       dimnames = list(if(inherits(query, "TxDb"))
                           sort(c(querynames, zeroquerynames))
                           else if(is.null(names(query)))
                               as.character(seq_len(length(query)))
                           else unique(names(query)),
                           colnames(res)))
        res2[rownames(res), ] <- res
        
        ## return results
        return(res2)
    }
}


## count junctions (with the C-function) for single bamfile and optionally selected target sequences
## return a named vector with junction elements (names of the form "chromosome:first_intronic_base:last_intronic_base:strand")
#' @keywords internal
countJunctionsOneBamfile <- function(bamfile, targets, allelic, 
                                     includeSecondary, mapqmin, mapqmax) {
    tryCatch({ # try catch block goes through the whole function
        # prepare region vectors
        bh <- scanBamHeader(bamfile)[[1]]$targets
        tid <- seq_along(bh) - 1L
        if (!is.null(targets)) {
            if (any(is.na(i <- match(targets,names(bh)))))
                stop(sprintf("some targets not found in bamfile '%s': %s",
                             bamfile, targets[which(is.na(i))]))
            bh <- bh[i]
            tid <- tid[i]
        }
        start <- rep(0L, length(tid)) # samtools library has 0-based inclusive start
        end <- unname(bh) ## samtool library has 0-based exclusiv end
        # count junctions
        count <- .Call(countJunctions, bamfile, tid, start, end, allelic,
                       includeSecondary, mapqmin, mapqmax)
        return(count)

    }, error = function(ex) {
        emsg <- paste("Internal error on", Sys.info()['nodename'], 
                      "query bamfile", bamfile, 
                      "\n Error message is:", ex$message)
        stop(emsg)
    })
}


## count alignments (with the C-function) for single bamfile, single shift, and single set of regions
## return a numeric vector with length(regions) elements (same order as regions)
#' @keywords internal
countAlignments <- function(bamfile, regions, shift, selectReadPosition, orientation,
                            useRead, broaden, allelic, includeSpliced, includeSecondary,
                            mapqmin, mapqmax, absisizemin, absisizemax) {
    tryCatch({ # try catch block goes through the whole function
        
        ## translate seqnames to tid and create region data.frame
        seqnamesBamHeader <- names(scanBamHeader(bamfile)[[1]]$targets)
        
        ## prepare region vectors
        #tid <- IRanges::as.vector(IRanges::match(seqnames(regions), seqnamesBamHeader)) - 1L
        tid <- as.vector(match(seqnames(regions), seqnamesBamHeader) - 1L) 
        start <- start(regions) - 1L ## samtool library has 0-based inclusiv start
        end <- end(regions) ## samtools library has 0-based exclusive end
        
        ## swap strand for 'orientation="opposite"' 
        if (orientation == "any")
            strand <- rep("*", length(regions))
        else if (orientation == "opposite")
            strand <- c("+"="-", "-"="+", "*"="*")[as.character(strand(regions))]
        else # orientation == "same"
            strand <- as.character(strand(regions))
        
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
        
        ## get counts
        if (!allelic) {
            count <- .Call(countAlignmentsNonAllelic, bamfile, tid, start, end, strand,
                           selectReadPosition, readBitMask, shift, broaden, includeSpliced,
                           mapqmin, mapqmax, absisizemin, absisizemax)
        } else {
            count <- as.matrix(as.data.frame(
                .Call(countAlignmentsAllelic, bamfile, tid, start, end, strand,
                      selectReadPosition, readBitMask, shift, broaden, includeSpliced,
                      mapqmin, mapqmax, absisizemin, absisizemax)
            ))
        }
        
        return(count)
    }, error = function(ex) {
        reg <- regions[c(1, length(regions))]
        emsg <- paste("Internal error on", Sys.info()['nodename'], 
                      "query bamfile", bamfile,"with regions\n", 
                      paste(seqnames(reg), start(reg), "-" , end(reg), 
                            strand(reg), collapse = "\n\t...\n"), 
                      "\n Error message is:", ex$message)
        stop(emsg)
    })
}



# ## Counts alignments in a given set of regions which are located in a subspace of the genome
# ## using a per-base-coverage vector approach
# ##     shift the read, broaden fetch region,
# ## return a numeric vector with length(regions) elements (same order as regions)
# countAlignmentsSubregionsC <- function(bamfile, regions, selectReadPosition, shift=0L, broaden=0L, includeSpliced=TRUE)
# {
#     ## check if region are located on one chromosme
#     seqName <- unique(seqnames(regions))
#     if(length(seqName) > 1L)
#         stop("regions should only be located on one chromosome")
#     
#     ## check broaden and shift parameter
#     if(broaden < 0)
#         stop("'broaden' should not be negative") 
#     #    if(shift > 0 && selectReadPosition="midwithin")
#     #        stop("'shift' parameter must be zero if 'selectReadPosition' is set to midwithin")
#     #    if(broaden > 0 && (selectReadPosition="startwithin" || selectReadPosition="endwithin"))
#     #        stop("'broaden' parameter must be zero if 'selectReadPosition' is set to startwithin or endwithin")
#     
#     ## translate seqName to tid
#     seqnamesList <- names(scanBamHeader(bamfile)[[1]]$targets)
#     tidList <- as.integer(seq_along(seqnamesList)-1)
#     tid <- tidList[ match(seqName, seqnamesList) ]
#     
#     ## convert grange to data.frame 
#     ## with 0-based start inclusive
#     ## with 0-based end exclusive
#     regionsTable <- data.frame(start=as.integer(start(regions)-1), ## samtool library has 0-based start
#                                end=as.integer(end(regions)),
#                                strand=as.character(strand(regions)),
#                                stringsAsFactors=FALSE
#     )
#     
#     ## call c-function
#     cnt <- .Call(countAlignmentsSubregions, bamfile, bamfile, tid, min(regionsTable$start), max(regionsTable$end),
#                  regionsTable, as.integer(shift), as.integer(broaden), selectReadPosition, includeSpliced)
#     
#     return(cnt)
# }
