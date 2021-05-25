#' Update qProject
#' 
#' Check if a qProject is compatible with the current class definition, and 
#' update it if necessary.
#' 
#' @param proj a \code{\link[=qProject-class]{qProject}} object.
#' @param quiet logical, whether to suppress update messages.
#' 
#' @return A \code{\link[=qProject-class]{qProject}} object.
#' 
#' @author Charlotte Soneson
#' 
#' @seealso \code{\link[=qProject-class]{qProject}}
#' 
#' @keywords utilities misc
#' 
#' @name qProjectUpdate
#' @aliases qProjectUpdate
#' 
#' @export
#' 
#' @importFrom methods .hasSlot
#' 
qProjectUpdate <- function(proj, quiet = TRUE) {
    mod <- FALSE
    if (!(methods::.hasSlot(proj, "geneAnnotation"))) {
        proj@geneAnnotation <- NA_character_
        if (!quiet) message("Adding geneAnnotation slot")
        mod <- TRUE
    } 
    if (!(methods::.hasSlot(proj, "geneAnnotationFormat"))) {
        proj@geneAnnotationFormat <- NA_character_
        if (!quiet) message("Adding geneAnnotationFormat slot")
        mod <- TRUE
    }
    if (!mod) {
        if (!quiet) message("qProject is up to date")
    }
    proj
}

### qProject class definition
#' qProject objects
#' 
#' The qProject class is a container for the meta-data (e.g. sample
#' names, paths and names of sequence and alignment files) associated
#' with a high-throughput sequencing experiment analyzed with \code{QuasR}.
#' 
#' The qProject class is returned by \link{qAlign} and stores all 
#' information on a high-throughput sequencing experiment analyzed with
#' \code{QuasR}. qProject objects can be conveniently passed to
#' \sQuote{q}-functions (function name starting with the letter
#' \sQuote{q}). The information is stored in the following slots:
#' \describe{
#'   \item{\code{reads}}{a 'data.frame' with sequence read files.}
#'   \item{\code{reads_md5subsum}}{a 'data.frame' with fingerprints for
#'   sequence read files.}
#'   \item{\code{alignments}}{a 'data.frame' with alignment files.}
#'   \item{\code{samplesFormat}}{a 'character(1)' specifying the format
#'   of input files.}
#'   \item{\code{genome}}{a 'character(1)' specifying the reference genome.}
#'   \item{\code{genomeFormat}}{a 'character(1)' specifying the format of
#'   the reference genome.}
#'   \item{\code{aux}}{a 'data.frame' with auxiliary reference sequence files.}
#'   \item{\code{auxAlignments}}{a 'data.frame' with alignment files for
#'   auxiliary reference sequence files.}
#'   \item{\code{aligner}}{a 'character(1)' specifying the aligner.}
#'   \item{\code{maxHits}}{a 'numeric(1)' specifying the maximum number
#'   of alignments per sequence.}
#'   \item{\code{paired}}{a 'character(1)' specifying the paired-type;
#'   one of "no", "fr", "rf", "ff".}
#'   \item{\code{splicedAlignment}}{a 'logical(1)'; \code{TRUE} when
#'   performing spliced-alignments.}
#'   \item{\code{snpFile}}{a 'character(1)' with a file name containing
#'   SNP information.}
#'   \item{\code{bisulfite}}{a 'character(1)' defining the bisulfite
#'   type; one of "no", "dir", "undir".}
#'   \item{\code{alignmentParameter}}{a 'character(1)' with aligner
#'   command line parameters.}
#'   \item{\code{projectName}}{a 'character(1)' with the project name.}
#'   \item{\code{alignmentsDir}}{a 'character(1)' with the directory to
#'   be used to store alignment files.}
#'   \item{\code{lib.loc}}{a 'character(1)' with the library directory to
#'   use for installing of alignment index packages.}
#'   \item{\code{cacheDir}}{a 'character(1)' with a directory to use for
#'   temporary files.}
#'   \item{\code{alnModeID}}{a 'character(1)' used internally to indicate
#'   the alignment mode.}
#' }
#' 
#' @section Accessors:
#' In the following code snippets, \code{x} is a qProject object.
#' \describe{
#'   \item{}{\code{length(x)}: Gets the number of input files.}
#'   \item{}{\code{genome(x)}: Gets the reference genome as a 'character(1)'. 
#'   The type of genome is stored as an attribute in
#'   \code{attr(genome(x),"genomeFormat")}: "BSgenome" indicates that
#'   \code{genome(x)} refers to the name of a BSgenome package, "file"
#'   indicates that it contains the path and filename of a genome in
#'   FASTA format.}
#'   \item{}{\code{auxiliaries(x)}: Gets a \code{data.frame} with auxiliary 
#'   target sequences, with one row per auxiliary target, and columns 
#'   "FileName" and "AuxName".}
#'   \item{}{ \code{alignments(x)}: Gets a list with two elements "genome" and
#'   "aux". \code{alignments(x)$genome} contains a \code{data.frame} 
#'   with \code{length(x)} rows and the columns "FileName" (containing
#'   the path to bam files with genomic alignments) and 
#'   "SampleName". \code{alignments(x)$aux} contains a 
#'   \code{data.frame} with one row per auxiliary target sequence (with
#'   auxiliary names as row names), and \code{length(x)} columns.}
#' }
#' 
#' @section Subsetting:
#' In the following code snippets, \code{x} is a qProject object.
#' \describe{
#'   \item{}{\code{x[i]}: Get \code{qProject} object instance with \code{i} 
#'   input files, where \code{i} can be an NA-free logical, numeric, or 
#'   character vector.}
#' }
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' 
#' @aliases show,qProject-method [,qProject,ANY,missing,missing-method 
#'   alignments,qProject-method alignments auxiliaries,qProject-method auxiliaries 
#'   genome,qProject-method length,qProject-method qProject qProject-class
#'   class:qProject
#'   
#' @name qProject-class
#' @docType class
#' 
#' @seealso \link{qAlign}
#' 
#' @export
#' 
#' @examples 
#' # copy example data to current working directory
#' file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)
#' 
#' # create alignments
#' sampleFile <- "extdata/samples_chip_single.txt"
#' genomeFile <- "extdata/hg19sub.fa"
#' auxFile <- "extdata/auxiliaries.txt"
#' 
#' proj <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile)
#' proj
#' 
#' # alignment statistics using a qProject
#' alignmentStats(proj)
#' 
#' # alignment statistics using bam files
#' alignmentStats(alignments(proj)$genome$FileName)
#' alignmentStats(unlist(alignments(proj)$aux))
#' 
setClass("qProject",
         slots = c(reads = "data.frame",
                   reads_md5subsum = "data.frame",
                   alignments = "data.frame",
                   samplesFormat = "character",
                   genome = "character",
                   genomeFormat = "character",
                   aux = "data.frame",
                   auxAlignments = "data.frame",
                   aligner = "character",
                   maxHits = "numeric",
                   paired = "character",
                   splicedAlignment = "logical",
                   snpFile = "character",
                   bisulfite = "character",
                   alignmentParameter = "character",
                   projectName = "character",
                   alignmentsDir = "character",
                   lib.loc = "character",
                   cacheDir = "character",
                   alnModeID = "character",
                   geneAnnotation = "character",
                   geneAnnotationFormat = "character")
)

### Methods
#' @export
setMethod("length", "qProject", function(x) nrow(x@reads))

#setGeneric("sampleNames", function(x) names(x))
#setMethod("sampleNames", "qProject", function(x) x@reads$SampleName)

#setGeneric("sampleNames<-", function(x, value) names(x) <- value)
#setReplaceMethod("sampleNames", "qProject", function(x, value) {
#    if(is.factor(value))
#        values <- as.character(value)
#    if(!is.character(value))
#        stop("sample names must be of type 'character' or 'factor'")
#    if(length(value) != nrow(x@reads))
#        stop(sprintf("length of sample names (%d) must be equal to the number of samples in qProject object (%d)",length(value),nrow(x@reads)))
#    x@reads$SampleName <- x@alignments$SampleName <- value
#    x
#})

#' @importFrom GenomeInfoDb genome
#' @export
setMethod("genome", signature(x = "qProject"), function(x) {
    y <- x@genome
    attr(y, "genomeFormat") <- x@genomeFormat
    return(y)
})

setGeneric("auxiliaries", function(x) return(NULL))
#' @export
setMethod("auxiliaries", signature(x = "qProject"), function(x) return(x@aux))

setGeneric("alignments", function(x) return(NULL))
#' @export
setMethod("alignments", signature(x = "qProject"), function(x) {
    return(list(genome = x@alignments, aux = x@auxAlignments))
})

#' @export
setMethod("[", signature(x = "qProject", i = "ANY", j = "missing",
                         drop = "missing"), 
          function(x, i) {
              y <- x
              if (is.character(i)) {
                  if (all(i %in% x@reads$SampleName)) {
                      if (any(inu <- vapply(unique(i), 
                                            function(ii) {
                                                sum(x@reads$SampleName == ii)
                                            }, 1) > 1))
                          warning("will only select first of multiple files ",
                                  "with identical sample name: ",
                                  paste(unique(i)[inu], collapse = ", "))
                      i <- match(i, x@reads$SampleName)
                  } else {
                      stop("undefined samples selected")
                  }
              } else if (is.logical(i) &&
                         length(i) != length(x@reads$SampleName)) {
                  stop("logical subsetting vector of length ", length(i), 
                       " does not match the number of sequence files (",
                       length(x@reads$SampleName), ")")
              } else if (is.numeric(i)) {
                  if (any(!is.finite(i)) ||
                      max(i) > length(x@reads$SampleName) ||
                      min(i) < 1)
                      stop("undefined samples selected")
              }
              y@reads <- y@reads[i, , drop = FALSE]
              y@reads_md5subsum <- y@reads_md5subsum[i, , drop = FALSE]
              y@alignments <- y@alignments[i, , drop = FALSE]
              if (ncol(y@auxAlignments))
                  y@auxAlignments <- y@auxAlignments[, i, drop = FALSE]
              return(y)
          })

#setGeneric("niceprint", function(x) print(x))
#setMethod("niceprint", "qProject", function(object) {
#' @importFrom methods show
#' @importFrom utils packageVersion
#' @export
setMethod("show", "qProject", function(object) {
    # project and global options
    cat("Project: " , object@projectName, "\n", sep = "")
    cat(" Options   : maxHits         : ", object@maxHits,
        "\n             paired          : ", object@paired,
        "\n             splicedAlignment: ", object@splicedAlignment,
        "\n             bisulfite       : ", object@bisulfite,
        "\n             snpFile         : ", if (length(object@snpFile) == 0 || is.na(object@snpFile)) "none" else truncPath(object@snpFile, getOption("width") - 32),
        "\n             geneAnnotation  : ", if (length(object@geneAnnotation) == 0 || is.null(object@geneAnnotation) || is.na(object@geneAnnotation)) "none" else paste0(truncPath(object@geneAnnotation, getOption("width") - 39), " (", object@geneAnnotationFormat, ")"), "\n", sep = "")
    if (length(object@aligner) == 0 || is.na(object@aligner))
        cat(" Aligner   : unknown\n")
    else
        cat(" Aligner   : ", object@aligner, " v", 
            as.character(utils::packageVersion(object@aligner)),
            " (parameters: ", object@alignmentParameter, ")\n", sep = "")
    cat(" Genome    : ", truncPath(object@genome, getOption("width") - 16 - 
                                       nchar(object@genomeFormat)),
        " (", object@genomeFormat, ")\n", sep = "")
    # reads
    nf <- nrow(object@reads)
    nfs <- if(nf > 1) "s" else ""
    ns <- length(unique(object@reads$SampleName))
    nss <- if(ns > 1) "s" else ""
    cat("\n")
    if (length(object@samplesFormat) == 0) {
        cat(" Reads     : none\n")
    } else if(object@samplesFormat == "bam") {
        cat(" Reads     : none (bam file project)\n")
    } else {
        cat(" Reads     : ", nf, if (object@paired != "no") paste(" pair", nfs,
                                                                  " of files, ", sep = "") 
            else paste(" file", nfs, ", ", sep = ""),
            ns, " sample", nss, " (", object@samplesFormat, " format):\n", sep = "")
        qual <- if(object@samplesFormat == "fastq") paste(" (phred",object@reads$phred, ")", 
                                                          sep = "") else ""
        fw <- max(nchar(basename(unlist(object@reads[, grep("FileName",colnames(object@reads))]))))
        sw <- max(nchar(object@reads$SampleName))
        if (object@paired != "no") {
            fw <- min(fw, floor((getOption("width") - 10 - sw - max(nchar(qual))) / 2))
            cat(sprintf(" %3d. %-*s  %-*s  %-*s%s\n", seq_len(nf),
                        fw, truncString(basename(object@reads$FileName1), fw),
                        fw, truncString(basename(object@reads$FileName2), fw),
                        sw, object@reads$SampleName, qual), sep = "")
        } else {
            fw <- min(fw, getOption("width") - 8 - sw - max(nchar(qual)))
            cat(sprintf(" %3d. %-*s  %-*s%s\n", seq_len(nf), fw, 
                        truncString(basename(object@reads$FileName),fw), 
                        sw, object@reads$SampleName, qual), sep = "")
        }
    }
    # alignments
    cat("\n")
    cat(" Genome alignments: directory: ", if (length(object@alignmentsDir) == 0) "" else if (is.na(object@alignmentsDir)) { if(object@samplesFormat == "bam") "not applicable (bam file project)" else "same as reads" } else truncPath(object@alignmentsDir,getOption("width") - 31), "\n", sep = "")
    #fw <- min(getOption("width") -8 -sw, max(nchar(basename(object@alignments[,"FileName"]))))
    #cat(sprintf(" %3d. %-*s  %-*s\n", 1:nf, fw, truncString(basename(object@alignments[,"FileName"]), fw), sw, object@reads$SampleName), sep="")
    if (length(object@alignments) > 0) {
      tmpmax <- suppressWarnings(max(nchar(basename(object@alignments[, "FileName"])), na.rm = TRUE))
      if (!is.finite(tmpmax)) tmpmax <- Inf
      fw <- min(getOption("width") - 6, tmpmax)
      cat(sprintf(" %3d. %-*s\n", seq_len(nf), fw, 
                  truncString(basename(object@alignments[, "FileName"]), fw)), sep = "")
      cat("\n")
    } 
    # auxiliaries
    cat(" Aux. alignments: ", if(nrow(object@aux) == 0) "none" else paste(nrow(object@aux), " file", if (nrow(object@aux) > 1) "s" else "", ", directory: ", if (is.na(object@alignmentsDir)) "same as reads" else truncPath(object@alignmentsDir, getOption("width") - 38), sep = ""),
        "\n", sep = "")
    if (length(object@aux)>0) {
        for (i in seq.int(nrow(object@aux))) {
            fw <- getOption("width") - 8 - nchar(object@aux[i, 'AuxName'])
            cat(sprintf(" %3s. %-*s  %-s\n", letters[i], fw, 
                        truncPath(object@aux[i, 'FileName'], fw), object@aux[i, 'AuxName']))
            #fw <- min(getOption("width") -10 -sw, max(nchar(basename(unlist(object@auxAlignments[i,])))))
            #cat(sprintf("   %3d. %-*s  %-*s\n", 1:nf, fw, truncString(basename(unlist(object@auxAlignments[i,])), fw), sw, object@reads$SampleName), sep="")
            fw <- min(getOption("width") - 8, 
                      max(nchar(basename(unlist(object@auxAlignments[i, ])))))
            cat(sprintf("   %3d. %-*s\n", seq_len(nf), fw, 
                        truncString(basename(unlist(object@auxAlignments[i,])), fw)), sep = "")
        }
    }
    cat("\n")
    return(invisible(NULL))
})
