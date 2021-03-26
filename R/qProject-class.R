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
qProjectUpdate <- function(proj, quiet = TRUE) {
    mod <- FALSE
    if (!(.hasSlot(proj, "geneAnnotation"))) {
        proj@geneAnnotation <- NA_character_
        if (!quiet) message("Adding geneAnnotation slot")
        mod <- TRUE
    } 
    if (!(.hasSlot(proj, "geneAnnotationFormat"))) {
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

setMethod("genome", signature(x = "qProject"), function(x) {
    y <- x@genome
    attr(y, "genomeFormat") <- x@genomeFormat
    return(y)
})

setGeneric("auxiliaries", function(x) return(NULL))
setMethod("auxiliaries", signature(x = "qProject"), function(x) return(x@aux))

setGeneric("alignments", function(x) return(NULL))
setMethod("alignments", signature(x = "qProject"), function(x) {
    return(list(genome = x@alignments, aux = x@auxAlignments))
})

setMethod("[", signature(x = "qProject", i = "ANY", j = "missing",
                         drop = "missing"), 
          function(x, i) {
              y <- x
              if (is.character(i)) {
                  if (all(i %in% x@reads$SampleName)) {
                      if (any(inu <- sapply(unique(i), 
                                            function(ii) sum(x@reads$SampleName == ii)) > 1))
                          warning(sprintf("will only select first of multiple files with identical sample name: %s", paste(unique(i)[inu], collapse = ", ")))
                      i <- match(i, x@reads$SampleName)
                  } else {
                      stop("undefined samples selected")
                  }
              } else if (is.logical(i) && length(i) != length(x@reads$SampleName)) {
                  stop(sprintf("logical subsetting vector of length %d does not match the number of sequence files (%d)", length(i), length(x@reads$SampleName)))
              } else if(is.numeric(i)) {
                  if (any(!is.finite(i)) || max(i) > length(x@reads$SampleName) || min(i) < 1)
                      stop("undefined samples selected")
              }
              y@reads <- y@reads[i, , drop = FALSE]
              y@reads_md5subsum <- y@reads_md5subsum[i, , drop = FALSE]
              y@alignments <- y@alignments[i, , drop = FALSE]
              if(ncol(y@auxAlignments))
                  y@auxAlignments <- y@auxAlignments[, i, drop = FALSE]
              return(y)
          })

#setGeneric("niceprint", function(x) print(x))
#setMethod("niceprint", "qProject", function(object) {
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
        cat(" Aligner   : ", object@aligner, " v", as.character(packageVersion(object@aligner)),
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
            cat(sprintf(" %3d. %-*s  %-*s  %-*s%s\n", 1:nf,
                        fw, truncString(basename(object@reads$FileName1), fw),
                        fw, truncString(basename(object@reads$FileName2), fw),
                        sw, object@reads$SampleName, qual), sep = "")
        } else {
            fw <- min(fw, getOption("width") - 8 - sw - max(nchar(qual)))
            cat(sprintf(" %3d. %-*s  %-*s%s\n", 1:nf, fw, 
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
      cat(sprintf(" %3d. %-*s\n", 1:nf, fw, 
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
            cat(sprintf("   %3d. %-*s\n", 1:nf, fw, 
                        truncString(basename(unlist(object@auxAlignments[i,])), fw)), sep = "")
        }
    }
    cat("\n")
    return(invisible(NULL))
})
