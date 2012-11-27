### qProject class definition
setClass("qProject",
         representation(reads="data.frame",
                        reads_md5subsum="data.frame",
                        alignments="data.frame",
                        samplesFormat="character",
                        genome="character",
                        genomeFormat="character",
                        aux="data.frame",
                        auxAlignments="data.frame",
                        aligner="character",
                        maxHits="numeric",
                        paired="character",
                        splicedAlignment="logical",
                        snpFile="character",
                        bisulfite="character",
                        alignmentParameter="character",
                        projectName="character",
                        alignmentsDir="character",
                        lib.loc="character",
                        cacheDir="character",
                        alnModeID="character")
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

setMethod("genome", "qProject", function(x) {
    y <- x@genome; attr(y, "genomeFormat") <- x@genomeFormat; return(y)
})

setGeneric("auxiliaries", function(x) return(NULL))
setMethod("auxiliaries", "qProject", function(x) return(x@aux))

setGeneric("alignments", function(x) return(NULL))
setMethod("alignments", "qProject", function(x) return(list(genome=x@alignments, aux=x@auxAlignments)))

setMethod("[", "qProject", function(x, i) {
    y <- x
    if(is.character(i)) {
        if(all(i %in% x@reads$SampleName)) {
            if(any(inu <- sapply(unique(i), function(ii) sum(x@reads$SampleName == ii)) > 1))
                warning(sprintf("will only select first of multiple files with identical sample name: %s",paste(unique(i)[inu], collapse=", ")))
            i <- match(i, x@reads$SampleName)
        } else {
            stop("undefined samples selected")
        }
    } else if(is.logical(i) && length(i) != length(x@reads$SampleName)) {
        stop(sprintf("logical subsetting vector of length %d does not match the number of sequence files (%d)",length(i),length(x@reads$SampleName)))
    } else if(is.numeric(i)) {
        if(any(!is.finite(i)) || max(i) > length(x@reads$SampleName) || min(i) < 1)
            stop("undefined samples selected")
    }
    y@reads <- y@reads[i, , drop=FALSE]
    y@reads_md5subsum <- y@reads_md5subsum[i, , drop=FALSE]
    y@alignments <- y@alignments[i, , drop=FALSE]
    if(ncol(y@auxAlignments))
        y@auxAlignments <- y@auxAlignments[, i, drop=FALSE]
    return(y)
})

#setGeneric("niceprint", function(x) print(x))
#setMethod("niceprint", "qProject", function(object) {
setMethod("show", "qProject", function(object) {
    # projet and global options
    cat("Project: " , object@projectName, "\n", sep="")
    cat(" Options   : maxHits         : ", object@maxHits,
        "\n             paired          : ", object@paired,
        "\n             splicedAlignment: ", object@splicedAlignment,
        "\n             bisulfite       : ", object@bisulfite,
        "\n             snpFile         : ", if(is.na(object@snpFile)) "none" else truncPath(object@snpFile,getOption("width")-32), "\n", sep="")
    if(is.na(object@aligner))
        cat(" Aligner   : unknown\n")
    else
        cat(" Aligner   : ", object@aligner, " v", installed.packages()[object@aligner, 'Version'],
            " (parameters: ", object@alignmentParameter, ")\n", sep="")
    cat(" Genome    : ", truncPath(object@genome, getOption("width")-16-nchar(object@genomeFormat)),
        " (", object@genomeFormat, ")\n", sep="")
    # reads
    nf <- nrow(object@reads)
    nfs <- if(nf>1) "s" else ""
    ns <- length(unique(object@reads$SampleName))
    nss <- if(ns>1) "s" else ""
    cat("\n")
    cat(" Reads     : ", nf, if(object@paired != "no") paste(" pair",nfs," of files, ",sep="") else paste(" file",nfs,", ",sep=""),
        ns, " sample", nss, " (", object@samplesFormat, " format):\n", sep="")
    qual <- if(object@samplesFormat=="fastq") paste(" (phred",object@reads$phred,")",sep="") else ""
    fw <- max(nchar(basename(unlist(object@reads[,grep("FileName",colnames(object@reads))]))))
    sw <- max(nchar(object@reads$SampleName))
    if(object@paired != "no") {
        fw <- min(fw, floor((getOption("width") -10 -sw -max(nchar(qual))) /2))
        cat(sprintf(" %3d. %-*s  %-*s  %-*s%s\n", 1:nf,
                    fw, truncString(basename(object@reads$FileName1),fw),
                    fw, truncString(basename(object@reads$FileName2),fw),
                    sw, object@reads$SampleName, qual), sep="")
    } else {
        fw <- min(fw, getOption("width") -8 -sw -max(nchar(qual)))
        cat(sprintf(" %3d. %-*s  %-*s%s\n", 1:nf, fw, truncString(basename(object@reads$FileName),fw), sw, object@reads$SampleName, qual), sep="")
    }
    # alignments
    cat("\n")
    cat(" Genome alignments: directory: ", if(is.na(object@alignmentsDir)) "same as reads" else truncPath(object@alignmentsDir,getOption("width")-31),"\n", sep="")
    #fw <- min(getOption("width") -8 -sw, max(nchar(basename(object@alignments[,"FileName"]))))
    #cat(sprintf(" %3d. %-*s  %-*s\n", 1:nf, fw, truncString(basename(object@alignments[,"FileName"]), fw), sw, object@reads$SampleName), sep="")
    fw <- min(getOption("width") -6, max(nchar(basename(object@alignments[,"FileName"]))))
    cat(sprintf(" %3d. %-*s\n", 1:nf, fw, truncString(basename(object@alignments[,"FileName"]), fw)), sep="")
    # auxiliaries
    cat("\n")
    cat(" Aux. alignments: ", if(nrow(object@aux)==0) "none" else paste(nrow(object@aux), " file", if(nrow(object@aux)>1) "s" else "", ", directory: ", if(is.na(object@alignmentsDir)) "same as reads" else truncPath(object@alignmentsDir,getOption("width")-38), sep=""), "\n", sep="")
    if(nrow(object@aux)>0) {
        for(i in seq.int(nrow(object@aux))) {
            fw <- getOption("width") -8 -nchar(object@aux[i,'AuxName'])
            cat(sprintf(" %3s. %-*s  %-s\n", letters[i], fw, truncPath(object@aux[i,'FileName'], fw), object@aux[i,'AuxName']))
            #fw <- min(getOption("width") -10 -sw, max(nchar(basename(unlist(object@auxAlignments[i,])))))
            #cat(sprintf("   %3d. %-*s  %-*s\n", 1:nf, fw, truncString(basename(unlist(object@auxAlignments[i,])), fw), sw, object@reads$SampleName), sep="")
            fw <- min(getOption("width") -8, max(nchar(basename(unlist(object@auxAlignments[i,])))))
            cat(sprintf("   %3d. %-*s\n", 1:nf, fw, truncString(basename(unlist(object@auxAlignments[i,])), fw)), sep="")
        }
    }
    cat("\n")
    return(invisible(NULL))
})
