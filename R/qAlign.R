qAlign <- function(qproject, lib=NULL, lib.loc=NULL)
{
    .progressReport("Starting alignments to genome", phase=-1)
    if(!is(qproject, "qProject"))
        stop("The object '", class(qproject), "' is not a 'qProject' object.")
    if(any(idx <- is.na(qproject@env$alignments$genome))){ 
        ## load genome index
        qproject@env$index$genome <- .loadIndex(qproject, lib=lib, lib.loc=lib.loc)
        ## align reads to genome
        genomeAlignments <- lapply(qproject@env$samples$filepath[idx],
                                                .align,
                                                index=qproject@env$index$genome,
                                                qproject=qproject)
        names(genomeAlignments) <- basename(qproject@env$samples$filepath[idx])
        qproject@env$alignments$genome[idx] <- unlist(lapply(genomeAlignments, "[[", "bamfile"))
        qproject@env$qc$mappingStats$genome <- t(as.data.frame(lapply(genomeAlignments, "[[", "mappingStats")))
        ## TODO get unmapped reads
        ## TODO align to exon-junction-db 
    }

    ## create index and align unmapped reads to annotation
    if(!is.null(qproject@env$annotations) && any(is.na(qproject@env$alignments))){
        ## get unmapped reads
        .progressReport("Get unmapped reads")
        unmapped <- unlist(lapply(qproject@env$alignments$genome, .unmappedToFasta))
        on.exit(unlink(unmapped))        
        ## create auxiliary index
        .progressReport("Creating index of auxiliaries")
        auxIndexes <- .createAuxiliaryIndex(qproject)
#         on.exit(unlink(auxIndexes$path)) ## TODO delete aux index file in the end
        qproject@env$index[names(auxIndexes)] <- auxIndexes
        auxIndexes <- names(auxIndexes)
        ## align unmapped reads auxiliary index
        .progressReport("Starting alignments to auxiliaries")
        lapply(auxIndexes, function(auxIndex){
            idx <- is.na(qproject@env$alignments[auxIndex])
            auxAlign <- lapply(unmapped[idx],
                                  .align,
                                  index=qproject@env$index[[auxIndex]],
                                  qproject=qproject)
            names(auxAlign) <- basename(qproject@env$samples$filepath[idx])
            qproject@env$alignments[idx, auxIndex] <- unlist(lapply(auxAlign, "[[", "bamfile"))
            qproject@env$qc$mappingStats[[auxIndex]] <- t(as.data.frame(lapply(auxAlign, "[[", "mappingStats")))
        })
    }

    # sort and index external bamfile
    if(any(idx <- qproject@env$samples$filetype == "bam")){
        .progressReport("Sort and index external bamfile")
        bamfile <- sortBam(as.character(qproject@env$samples$filepath[idx]), tempfile(tmpdir=qproject@env$cacheDir))
        file.rename(bamfile, qproject@env$samples$filepath[idx])
        indexBam(as.character(qproject@env$samples$filepath[idx]))
    }
#     .progressReport("Weight alignments")
#     lapply(t(qproject@env$alignments$genome),
#            function(elem){
#                .weightAlignments(elem, elem, qproject=qproject, overwrite=TRUE) #TODO overwrite as qproject parameter
#            })
    qproject@env$index <- qproject@env$index["genome"] ## remove aux indexes
    .progressReport("Successfully finished the quasr alignment.", phase=1)
    return(qproject)
}

.unmappedToFasta <- function(bamFile, destFile){
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE),
                          what=c("qname", "seq", "qual"))
    unmapped <- scanBam(bamFile, param=param)
    names(unmapped[[1]]$seq) <- unmapped[[1]]$qname

    if(length(IRanges::unique(unmapped[[1]]$qual)) <= 1L){
        if(missing(destFile))
            destfile <- file.path(tempdir(), sprintf("%s.fasta", .baseFileName(bamFile)))
        write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fasta")
    } else {
        if(missing(destFile))
            destfile <- file.path(tempdir(), sprintf("%s.fastq", .baseFileName(bamFile)))
        write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fastq", qualities=unmapped[[1]]$qual)
    }
    return(destfile)
}