qAlign <- function(qproject, lib=NULL, lib.loc=NULL)
{
    .progressReport("Starting alignments to genome", phase=-1)
    if(!is(qproject, "qProject"))
        stop("The object '", class(qproject), "' is not a 'qProject' object.")
    if(any(idx <- is.na(qproject@alignments$genome))){
        ## load genome index
        qproject@index <- .loadIndex(qproject, lib=lib, lib.loc=lib.loc)
        ## align reads to genome
        genomeAlignments <- lapply(qproject@samples$filepath[idx],
                                                .align,
                                                index=qproject@index,
                                                qproject=qproject)
        names(genomeAlignments) <- basename(qproject@samples$filepath[idx])
        qproject@alignments$genome[idx] <- unlist(lapply(genomeAlignments, "[[", "bamfile"))
        qproject@qc$mappingStats$genome <- t(as.data.frame(lapply(genomeAlignments, "[[", "mappingStats")))
#         qproject@alignments$genome[idx] <- unlist(lapply(qproject@samples$filepath[idx],
#                                                 .align,
#                                                 index=qproject@index,
#                                                 qproject=qproject))
        ## TODO get unmapped reads
        ## TODO align to exon-junction-db 
    }

    ## create index and align unmapped reads to annotation
    if(!is.null(qproject@annotations) && any(is.na(qproject@alignments))){
        ## get unmapped reads
        .progressReport("Get unmapped reads")
        unmapped <- unlist(lapply(qproject@alignments$genome, .unmappedToFasta))
        on.exit(unlink(unmapped))        
        ## create auxiliary index
        .progressReport("Creating index of auxiliaries")
        auxIndexes <- .createAuxiliaryIndex(qproject)
        names(auxIndexes) <- qproject@annotations$feature[qproject@annotations$filetype == "fasta"]
        on.exit(unlink(auxIndexes$path))
        ## align unmapped reads auxiliary index
        .progressReport("Starting alignments to auxiliaries")
        auxAlignment <- lapply(auxIndexes, function(auxIndex){ # TODO why index creation not in this loop
            idx <- is.na(qproject@alignments[auxIndex$shortname])
            auxAlign <- lapply(unmapped[idx],
                                  .align,
                                  index=auxIndex,
                                  qproject=qproject)
            names(auxAlign) <- basename(qproject@samples$filepath[idx])
            qproject@alignments[idx, auxIndex$shortname] <- unlist(lapply(auxAlign, "[[", "bamfile"))
            qproject@qc$mappingStats[[auxIndex$shortname]] <- t(as.data.frame(lapply(auxAlign, "[[", "mappingStats")))         
#             qproject@alignments[idx, auxIndex$shortname] <- unlist(lapply(unmapped[idx],
#                           .align,
#                           index=auxIndex,
#                           qproject=qproject))
            return(qproject@alignments[, auxIndex$shortname])
        })
        qproject@alignments <- cbind.data.frame(genome=qproject@alignments$genome, auxAlignment, stringsAsFactors=FALSE)
    }

    # sort and index external bamfile
    if(any(idx <- qproject@samples$filetype == "bam")){
        .progressReport("Sort and index external bamfile")
        bamfile <- sortBam(as.character(qproject@samples$filepath[idx]), tempfile(tmpdir=qproject@cacheDir))
        file.rename(bamfile, qproject@samples$filepath[idx])
        indexBam(as.character(qproject@samples$filepath[idx]))
    }
#     .progressReport("Weight alignments")
#     lapply(t(qproject@alignments),
#            function(elem){
#                .weightAlignments(elem, elem, qproject=qproject, overwrite=TRUE) #TODO overwrite as qproject parameter
#            })
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