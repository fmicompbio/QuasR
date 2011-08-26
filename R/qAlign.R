qAlign <- function(qProject, lib=NULL, lib.loc=NULL)
{
    .progressReport("Starting alignments to genome", phase=-1)
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")

    ## load genome index
    qProject@index <- .loadIndex(qProject, lib=lib, lib.loc=lib.loc)

    ## align to genome
    qProject@alignments$genome <- unlist(lapply(qProject@samples$filepath,
                                                .align,
                                                qProject@aligner,
                                                qProject@index,
                                                qProject@path,
                                                maxHits=qProject@maxHits))
    ## get unmapped reads
    genomeAlignmentUnmapped <- lapply(qProject@alignments$genome,
                                      .unmappedToFasta)
    on.exit(unlink(genomeAlignmentUnmapped))
    ## TODO align to exon-junction-db
    ## qProject@alignments$exonjunction <- unlist(lapply(qProject@samples$filepath,
    #                                           .align,
    #                                           qProject@aligner,
    #                                           exonjunctionIndex,
    #                                           qProject@path))
    ## TODO get unmapped reads
    #exonjunctionAlignmentUnmapped <- lapply(qProject@alignments$exonjunction,
    #                                  .unmappedToFasta)
    
    ## create index and align to annotation
    .progressReport("Creating index of auxiliaries")
    auxIndexes <- .createAuxiliaryIndex(qProject)
    on.exit(unlink(auxIndexes$path))
    
    .progressReport("Starting alignments to auxiliaries")
    auxAlignment <- lapply(auxIndexes, function(auxIndex){
        unlist(lapply(genomeAlignmentUnmapped,
                      .align,
                      qProject@aligner,
                      auxIndex,
                      qProject@path,
                      maxHits=qProject@maxHits))
    })
    qProject@alignments <- cbind.data.frame(qProject@alignments, auxAlignment, stringsAsFactors=FALSE)

    .progressReport("Weight alignments")
    lapply(t(qProject@alignments),
           function(elem){
               .weightAlignments(elem, elem, qProject@maxHits, overwrite=TRUE)
           })
    .progressReport("Successfully terminated the quasr alignment.", phase=1)
    return(qProject)
}

.unmappedToFasta <- function(bamFile, destFile){
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE),
                          what=c("qname", "seq", "qual"))
    unmapped <- scanBam(bamFile, param=param)
    names(unmapped[[1]]$seq) <- unmapped[[1]]$qname

    if(length(IRanges::unique(unmapped[[1]]$qual)) <= 1L){
        if(missing(destFile))
            destfile <- file.path(tempdir(), sprintf("%s-unmapped.fasta", .baseFileName(bamFile)))
        write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fasta")
    } else {
        if(missing(destFile))
            destfile <- file.path(tempdir(), sprintf("%s-unmapped.fastq", .baseFileName(bamFile)))
        write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fastq", qualities=unmapped[[1]]$qual)
    }
    return(destfile)
}