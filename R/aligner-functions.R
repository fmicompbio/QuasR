.execute <- function(object, callstr){
    .requirePkg(object)
    fun <- sprintf("%s:::.execute('%s')", object, callstr)
    return(eval(parse(text=fun)))
}

.index <- function(aligner, fastaFilepath, indexName)
{
    outputpath <- switch(aligner$pkgname,
                         Rbowtie = .indexBowtie(fastaFilepath, indexName),
                         Rbwa = .indexBWA(fastaFilepath, indexName),
                         stop("The '", aligner$pkgname, "' Aligner is not supported.")
                         )
    return(outputpath)
}

.indexBowtie <- function(references, outdir){
    references <- paste(references, collapse=",")
    out <- .execute("Rbowtie", sprintf("bowtie-build %s %s", references, outdir))
    return(outdir)
}

.indexBWA <- function(references, outdir){
    references <- .multiToSingleFasta(references)
    out <- .execute("Rbwa", sprintf("bwa index -p %s %s", outdir, references))
    return(outdir)
}

.createAlignmentParameters <- function(qproject){
    alignmentParameters <- switch(qproject@aligner$pkgname,
                             Rbowtie = .createBowtieAlignmentParameters(qproject),
                             Rbwa = .createBwaAlignmentParameters(qproject),
                             stop("The '", aligner$pkgname, "' Aligner is not supported.")
                             )
    return(alignmentParameters)
}

.createBowtieAlignmentParameters <- function(qproject){
    seqFormat <- ifelse(all(qproject@samples$filetype %in% c("fasta", "bam")), seqFormat <- "-f", "-q")
    maxHits <- ifelse(is.null(qproject@maxHits), "", sprintf("-v 2 -k %s -m %s --best --strata", qproject@maxHits, qproject@maxHits))
    #samFilename <- sprintf("-S %s", tempfile(pattern=.baseFileName(qproject@samples$filepath), fileext=".sam"))
    alignmentParameters <- paste("bowtie", maxHits, seqFormat)
    return(alignmentParameters)
}

.createBwaAlignmentParameters <- function(qproject){
    maxHits <- ifelse(is.null(qproject@maxHits), "", sprintf("-n %s", qproject@maxHits-1))
    #saiFilename <- sprintf("%s.sai", tempfile())
    #samFilename <- sprintf("%s.sam", tempfile())
    #out <- .execute("Rbwa", paste("bwa bwasw", numThreads, index, sequences, "-f", samFilename))
    alignmentParameters <- c(paste("bwa aln"),
                             paste("bwa samse", maxHits))
    return(alignmentParameters)
}

.align <- function(readsFilepath, index, qproject)
{
#     bamFilename <- file.path(outpath,
#                              sprintf("%s-%s-%s",
#                                      .baseFileName(readsFilepath),
#                                      index$name,
#                                      aligner$pkgname))
    bamFilename <- .createBamFilename(qproject@path, readsFilepath, index$shortname)
    .progressReport(sprintf("Aligning reads to index %s-%s for sample '%s'", 
                            index$shortname, qproject@aligner$pkgname, basename(readsFilepath)))
    outputFilename <- switch(qproject@aligner$pkgname,
                             Rbowtie = .alignBowtie(readsFilepath, index$path, bamFilename, 
                                                    alignmentParameter=qproject@alignmentParameter),
                             Rbwa = .alignBWA(readsFilepath, index$path, bamFilename, 
                                              alignmentParameter=qproject@alignmentParameter),
                             stop("The '", qproject@aligner$pkgname, "' Aligner is not supported.")
                             )
    .weightAlignments(outputFilename, outputFilename, index=index, qproject=qproject, overwrite=TRUE) #TODO overwrite as qproject parameter
    return(outputFilename)
}

.alignBowtie <- function(sequences, index, outfile, alignmentParameter, ..., indexDestination=FALSE){
    numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-p %s", getOption("quasr.cores")), "")
    samFilename <- sprintf("%s.sam", tempfile())
    on.exit(unlink(samFilename))
    bamFilename <- unlist(strsplit(outfile, "\\.bam$"))
#     out <- .execute("Rbowtie", paste("bowtie -k 99 -m 99 --best --strata -f", numThreads, "--quiet", index, sequences, "-S", samFilename))
    out <- .execute("Rbowtie", paste(alignmentParameter, numThreads, "--quiet", index, sequences, "-S", samFilename))
    outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination)
    return(outfile)
}

.alignBWA <- function(sequences, index, outfile, alignmentParameter, ..., indexDestination=FALSE, force=FALSE){
    numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-t %s", getOption("quasr.cores")), "")
#     maxHits <- ifelse(is.null(maxHits), "", sprintf("-n %s", maxHits-1))
    saiFilename <- sprintf("%s.sai", tempfile())
    samFilename <- sprintf("%s.sam", tempfile())
    on.exit(unlink(samFilename))
    bamFilename <- unlist(strsplit(outfile, "\\.bam$"))
    #out <- .execute("Rbwa", paste("bwa bwasw", numThreads, index, sequences, "-f", samFilename))
#     out <- .execute("Rbwa", paste("bwa aln", numThreads, index, sequences, "-f", saiFilename))
#     out <- .execute("Rbwa", paste("bwa samse", maxHits, index, saiFilename, sequences, "-f", samFilename))
    out <- .execute("Rbwa", paste( alignmentParameter[1], numThreads, index, sequences, "-f", saiFilename))
    out <- .execute("Rbwa", paste( alignmentParameter[2], index, saiFilename, sequences, "-f", samFilename))
    outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination, overwrite=force)
    return(outfile)
}
