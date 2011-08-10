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

## Again, we the separate function?
.align <- function(readsFilepath, aligner, index, outpath, overwrite=TRUE)
{
    bamFilename <- file.path(outpath,
                             sprintf("%s-%s",
                                     .baseFileName(readsFilepath),
                                     index$name))
    .progressReport(sprintf("Aligning reads to index %s for sample '%s'", index$name, basename(readsFilepath)))
    outputFilename <- switch(aligner$pkgname,
                             Rbowtie = .alignBowtie(readsFilepath, index$path, bamFilename, force=overwrite),
                             Rbwa = .alignBWA(readsFilepath, index$path, bamFilename, force=overwrite),
                             stop("The '", aligner$pkgname, "' Aligner is not supported.")
                             )
    return(outputFilename)
}

.alignBowtie <- function(sequences, index, outfile, ..., indexDestination=FALSE, force=FALSE){
    seqFormat <- ifelse(.fileExtension(sequences) %in% c("fa","fna","mfa","fasta"),
                        seqFormat <- "-f", "")
    numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-p %s", getOption("quasr.cores")), "")
    maxHits <- ifelse(getOption("quasr.maxhits") > 0 , sprintf("-k %s", getOption("quasr.maxhits")), "")
    samFilename <- sprintf("%s.sam", tempfile())
    on.exit(unlink(samFilename))
    bamFilename <- unlist(strsplit(outfile, "\\.bam$"))
    out <- .execute("Rbowtie", paste("bowtie", numThreads, maxHits, index, seqFormat, sequences, "-S", samFilename, "--quiet"))
    outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination, overwrite=force)
    return(outfile)
}

.alignBWA <- function(sequences, index, outfile, ..., indexDestination=FALSE, force=FALSE){
    numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-t %s", getOption("quasr.cores")), "")
    saiFilename <- sprintf("%s.sai", tempfile())
    samFilename <- sprintf("%s.sam", tempfile())
    on.exit(unlink(samFilename))
    bamFilename <- unlist(strsplit(outfile, "\\.bam$"))
    #out <- .execute("Rbwa", paste("bwa bwasw", numThreads, index, sequences, "-f", samFilename))
    out <- .execute("Rbwa", paste("bwa aln", numThreads, index, sequences, "-f", saiFilename))
    out <- .execute("Rbwa", paste("bwa samse", index, saiFilename, sequences, "-f", samFilename))
    outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination, overwrite=force)
    return(outfile)
}