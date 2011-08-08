loadAligner <- function(pkgname, bisulfiteCoversion=FALSE, ...){
    .progressReport("Loading aligner")
    .requirePkg(pkgname, ...)
    aligner <- list(pkgname=pkgname,
                    pkgversion=installed.packages()[pkgname, 'Version'],
                    bisulfiteCoversion=bisulfiteCoversion
                    )
    return(aligner)
}

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
###

qAlign <- function(qProject, lib=NULL, ...)
{
    .progressReport("Starting alignments to genome", phase=-1)
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")

    ## load genome index
    qProject@index <- loadIndex(qProject, lib=lib, ...)

    ## align to genome
    qProject@alignments$genome <- unlist(lapply(qProject@samples$filepath,
                                                   .align,
                                                   qProject@aligner,
                                                   qProject@index,
                                                   qProject@path))
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
    auxIndexes <- createAuxiliaryIndex(qProject)
    #on.exit(unlink(auxIndexes$path))
    .progressReport("Starting alignments to auxiliaries")
    auxAlignment <- lapply(auxIndexes, function(auxIndex){
        unlist(lapply(genomeAlignmentUnmapped,
                      .align,
                      qProject@aligner,
                      auxIndex,
                      qProject@path))
    })
    qProject@alignments <- cbind.data.frame(qProject@alignments, auxAlignment, stringsAsFactors=FALSE)

    .progressReport("Count alignments")
    lapply(t(qProject@alignments),
           function(elem){
               countAlignments(elem, elem, overwrite=TRUE)
           })
    .progressReport("Successfully terminated the hierarchical alignment.", phase=1)
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
###

loadIndex <- function(qProject, lib=NULL, lib.loc=NULL, ...)
{
    .progressReport("Loading genome index")
    ## TODO check if genome and aligner is not NULL
    if(qProject@genome$bsgenome)
    {
        #Genome is BSgenome
        pkgname <- paste(qProject@aligner$pkgname,
                         qProject@genome$name,
                         sep=".")
        if(suppressWarnings(require(pkgname, character.only=TRUE, quietly=TRUE)))
        {
            index <- eval(parse(text=ls(sprintf("package:%s", pkgname))))
        } else {
            .progressReport("No index found. Create index now")
            ## create and install index
            srcPkgDir <- createIndexPackage(qProject)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE)
            index <- get(ls(sprintf("package:%s", pkgname))[1])
        }
    } else {
        indexDir <- file.path(qProject@genome$dir, sprintf("%sIndex", qProject@aligner$pkgname))
        if(!file.exists(indexDir)){
            .progressReport("No index found. Create index now")
            index <- createGenomeIndex(qProject)
        } else {
            index <- read.table(file=file.path(indexDir,"index.tab"), sep="\t", header=TRUE)
        }
    }
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to create the index
###

createGenomeIndex <- function(qProject, destDir)
{
    if(missing(destDir))
        destDir <- qProject@genome$dir
    .progressReport("Creating index")
    fastaFiles <- file.path(qProject@genome$dir, qProject@genome$files)
    index <- .createIndex(fastaFiles, qProject@aligner, qProject@genome$name, destDir)
    return(index)
}

createAuxiliaryIndex <- function(qProject)
{
    .progressReport("Creating auxiliary index")
    isFastaFormat <- .fileExtension(qProject@annotations$filepath) %in% c("fa","fna","mfa","fasta")
    indexes <- mapply(.createIndex,
                      fastaFiles=qProject@annotations[isFastaFormat,]$filepath,
                      name=as.character(qProject@annotations[isFastaFormat,]$feature),
                      MoreArgs=list(aligner=qProject@aligner),
                      SIMPLIFY=FALSE,
                      USE.NAMES = TRUE)
    names(indexes) <- .baseFileName(names(indexes)) ## FIXME: maybe use name or name should be basefilename
    return(indexes)
}

.createIndex <- function(fastaFiles, aligner, name, destDir=tempfile())
{
    if(missing(name))
        name <- .baseFileName(fastaFiles)
    indexName <- file.path(destDir, sprintf("%sIndex", aligner$pkgname), name)
    dir.create(dirname(indexName), showWarnings=FALSE, recursive=TRUE)
    .index(aligner, fastaFiles, indexName)
    index <- list(name=paste(name, aligner$pkgname, sep="-"),
                  path=indexName,
                  aligner=aligner$pkgname,
                  alignerversion=aligner$pkgversion,
                  organism=name,
                  sourceurl=paste(fastaFiles, collapse=",")
                  )
    write.table(index, file=file.path(destDir, sprintf("%sIndex", aligner$pkgname), "index.tab"), sep="\t", col.names=TRUE, row.names=FALSE)
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All Function to create an index package of a genome for a specific aligner
###

createIndexPackage <- function(qProject, sourcePackageFilepath=tempdir())
{
    .progressReport("Creating index package")
    .requirePkg(qProject@aligner$pkgname)
    suppressPackageStartupMessages(require(Biobase, quietly=TRUE))
    genome <- loadBSgenome(qProject@genome$name)
    dir.create(sourcePackageFilepath, showWarnings=FALSE, recursive=TRUE)
    fastaFilepath <- .BSgenomeSeqToFasta(genome)
    seedList <- .createSeedList(genome, qProject@aligner)
    ## create package
    pkgname <- seedList$ALIGNERINDEXNAME
    destinationDir <- sourcePackageFilepath
    templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
    pkgDir <- Biobase::createPackage(pkgname, destinationDir, templatePath, seedList, quiet=TRUE)
    ## create index files
    indexDir <- file.path(pkgDir, "inst", "alignerIndex", seedList$GENOMENAME)
    .index(qProject@aligner, fastaFilepath, indexDir)
    return(pkgDir$pkgdir)
}

.BSgenomeSeqToFasta <- function(bsGenome, outFile=NULL)
{
    if(!is(bsGenome, "BSgenome"))
        stop("The variable bsGenome is not a BSgenome")
    if(missing(outFile) || is.null(outFile)){
        td <- tempdir()
        outFile <- file.path(td, paste(bsgenomeName(bsGenome),"fa",sep="."))
    }
    append <- FALSE
    for(chrT in seqnames(bsGenome)){
        if(is.null(masks(bsGenome[[chrT]])))
            chrSeq <- DNAStringSet(bsGenome[[chrT]])
        else
            chrSeq <- DNAStringSet(unmasked(bsGenome[[chrT]]))
        names(chrSeq) <- chrT
        ## FH: I wonder why filepath can't be a connection. Would have made this appending nightmare a lot cleaner...
        write.XStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
        append <- TRUE
    }
    return(outFile)
}

.installIndexPackage <- function(pkgdir, lib=NULL)
{
    curwd <- setwd(dirname(pkgdir))
    on.exit(setwd(curwd))
    .progressReport("Installing new index")
    buildFun <- tools:::.build_packages
    fbody <- body(buildFun)
    body(buildFun) <- fbody[1:(length(fbody)-1)]
    out <- capture.output(buildFun(basename(pkgdir)))
    newPack <- dir(pattern="\\.tar.gz$")
    if(is.null(lib))
        lib <- .libPaths()[1]
    out <- c(out, capture.output(install.packages(file.path(getwd(), newPack), repos=NULL, dependencies=FALSE, lib=lib)))
    return(out)
}

.defaultSeedList <- function(genome, aligner)
{
    seed <- list(#package seeds
                 PKGTITLE="",
                 PKGDESCRIPTION="",
                 PKGVERSION="0.1.1",
                 DATE=format(Sys.time(), "%Y-%M-%d"),
                 AUTHOR="Who wrote it",
                 MAINTAINER="Who to complain to <yourfault@somewhere.net>",
                 LIC="What license is it under?",
                 PKGDETAILS="An overview of how to use the package, including the most important functions.",
                 PKGEXAMPLES="Simple examples of the most important functions.",
                 ##genome seeds
                 GENOMENAME="",
                 PROVIDER="",
                 PROVIDERVERSION="",
                 RELEASEDATE="",
                 RELEASENAME="",
                 ORGANISM="",
                 SPECIES="",
                 SRCDATAFILES="",
                 ORGANISMBIOCVIEW="",
                 ##aligner seeds
                 ALIGNER=aligner$pkgname,
                 ALIGNERVERSION=aligner$pkgversion,
                 ALIGNERINDEXNAME=""
                 #OBJECTNAME=""
                 )
    if(is(genome,"BSgenome")){
        seed$GENOMENAME <- bsgenomeName(genome)
        seed$PROVIDER <- provider(genome)
        seed$PROVIDERVERSION <- providerVersion(genome)
        seed$RELEASEDATE <- releaseDate(genome)
        seed$RELEASENAME <- releaseName(genome)
        seed$ORGANISM <- organism(genome)
        seed$SPECIES <- species(genome)
        seed$SRCDATAFILES <- bsgenomeName(genome)
        seed$ORGANISMBIOCVIEW <- gsub(" ", "_", organism(genome))
    }else{
        seed$GENOMENAME <- genome$name #.baseFileName(genome)
        seed$SRCDATAFILES <- genome$dir
    }
    seed$PKGTITLE <- paste(seed$ALIGNER, "Index of", seed$GENOME)
    seed$PKGDESCRIPTION <- seed$PKGTITLE
    seed$ALIGNERINDEXNAME <- paste(seed$ALIGNER, gsub("_", "", seed$GENOME), sep=".")
    return(seed)
}

.createSeedList <- function(genome, aligner, pkgtitle, pkgdescription,
                            pkgversion, author, maintainer, license)
{
    seed <- .defaultSeedList(genome, aligner)
    if(!missing(pkgtitle))
        seed$PKGTITLE <- pkgtitle
    if(!missing(pkgdescription))
        seed$PKGDESCRIPTION <- pkgdescription
    if(!missing(pkgversion))
        seed$PKGVERSION <- pkgversion
    if(!missing(author))
        seed$AUTHOR <- author
    if(!missing(maintainer))
        seed$MAINTAINER <- maintainer
    if(!missing(license))
        seed$LIC <- license
    ##  seed$PKGDETAILS <- if(!missing(pkgdetails)) pkgdetails
    ##  seed$PKGEXAMPLES <- if(!missing(pkgexamples)) pkgexamples
    ##  if(!missing(indexname))
    ##    seed$ALIGNERINDEXNAME <- indexname
    return(seed)
}

