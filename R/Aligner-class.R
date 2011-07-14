# Aligner Class
setClass("Aligner", representation(pkgname="character",
                                   pkgversion="character"))

setMethod("initialize", "Aligner", function(.Object, pkgname)
      {
          .requirePkg(pkgname)
          pkgname <- pkgname
          pkgversion <- installed.packages()[pkgname, 'Version']
          callNextMethod(.Object, pkgname=pkgname, pkgversion=pkgversion)
      })

Aligner <- function(aligner) new("Aligner", aligner)


setMethod("show","Aligner", function(object) cat(object@pkgname, "Aligner Version", object@pkgversion, "\n"))

setGeneric("execute", function(object, callstr) standardGeneric("execute"))

setMethod("execute", "Aligner", function(object, callstr)
      {
          .requirePkg(object@pkgname) 
          fun <- sprintf("%s:::.execute('%s')", object@pkgname, callstr)
          return(eval(parse(text=fun)))
      })

setMethod("execute", "character", function(object, callstr)
      {
          .requirePkg(object) 
          fun <- sprintf("%s:::.execute('%s')", object, callstr)
          return(eval(parse(text=fun)))
      })

## TODO display the usage page of the aligner
##help <- function(alinger) return(character())

.index <- function(aligner, fastaFilepath, indexName)
{
  outputpath <- switch(aligner@pkgname,
                Rbowtie = .indexBowtie(fastaFilepath, indexName),
                Rbwa = .indexBWA(fastaFilepath, indexName),
                stop("The '", aligner@pkgname, "' Aligner is not supported.")
                )
  return(outputpath)
}

.indexBowtie <- function(references, outdir){
  references <- paste(references, collapse=",")
  out <- execute("Rbowtie", sprintf("bowtie-build %s %s", references, outdir))
  return(outdir)
}

.indexBWA <- function(references, outdir){
  references <- .multiToSingleFasta(references)
  out <- execute("Rbwa", sprintf("bwa index -p %s %s", outdir, references))
  return(outdir)
}

## Again, we the separate function?
.align <- function(readsFilepath, aligner, index, outpath, overwrite=TRUE)
{
    bamFilename <- file.path(outpath,
                                sprintf("%s-%s",
                                        .baseFileName(readsFilepath),
                                        index$name))   
    .progressReport(sprintf("\nAligning reads to index %s for sample '%s'", index$name, basename(readsFilepath)), phase=-1)
    outputFilename <- switch(aligner@pkgname,
           Rbowtie = .alignBowtie(readsFilepath, index$path, bamFilename, force=overwrite),
           Rbwa = .alignBWA(readsFilepath, index$path, bamFilename, force=overwrite),
           stop("The '", aligner@pkgname, "' Aligner is not supported.")
           )
    return(outputFilename)
}

.alignBowtie <- function(sequences, index, outfile, ..., indexDestination=FALSE, force=FALSE){
  seqFormat <- ifelse(.fileExtension(sequences) %in% c("fa","fna","mfa","fasta"),
                      seqFormat <- "-f", "")
  numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-p %s", getOption("quasr.cores")), "")
  maxHits <- ifelse(getOption("quasr.maxhits") > 0 , sprintf("-k %s", getOption("quasr.maxhits")), "")
  samFilename <- sprintf("%s.sam", tempfile())
  bamFilename <- unlist(strsplit(outfile, "\\.bam$")) 
  out <- execute("Rbowtie", paste("bowtie", numThreads, maxHits, index, seqFormat, sequences, "-S", samFilename))
  outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination, overwrite=force)
  return(outfile)
}

.alignBWA <- function(sequences, index, outfile, ..., indexDestination=FALSE, force=FALSE){
  numThreads <- ifelse(getOption("quasr.cores") > 1 , sprintf("-t %s", getOption("quasr.cores")), "")
  samFilename <- sprintf("%s.sam", tempfile())
  bamFilename <- unlist(strsplit(outfile, "\\.bam$"))  
  out <- execute("Rbwa", paste("bwa bwasw", numThreads, index, sequences, ">", samFilename))
  outfile <- asBam(samFilename, bamFilename, indexDestination=indexDestination, overwrite=force)
  return(outfile)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 
###

alignHierarchical <- function(projectInfo)
{
    if(!is(projectInfo, "ProjectInfo"))
        stop("The object '", class(projectInfo), "' is not a 'ProjectInfo' object.")
    .progressReport("Starting alignments to genome", phase=-1)
    ## align to genome
    projectInfo@alignments$genome <- unlist(lapply(projectInfo@samples$filepath,
                                               .align,
                                               projectInfo@aligner,
                                               projectInfo@index,
                                               projectInfo@path))
    ## get unmapped reads
    genomeAlignmentUnmapped <- lapply(projectInfo@alignments$genome,
                                      .unmappedToFasta)
    ## TODO align to exon-junction-db
    ## projectInfo@alignments$exonjunction <- unlist(lapply(projectInfo@samples$filepath,
    #                                           .align,
    #                                           projectInfo@aligner,
    #                                           exonjunctionIndex,
    #                                           projectInfo@path))
    ## TODO get unmapped reads
    #exonjunctionAlignmentUnmapped <- lapply(projectInfo@alignments$exonjunction,
    #                                  .unmappedToFasta)
    ## create index and align to annotation
    .progressReport("Creating index of auxiliaries")
    auxIndexes <- createAuxiliaryIndex(projectInfo)
    .progressReport("Starting alignments to auxiliaries")  
    auxAlignment <- lapply(auxIndexes, function(auxIndex){
      unlist(lapply(genomeAlignmentUnmapped,
                    .align,
                    projectInfo@aligner,
                    auxIndex,
                    projectInfo@path))
    })
    projectInfo@alignments <- cbind(projectInfo@alignments, auxAlignment)

    #TODO calc weights
    browser()
    # TODO all alignments with same filename
    projectInfo@alignments[1,1] <- countAlignments(projectInfo@alignments[1,1], projectInfo@alignments[1,1], overwrite=T)
    ##Testcode
    ##param <- ScanBamParam(what=c("qname","rname"), tag="IH")
    ##list <- scanBam(projectInfo@alignments[1,1], param=param)
    ##list[[1]]$tag$IH[list[[1]]$tag$IH>1]
    ##summary(list[[1]]$tag$IH)
    .progressReport("Successfully terminated the hierarchical alignment.", phase=1)
    return(projectInfo)
p}

.unmappedToFasta <- function(bamFile, destFile){
    ## get unmapped # list of qname in fasta minus qname in bam
    #param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), what="qname")
    #mapped <- scanBam(bamFilename, param=param)[[1]]$qname
    #fa <- open(FaFile(fastaFile))
    #faidx <- scanFaIndex(fa)
    #unmapped <- seqnames(faidx) != mapped
    #tempFasta <- tempfile(tmpdir=tempdir())
    #write.XStringSet(scanFa(fa,faidx[unmapped]), file=tempFasta, format="fasta")
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE),
                        what=c("qname", "seq", "qual"))
  unmapped <- scanBam(bamFile, param=param)
  names(unmapped[[1]]$seq) <- unmapped[[1]]$qname
  #if(length(unique(unmapped[[1]]$qual)) == 1L){
    if(missing(destFile))
      destfile <- file.path(tempdir(), sprintf("%s-unmapped.fasta", .baseFileName(bamFile)))
    write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fasta")
  # TODO
  #} else {
  #  if(missing(destFile))
  #    destfile <- file.path(tempdir(), sprintf("%s-unmapped.fastq", .baseFileName(bamFile)))
  #  write.XStringSet(unmapped[[1]]$seq, file=destfile, format="fastq", qualities=unmapped[[1]]$qual) 
  #}
  return(destfile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 
###

loadIndex <- function(genome, aligner, lib=NULL)
{
    .progressReport("Loading aligner index")
    if(is(genome,"BSgenome"))
    {
        pkgname <- paste(aligner@pkgname, bsgenomeName(genome), sep=".")
        if(suppressWarnings(require(pkgname, character.only=TRUE, quietly=TRUE)))
        {
            index <- eval(parse(text=ls(sprintf("package:%s", pkgname))))
        } else {
            .progressReport("No index found. Creating one now")
            srcPkgDir <- createIndexPackage(genome, aligner)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE)
            index <- get(ls(sprintf("package:%s", pkgname))[1]) 
        }
    } else {
        indexDir <- file.path(genome$dir, sprintf("%sIndex", aligner@pkgname))
        if(!file.exists(indexDir))
          index <- createGenomeIndex(genome, aligner)
        else{ #TODO what if the index is there, how to fix
          name <- paste(aligner@pkgname, gsub("_", "", genome$name), sep=".") 
          index <- list(name=name,
                        path=indexDir,
                        aligner=aligner@pkgname,
                        organism=genome$name)
        }
    }
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to create the index
###

createGenomeIndex <- function(genome, aligner, destDir)
{
  if(missing(destDir))
    destDir <- genome$dir
  .progressReport("Creating index")
  fastaFiles <- file.path(genome$dir, genome$files)
  index <- createIndex(fastaFiles, aligner, genome$name, destDir)
  #indexName <- file.path(destDir, sprintf("%sIndex", aligner@pkgname), genome$name)
  #dir.create(dirname(indexName), showWarnings=FALSE, recursive=TRUE)
  #.index(aligner, fastaFiles, indexName)
  #return(indexName)
  return(index)
}

createAuxiliaryIndex <- function(projectInfo)
{
  .progressReport("Creating auxiliary index")
  isFastaFormat <- .fileExtension(projectInfo@annotations$filepath) %in% c("fa","fna","mfa","fasta")
  indexes <- mapply(createIndex,
                    fastaFiles=projectInfo@annotations[isFastaFormat,]$filepath,
                    name=as.character(projectInfo@annotations[isFastaFormat,]$feature),
                    MoreArgs=list(aligner=projectInfo@aligner),
                    SIMPLIFY=FALSE,
                    USE.NAMES = TRUE)
  names(indexes) <- .baseFileName(names(indexes)) ## FIXME: maybe use name or name should be basefilename
  return(indexes)         
}

createIndex <- function(fastaFiles, aligner, name, destDir=tempfile())
{
  if(missing(name))
    name <- .baseFileName(fastaFiles)
  indexName <- file.path(destDir, sprintf("%sIndex", aligner@pkgname), name)
  dir.create(dirname(indexName), showWarnings=FALSE, recursive=TRUE)
  .index(aligner, fastaFiles, indexName)
  ## create index readme file
  logFile <- file(paste(indexName, "readme.txt", sep="."), "w")
  cat("Name:", name, "\n", file=logFile)
  cat("Index base path:", indexName, "\n", file=logFile)
  cat("Aligner:", aligner@pkgname, "\n", file=logFile)
  cat("Alinger version:", aligner@pkgversion, "\n", file=logFile)
  cat("Creation date:", date(), "\n", file=logFile)  
  cat("Organism:", name, "\n", file=logFile)
  cat("Source:", fastaFiles, file=logFile, sep="\n")  
  close(logFile)
  index <- list(name=paste(name, aligner@pkgname, sep="-"),
                path=indexName,
                aligner=aligner@pkgname,
                alignerversion=aligner@pkgversion,
                organism=name,
                sourceurl=fastaFiles
                )
  #write.table(index, file = paste(indexName,"index.tab", sep="."), row.names = FALSE)
  #return(indexName)
  return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All Function to create an index package of a genome for a specific aligner
###

setGeneric("createIndexPackage", function(genome, aligner, ...) standardGeneric("createIndexPackage"))

## From fasta file
#' Description of Function
#'
#' @param xx xxx
#' @return xxxx

setMethod("createIndexPackage", c("list", "Aligner"), function(genome, aligner, sourcePackageFilepath=tempdir())
      {
          .progressReport("Creating index package")
          fastaFilepath <- file.path(genome$dir, genome$files)
          seedList <- .createSeedList(genome, aligner)
          sourcePackageFilepath <- .createIndexPackage(aligner,
                                                       fastaFilepath,
                                                       seedList,
                                                       sourcePackageFilepath)$pkgdir
          return(sourcePackageFilepath)
      })

### From BSgenome
setMethod("createIndexPackage", c("BSgenome", "Aligner"), function(genome, aligner, sourcePackageFilepath=tempdir())
      {
          .progressReport("Creating index package")
          dir.create(sourcePackageFilepath, showWarnings=FALSE, recursive=TRUE)
          fastaFilepath <- .BSgenomeSeqToFasta(genome)
          seedList <- .createSeedList(genome, aligner)
          sourcePackageFilepath <- .createIndexPackage(aligner,
                                                       fastaFilepath,
                                                       seedList,
                                                       sourcePackageFilepath)
          return(sourcePackageFilepath$pkgdir)
      })

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
        chrSeq <- DNAStringSet(unmasked(bsGenome[[chrT]]))
        names(chrSeq) <- chrT
        ## FH: I wonder why filepath can't be a connection. Would have made this appending nightmare a lot cleaner...
        write.XStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
        append <- TRUE
    }
    return(outFile)
}

## Same here, unless you plan to re-use this I'd prefer to put it in the method body. One less function in the call stack to find back in the source code...
.createIndexPackage <- function(aligner, fastaFilepath, seedList, sourcePackageFilepath=".")
{
    .requirePkg(aligner@pkgname)
    require(Biobase, quietly=TRUE)
    ## create package
    pkgname <- seedList$ALIGNERINDEXNAME
    destinationDir <- sourcePackageFilepath
    templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
    pkgDir <- createPackage(pkgname, destinationDir, templatePath, seedList, quiet=TRUE) 
    ## create index files
    indexDir <- file.path(pkgDir, "inst", "alignerIndex", seedList$GENOMENAME)
    .index(aligner, fastaFilepath, indexDir)
    return(pkgDir)
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
                 ALIGNER=aligner@pkgname,
                 ALIGNERVERSION=aligner@pkgversion,
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
