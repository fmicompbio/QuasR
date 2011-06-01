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

##index <- function(projectInfo){
##  .index(projectInfo@aligner, projectInfo@genome$path, projectInfo@genome$path)
##  project
##}

index <- function(aligner, genome, indexDir){
    .index(aligner, genome, indexDir)
    ## FH: I don't think you want to return 'project' here. This might not throw an error but only because of
    ## R's lexical scoping...
    ## project
}


.index <- function(aligner, fastaFilepath, indexDir)
{
    dir.create(indexDir, recursive=TRUE, showWarnings=FALSE)
    if(aligner@pkgname == "Rbowtie")
    {
        out <- execute(aligner, sprintf("bowtie-build %s %s", paste(fastaFilepath, collapse=","), indexDir))
    } else if(aligner@pkgname == "Rbwa") {
        out <- execute(aligner,
                       paste("bwa index -p",
                             indexDir,
                             fastaFilepath))
    } else {
        stop("The '", aligner@pkgname, "' Aligner is not supported.")
    }
}

## Again, we the separate function?
.align <- function(readsFilepath, aligner, index, outpath)
{
    outputFilename <- file.path(outpath,
                                sprintf("%s-%s.sam",
                                        .baseFileName(readsFilepath),
                                        index$name))
    .progressReport(sprintf("\nAligning reads for sample '%s'", basename(readsFilepath)), phase=-1)
    if(aligner@pkgname == "Rbowtie")
        execute(aligner, paste("bowtie", index$path, "-f", readsFilepath, "-S", outputFilename))
    else if(aligner@pkgname == "Rbwa")
        execute(aligner, paste("bwa bwasw", index$path, readsFilepath, ">", outputFilename))
    else
        stop("The '", aligner@pkgname, "' Aligner is not supported.")
    return(outputFilename)
}

align <- function(projectInfo)
{
    if(!is(projectInfo, "ProjectInfo"))
        stop("The object '", class(projectInfo), "' is not a 'ProjectInfo' object.")
    .progressReport("Starting alignments", phase=-1)
    projectInfo@samples$rawAlignment <- lapply(projectInfo@samples$filepath,
                                               .align,
                                               projectInfo@aligner,
                                               projectInfo@index,
                                               projectInfo@path)
    return(projectInfo)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 
###

##loadAligner <- function(pkgname){
##  aligner <- Aligner(pkgname)
##  return(aligner)
##}

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
            srcPkgDir <- createIndexPackage(genome, aligner)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE)
            index <- get(ls(sprintf("package:%s", pkgname))[1]) 
        }
    } else {
        ##TODO
        indexDir <- file.path(genome$dir, sprintf("%sIndex", aligner@pkgname))
        if(!file.exists(indexDir))         
            indexDir <- createIndexPackage(genome, aligner)
        name <- paste(aligner@pkgname, gsub("_", "", genome$name), sep=".") 
        ##alignerversion <- aligner@pkgversion
        sourceurl <- file.path(genome$path, genome$files)
        index <- list(name=name,
                      path=indexDir,
                      aligner=aligner@pkgname,
                      organism=genome$name,
                      sourceurl=sourceurl
                      )
    }
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All Function to create the index of a genome for a specific aligner
###

setGeneric("createIndexPackage", function(genome, aligner, ...) standardGeneric("createIndexPackage"))

## From fasta file
#' Description of Function
#'
#' @param xx xxx
#' @return xxxx
setMethod("createIndexPackage", c("list", "Aligner"), function(genome, aligner, sourcePackageFilepath=tempfile())
      {
          ##fastaFilepath <- file.path(genome$dir, genome$files)
          ##seedList <- .createSeedList(genome, aligner)
          ##sourcePackageFilepath <- .createIndexPackage(aligner,
          ##                                             fastaFilepath,
          ##                                             seedList,
          ##                                             sourcePackageFilepath)$pkgdir
           .progressReport("No index found. Creating one now")
          dir.create(sourcePackageFilepath, showWarnings=FALSE, recursive=TRUE)
          sourcePackageFilepath <- createIndex(genome, aligner)
          return(sourcePackageFilepath)
      })


## FH: Why not put this in the previous method as the method body? Do you expect to re-use this?
createIndex <- function(genome, aligner)
{
    fastaFilepath <- file.path(genome$dir, genome$files)
    indexName <-  sprintf("%sIndex", aligner@pkgname)
    indexDir <- file.path(genome$dir, indexName, genome$name)
    dir.create(dirname(indexDir), showWarnings=FALSE, recursive=TRUE)
    .index(aligner, fastaFilepath, indexDir)
    return(indexDir)
}

### From BSgenome
setMethod("createIndexPackage", c("BSgenome", "Aligner"), function(genome, aligner, sourcePackageFilepath=tempfile())
      {
          .progressReport("No index found. Creating one now")
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
    ## FIXME: only the first chromsome while debugging
    ##for(chrT in seqnames(bsGenome)){
    for(chrT in seqnames(bsGenome)[1]){
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
    ## FH: this should be in the Depends field or we should import what we need
    require(Biobase, quietly=TRUE)
    ## create package
    pkgname <- seedList$ALIGNERINDEXNAME
    destinationDir <- sourcePackageFilepath
    templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
    pkgDir <- createPackage(pkgname, destinationDir, templatePath, seedList, quiet=TRUE) 
    ## create index files
    indexDir <- file.path(pkgDir, "inst", "alignerIndex", seedList$GENOMENAME)
    index(aligner, fastaFilepath, indexDir)
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
