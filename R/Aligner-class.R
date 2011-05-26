# Aligner Class
setClass("Aligner", 
    representation(
        pkgname="character",
        pkgversion="character"
    )
)

setMethod("initialize", "Aligner", 
    function(.Object, pkgname) {
        .requirePkg(pkgname)
        pkgname <- pkgname
        pkgversion <- installed.packages()[pkgname, 'Version']
        callNextMethod(.Object, pkgname=pkgname, pkgversion=pkgversion)
    }
)

Aligner <- function(aligner){
  new("Aligner", aligner)
}

setMethod("show","Aligner",
    function(object){
        cat(object@pkgname, "Aligner Version", object@pkgversion, "\n")
    }
)

setGeneric("execute", function(object, callstr){
         standardGeneric("execute")
})

setMethod("execute", "Aligner",
function(object, callstr){
  .requirePkg(object@pkgname) 
  fun <- sprintf("%s:::.execute('%s')", object@pkgname, callstr)
  return(eval(parse(text=fun)))
})

#index <- function(projectInfo){
#  .index(projectInfo@aligner, projectInfo@genome$path, projectInfo@genome$path)
#  project
#}

index <- function(aligner, genome, indexDir){
  .index(aligner, genome, indexDir)
  project
}

.index <- function(aligner, fastaFilepath, indexDir){
  cat("Create index in", indexDir, "\n")
  if(aligner@pkgname == "Rbowtie")
    out <- QuasR::execute(aligner,
                   paste("bowtie-build",
                         paste(fastaFilepath, collapse=","),
                         indexDir))
  else if(aligner@pkgname == "Rbwa")
    out <- QuasR::execute(aligner,
                   paste("bwa index -p",
                         indexDir,
                         fastaFilepath))
  else
    stop("The '", aligner@pkgname, "' Aligner is not supported.")
}

.align <- function(readsFilepath, aligner, index, outpath){
  outputFilename <- file.path(outpath,
                              sprintf("%s-%s.sam",
                                      .baseFileName(readsFilepath),
                                      index$name))
  if(aligner@pkgname == "Rbowtie")
    QuasR::execute(aligner, paste("bowtie", index$path, "-f", readsFilepath, "-S", outputFilename))
  else if(aligner@pkgname == "Rbwa")
    QuasR::execute(aligner, paste("bwa bwasw", index$path, readsFilepath, ">", outputFilename))
  else
    stop("The '", aligner@pkgname, "' Aligner is not supported.")
  return(outputFilename)
}

align <- function(projectInfo){
  if(!is(projectInfo,"ProjectInfo"))
    stop("The object '", class(projectInfo), "' is not a 'ProjectInfo' object.")
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

#loadAligner <- function(pkgname){
#  aligner <- Aligner(pkgname)
#  return(aligner)
#}

loadIndex <- function(genome, aligner){
  if(is(genome,"BSgenome")){
    pkgname <- paste(aligner@pkgname, bsgenomeName(genome), sep=".")
    if(length(grep(pkgname, installed.packages()[,'Package'])) >= 1L){
      require(pkgname, character.only=TRUE)
      index <- eval(parse(text=ls(sprintf("package:%s", pkgname))))
    }else{
      srcPkgDir <- createIndexPackage(genome, aligner)
      QuasR:::.installIndexPackage(srcPkgDir)
      pkgname <- basename(srcPkgDir)
      require(pkgname, character.only=TRUE)
      index <- eval(parse(text=ls(sprintf("package:%s", pkgname)))) 
    }
  }else{
       #TODO
    indexDir <- file.path(genome$dir, sprintf("%sIndex", aligner@pkgname))
    if(!file.exists(indexDir))         
      indexDir <- createIndexPackage(genome, aligner)
    name <- paste(aligner@pkgname, gsub("_", "", genome$name), sep=".") 
    #alignerversion <- aligner@pkgversion
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

setGeneric("createIndexPackage", function(genome, aligner, ...){
         standardGeneric("createIndexPackage")
})

## From fasta file
#' Description of Function
#'
#' @param xx xxx
#' @return xxxx
setMethod("createIndexPackage", c("list", "Aligner"),
function(genome, aligner, sourcePackageFilepath="."){
  #fastaFilepath <- file.path(genome$dir, genome$files)
  #seedList <- .createSeedList(genome, aligner)
  #sourcePackageFilepath <- .createIndexPackage(aligner,
  #                                             fastaFilepath,
  #                                             seedList,
  #                                             sourcePackageFilepath)$pkgdir
  sourcePackageFilepath <- createIndex(genome, aligner)
  return(sourcePackageFilepath)
})

createIndex <- function(genome, aligner){
  fastaFilepath <- file.path(genome$dir, genome$files)
  indexName <-  sprintf("%sIndex", aligner@pkgname)
  indexDir <- file.path(genome$dir, indexName, genome$name)
  dir.create(dirname(indexDir))
  .index(aligner, fastaFilepath, indexDir)
  return(indexDir)
}

### From BSgenome
setMethod("createIndexPackage", c("BSgenome", "Aligner"),
function(genome, aligner, sourcePackageFilepath="."){
  fastaFilepath <- .BSgenomeSeqToFasta(genome)
  seedList <- .createSeedList(genome, aligner)
  sourcePackageFilepath <- .createIndexPackage(aligner,
                                               fastaFilepath,
                                               seedList,
                                               sourcePackageFilepath)
  return(sourcePackageFilepath$pkgdir)
})

.BSgenomeSeqToFasta <- function(bsGenome, outFile=NULL){
  if(!is(bsGenome, "BSgenome"))
    stop("The variable bsGenome is not a BSgenome")
  if(missing(outFile) || is.null(outFile)){
    td <- tempdir()
    outFile <- file.path(td, paste(bsgenomeName(bsGenome),"fa",sep="."))    
  }
  append=FALSE
  for(chrT in seqnames(bsGenome)){
    chrSeq <- DNAStringSet(unmasked(bsGenome[[chrT]]))
    names(chrSeq) <- chrT
    write.XStringSet(chrSeq, file=outFile, format="fasta", append=append)
    append=TRUE
  }
  return(outFile)
}

.createIndexPackage <- function(aligner, fastaFilepath, seedList, sourcePackageFilepath="."){
  .requirePkg(aligner@pkgname)
  require(Biobase)
  # create package
  pkgname <- seedList$ALIGNERINDEXNAME
  destinationDir <- sourcePackageFilepath
  templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
  pkgDir <- createPackage(pkgname, destinationDir, templatePath, seedList) 
  #create index files
  indexDir <- file.path(pkgDir, "inst", "index", seedList$GENOMENAME)
  index(aligner, fastaFilepath, indexDir)
  return(pkgDir)
}

.installIndexPackage <- function(pkgdir){
    #build <- paste("R CMD build", pkgdir)
    #check <- paste("R CMD check", pkgdir)
    install <- paste("R CMD INSTALL", pkgdir)
    output <- system(install, intern=TRUE)
    return(output)
}

.defaultSeedList <- function(genome, aligner){
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
               #genome seeds
               GENOMENAME="",
               PROVIDER="",
               PROVIDERVERSION="",
               RELEASEDATE="",
               RELEASENAME="",
               ORGANISM="",
               SPECIES="",
               SRCDATAFILES="",
               ORGANISMBIOCVIEW="",
               #aligner seeds
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
                            pkgversion, author, maintainer, license){
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
#  seed$PKGDETAILS <- if(!missing(pkgdetails)) pkgdetails
#  seed$PKGEXAMPLES <- if(!missing(pkgexamples)) pkgexamples
#  if(!missing(indexname))
#    seed$ALIGNERINDEXNAME <- indexname
  return(seed)
}
