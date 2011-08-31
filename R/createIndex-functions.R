### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to create the index
###

.createGenomeIndex <- function(qProject)
{
    if(is.null(qProject@indexLocation))
        destDir <- qProject@genome$dir
    else
        destDir <- qProject@indexLocation
    .progressReport("Creating index")
    fastaFiles <- file.path(qProject@genome$dir, qProject@genome$files)
    index <- .createIndex(fastaFiles, qProject@aligner, qProject@genome$name, destDir)
    return(index)
}

.createAuxiliaryIndex <- function(qProject)
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
    index <- list(name=name,
                  path=indexName,
                  aligner=aligner$pkgname,
                  alignerversion=aligner$pkgversion,
                  #organism=name,
                  sourceurl=paste(fastaFiles, collapse=",")
                  )
    write.table(index, file=file.path(destDir, sprintf("%sIndex", aligner$pkgname), "index.tab"), 
                sep="\t", col.names=TRUE, row.names=FALSE)
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All Function to create an index package of a genome for a specific aligner
###

.createIndexPackage <- function(qProject, sourcePackageFilepath=tempdir(), lib.loc=NULL)
{
    .progressReport("Load genome")
    .requirePkg(qProject@aligner$pkgname, lib.loc=NULL)
    suppressPackageStartupMessages(require(Biobase, quietly=TRUE, lib.loc=lib.loc))
    genome <- .loadBSgenome(qProject@genome$name, lib.loc=lib.loc)
    .progressReport("Creating index package")
    fastaFilepath <- .BSgenomeSeqToFasta(genome)
    on.exit(unlink(fastaFilepath))
    seedList <- .createSeedList(genome, qProject@aligner)
    ## create package
    pkgname <- seedList$ALIGNERINDEXNAME
    dir.create(sourcePackageFilepath, showWarnings=FALSE, recursive=TRUE)
    destinationDir <- sourcePackageFilepath
    templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
    pkgDir <- Biobase::createPackage(pkgname, destinationDir, templatePath, seedList, quiet=TRUE)
    ## create index files
    indexDir <- file.path(pkgDir, "inst", "alignerIndex", seedList$GENOMENAME)
    .index(qProject@aligner, fastaFilepath, indexDir)
    return(pkgDir$pkgdir)
}

.BSgenomeSeqToFasta <- function(bsGenome, outFile=tempfile(fileext=".fa"))
{
    if(!is(bsGenome, "BSgenome"))
        stop("The variable bsGenome is not a BSgenome")
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
        ## seed$GENOMEVERSION <- installed.packages()[seed$GENOMENAME,'Version']
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
    seed$PKGTITLE <- paste(seed$ALIGNER, "Index of", seed$GENOMENAME)
    seed$PKGDESCRIPTION <- seed$PKGTITLE
    seed$ALIGNERINDEXNAME <- paste(seed$ALIGNER, gsub("_", "", seed$GENOMENAME), sep=".")
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
