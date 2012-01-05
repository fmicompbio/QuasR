### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to create the index
###

.createGenomeIndex <- function(qproject)
{
    if(is.null(qproject@env$indexLocation))
        destDir <- tools::file_path_as_absolute(qproject@env$genome$dir)
    else
        destDir <- qproject@env$indexLocation
    .progressReport("Creating index")
    fastaFiles <- file.path(qproject@env$genome$dir, qproject@env$genome$files)
    index <- .createIndex(fastaFiles, qproject@env$aligner, qproject@env$genome$name, qproject@env$genome$shortname, destDir)
    return(index)
}

.createAuxiliaryIndex <- function(qproject)
{
    .progressReport("Creating auxiliary index")
    isFastaFormat <- qproject@env$auxiliaries$filetype == "fasta"
    indexes <- mapply(.createIndex,
                      fastaFiles=qproject@env$auxiliaries[isFastaFormat,]$filepath,
                      name=qproject@env$auxiliaries[isFastaFormat,]$filepath,
                      shortname=as.character(qproject@env$auxiliaries[isFastaFormat,]$feature),
                      MoreArgs=list(aligner=qproject@env$aligner, 
                                    destDir=tempfile(pattern="auxIndex_", tmpdir=qproject@env$cacheDir)),
                      SIMPLIFY=FALSE,
                      USE.NAMES = TRUE)
#     names(indexes) <- .baseFileName(names(indexes)) ## FIXME: maybe use name or name should be basefilename
    names(indexes) <- as.character(qproject@env$auxiliaries[isFastaFormat,]$feature)
    return(indexes)
}

.createIndex <- function(fastaFiles, aligner, name, shortname, destDir=tempfile())
{
    indexName <- file.path(destDir, sprintf("%sIndex_%s", aligner$pkgname, shortname), shortname)
    if(!dir.create(dirname(indexName), showWarnings=TRUE, recursive=TRUE))
        stop("Could not create index directory '", dirname(indexName),"'.")
    .index(aligner, fastaFiles, indexName)
    index <- list(name=name,
                  shortname=shortname,
                  path=indexName,
                  aligner=aligner$pkgname,
                  alignerversion=aligner$pkgversion,
                  #organism=name,
                  sourceurl=paste(fastaFiles, collapse=","),
                  md5sum=paste(tools::md5sum(fastaFiles), collapse=",")
                  )
    write.table(index, file=sprintf("%s.tab", indexName), 
                sep="\t", col.names=TRUE, row.names=FALSE)
    return(index)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### All Function to create an index package of a genome for a specific aligner
###

.createIndexPackage <- function(qproject, lib.loc=NULL)
{
    .progressReport("Load genome")
    .requirePkg(qproject@env$aligner$pkgname, lib.loc=NULL)
#     suppressPackageStartupMessages(require(Biobase, quietly=TRUE, lib.loc=lib.loc))
    genome <- .loadBSgenome(qproject@env$genome$name, lib.loc=lib.loc)
    .progressReport("Creating index package")
    fastaFilepath <- .BSgenomeSeqToFasta(genome, tempfile(tmpdir=qproject@env$cacheDir, fileext=".fa"))
    on.exit(unlink(fastaFilepath))
    seedList <- .createSeedList(genome, qproject@env$aligner)
    seedList$MD5SUM <- tools::md5sum(fastaFilepath)
    ## create package
    pkgname <- seedList$ALIGNERINDEXNAME
    destinationDir <- qproject@env$cacheDir  
    dir.create(destinationDir, showWarnings=FALSE, recursive=TRUE)
    templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
    pkgDir <- Biobase::createPackage(pkgname, destinationDir, templatePath, seedList, quiet=TRUE)
    ## create index files
    indexDir <- file.path(pkgDir, "inst", "alignerIndex", seedList$GENOMENAME)
    .index(qproject@env$aligner, fastaFilepath, indexDir)
    return(pkgDir$pkgdir)
}

.BSgenomeSeqToFasta <- function(bsgenome, outFile=tempfile(fileext=".fa"))
{
    if(!is(bsgenome, "BSgenome"))
        stop("The variable 'bsgenome'' is not a BSgenome")
    append <- FALSE
    for(chrT in seqnames(bsgenome)){
        if(is.null(masks(bsgenome[[chrT]])))
            chrSeq <- DNAStringSet(bsgenome[[chrT]])
        else
            chrSeq <- DNAStringSet(unmasked(bsgenome[[chrT]]))
        names(chrSeq) <- chrT
        write.XStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
        append <- TRUE
    }
    return(outFile)
}

.installIndexPackage <- function(pkgdir, lib=NULL)
{
    .progressReport("Build index packages")
    curwd <- setwd(dirname(pkgdir)) # crash if wd not changed 
    on.exit(setwd(curwd))
    buildFun <- tools:::.build_packages
    fbody <- body(buildFun)
    body(buildFun) <- fbody[1:(length(fbody)-1)]
    out <- capture.output(buildFun(basename(pkgdir)))
    .progressReport("Installing index packages")
    srcpkg <- file.path(getwd(), dir(pattern="\\.tar.gz$"))
    if(is.null(lib))
        lib <- .libPaths()[1]
    out <- c(out, capture.output(install.packages(srcpkg, repos=NULL, dependencies=FALSE, lib=lib, type="source")))
    unlink(srcpkg, recursive=TRUE)
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
