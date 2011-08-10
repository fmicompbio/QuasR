readSamples <- function(file="samples.txt", sep="\t", row.names=NULL,  quote="\"", ...)
{
    .progressReport("Read sample file")
    tab <- read.table(file, header=TRUE, as.is=TRUE, sep=sep, quote=quote, fill=TRUE, ...)
    if(dirname(file) != ".")
        tab$FileName[dirname(tab$FileName) == "."] <- file.path(dirname(file), tab$FileName)
    checkFile <- file.exists(tab$FileName)
    if(any(!checkFile))
        stop("File not found: ", paste(tab$FileName[!checkFile], collapse=", "))
    return(data.frame(name=tab$SampleName, filepath=I(tab$FileName)))
}

readAnnotations <- function(file="annotations.txt", sep="\t", row.names=NULL, quote="\"", ...)
{
    .progressReport("Read annotation file")
    if(!file.exists(file))
        stop("File '", file, "'not found.")
    tab <- read.table(file, header=TRUE, as.is=TRUE, sep=sep, quote=quote, fill=TRUE, ...)
    if(dirname(file) != ".")
        tab$FileName[dirname(tab$FileName) == "."] <- file.path(dirname(file), tab$FileName)
    checkFile <- file.exists(tab$FileName)
    if(any(!checkFile))
        stop("File not found: ", paste(tab$FileName[!checkFile], collapse=", "))
    return(data.frame(feature=tab$Feature, filepath=I(tab$FileName)))
}

loadBSgenome <- function(pkgname, lib.loc=NULL)
{
    if(!pkgname %in% installed.packages()[,'Package']){
        ## download BSgenome
        if(require(BiocInstaller, lib.loc=lib.loc)){ 
            ##source("http://bioconductor.org/scratch-repos/biocLite.R")
            BiocInstaller::biocLite(pkgname)
        } else {
            source("http://www.bioconductor.org/biocLite.R")       
            biocLite(pkgname)
        }
    }
    ## load BSgenome
    require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)
    ## get BSgenome object
    ## genome <- eval(parse(text=strsplit(pkgname,"\\.")[[1]][2]))
    genome <- eval(parse(text=ls(sprintf("package:%s", pkgname))))
    if(!is(genome,"BSgenome"))
        stop("'", pkgname, "' is not a BSgenome package")
    return(genome)
}

loadFastaGenome <- function(dirname)
{
    if(file.info(dirname)$isdir){
        ##dirname is a directory
        dir <- dirname
        files <- list.files(dirname, pattern="\\.fa$|\\.fna$|\\.fasta$")
    }else{
        ##dirname is a file
        dir <- dirname(dirname)
        files <- basename(dirname)
      }
    name <- .baseFileName(gsub("_", "",dirname))
    return(list(name=name, dir=dir, files=files, bsgenome=FALSE))
}

checkGenome <- function(genomeName, lib.loc=NULL)
{
    .progressReport("Check genome name")
    ## check if fasta file or directory
    if(file.exists(genomeName))
        return(loadFastaGenome(genomeName))
    ## check for installed BSgenome
    if(genomeName %in% installed.packages()[,'Package'])
        return(list(name=genomeName, bsgenome=TRUE))
    ## check if there is a BSgenome available with this name
    require(BSgenome, quietly=TRUE, lib.loc=lib.loc)
    if(genomeName %in% available.genomes()){
        warning("Genome '", genomeName, "' is not installed. It will be downloaded and installed during the alignment process.")
        return(list(name=genomeName, bsgenome=TRUE))
    } else
        stop("Genome '", genomeName, "' not found.\nChoose a 'fasta' file, a directory containing 'fasta' files or one of the following BSgenomes:\n", paste(available.genomes(), "\n", collapse=""), sep="")
}

loadAligner <- function(pkgname, lib.loc=NULL){
    .progressReport("Loading aligner")
    .requirePkg(pkgname, lib.loc=lib.loc)
    aligner <- list(pkgname=pkgname,
                    pkgversion=installed.packages()[pkgname, 'Version']
                    )
    return(aligner)
}

loadIndex <- function(qProject, lib=NULL, lib.loc=NULL)
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
            srcPkgDir <- createIndexPackage(qProject, lib.loc=lib.loc)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)
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
