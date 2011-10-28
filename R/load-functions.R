.readSamples <- function(file="samples.txt", paired=FALSE, sep="\t", row.names=NULL,  quote="\"", ...)
{
    .progressReport("Read sample file")
    tab <- read.table(file, header=TRUE, as.is=TRUE, sep=sep, quote=quote, fill=TRUE, ...)
    if(!paired){
        if(dirname(file) != ".")
            tab$FileName[dirname(tab$FileName) == "."] <- file.path(dirname(file), tab$FileName)
        checkFile <- file.exists(tab$FileName)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileName[!checkFile], collapse=", "))
        return(data.frame(name=tab$SampleName, filepath=I(tab$FileName), filetype=.fileType(tab$FileName)))
    }else{
        if(dirname(file) != "."){
            tab$FileName[dirname(tab$FileName) == "."] <- file.path(dirname(file), tab$FileName)
            tab$FileNameMate[dirname(tab$FileNameMate) == "."] <- file.path(dirname(file), tab$FileNameMate)
        }
        checkFile <- file.exists(tab$FileName)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileName[!checkFile], collapse=", "))
        checkFile <- file.exists(tab$FileNameMate)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileNameMate[!checkFile], collapse=", "))
        return(data.frame(name=tab$SampleName, filepath1=I(tab$FileName),
                          filepath2=I(tab$FileNameMate), filetype=.fileType(tab$FileName)))
    }
}

.loadAlignments <- function(qproject){
    ## emtpy alignment data.frame
    alignments <- as.data.frame(matrix(NA,
                                       nrow=nrow(qproject@samples),
                                       ncol=0,
                                       dimnames=dimnames(qproject@samples)[1]))
    ## load genomic alignment
    alignments$genome <- unlist(lapply(.baseFileName(qproject@samples$filepath), 
           function(f){
               lst <- list.files(qproject@path, pattern=sprintf("%s.*\\.bam$", f), full.names=TRUE)
               ifelse(length(lst),
                      .getBamFile(lst, qproject@genome$name, qproject@alignmentParameter),
                      NA)
           }))
           #TODO what when more than one results
    ## load auxiliary alignment
    aux <- qproject@annotations$filepath[qproject@annotations$filetype == "fasta"]
    names(aux) <- qproject@annotations$feature[qproject@annotations$filetype == "fasta"]
    auxAlignment <- lapply(aux, function(a){
        unlist(lapply(.baseFileName(qproject@samples$filepath), 
                  function(f){
                      lst <- list.files(qproject@path, pattern=sprintf("%s.*\\.bam$", f), full.names=TRUE)
                      ifelse(length(lst),
                             .getBamFile(lst, a, qproject@alignmentParameter),
                             NA)
                        }))
     })
     #names(auxAlignment) <- name(aux)
     alignments <- cbind.data.frame(alignments, 
                                    as.data.frame(auxAlignment, stringsAsFactors=FALSE), 
                                    stringsAsFactors=FALSE)
    ## set external bamfiles
    if(any(idx <- qproject@samples$filetype == "bam")){
        alignments[idx,] <- ""
        alignments$genome[idx] <- qproject@samples$filepath[idx]
    }
    return(alignments)
}

.readAnnotations <- function(file="annotations.txt", sep="\t", row.names=NULL, quote="\"", ...)
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
    return(data.frame(feature=tab$Class, filepath=I(tab$FileName), filetype=.fileType(tab$FileName)))
}

.loadBSgenome <- function(pkgname, lib.loc=NULL)
{
    if(!pkgname %in% installed.packages()[,'Package']){
        ## download BSgenome
        if(suppressWarnings(require(BiocInstaller, quietly=TRUE, lib.loc=lib.loc)))
            BiocInstaller::biocLite(pkgname, suppressUpdates=TRUE, lib=lib.loc)
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

.loadFastaGenome <- function(dirname)
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
    name <- gsub("_", "", dirname) # TODO maybe remove gsub
    shortname <- gsub("_", "", basename(name))
    return(list(name=name, shortname=shortname, dir=dir, files=files, bsgenome=FALSE))
}

.checkGenome <- function(genomeName, lib.loc=NULL)
{
    .progressReport("Check genome name")
    ## check if fasta file or directory
    if(file.exists(genomeName)){
        genomeName <- tools::file_path_as_absolute(genomeName)
        return(.loadFastaGenome(genomeName))
    }
    ## check for installed BSgenome
    if(genomeName %in% installed.packages()[,'Package'])
        return(list(name=genomeName, bsgenome=TRUE))
    ## check if there is a BSgenome available with this name
    require(BSgenome, quietly=TRUE, lib.loc=lib.loc)
    if(genomeName %in% available.genomes()){
        return(list(name=genomeName, bsgenome=TRUE))
    } else
        stop("Genome '", genomeName, "' not found.\nChoose a 'fasta' file, a directory containing 'fasta' files or one of the following BSgenomes:\n", paste(available.genomes(), "\n", collapse=""), sep="")
}

.loadAligner <- function(pkgname, lib.loc=NULL){
    .progressReport("Loading aligner")
    .requirePkg(pkgname, lib.loc=lib.loc)
    aligner <- list(pkgname=pkgname,
                    pkgversion=installed.packages()[pkgname, 'Version']
                    )
    return(aligner)
}

.loadIndex <- function(qProject, lib=NULL, lib.loc=NULL)
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
            srcPkgDir <- .createIndexPackage(qProject, lib.loc=lib.loc)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)
            index <- get(ls(sprintf("package:%s", pkgname))[1])
        }
    } else {
        if(is.null(qProject@indexLocation))
            indexDir <- file.path(qProject@genome$dir, sprintf("%sIndex", qProject@aligner$pkgname))
        else
            indexDir <- file.path(qProject@indexLocation, sprintf("%sIndex", qProject@aligner$pkgname))
        if(!file.exists(indexDir)){
            .progressReport("No index found. Create index now")
            index <- .createGenomeIndex(qProject)
        } else {
            index <- read.table(file=file.path(indexDir,"index.tab"), sep="\t", header=TRUE)
            ## set index path
            index$path <- file.path(indexDir, index$shortname) # TODO or from index.tab file
        }
    }
    return(index)
}
