.readSamples <- function(file="samples.txt", paired=FALSE, sep="\t", row.names=NULL,  quote="\"", ...)
{
    .progressReport("Read sample file")
    tab <- read.table(file, header=TRUE, as.is=TRUE, sep=sep, quote=quote, fill=TRUE, ...)
    if(!paired){
        ## convert to absolut filename
        tab$FileName[dirname(tab$FileName) == "."] <- 
            file.path(dirname(file), tab$FileName[dirname(tab$FileName) == "."])
        tab$FileName[substr(tab$FileName, 1, 1) != "/"]  <- 
            file.path(dirname(file), tab$FileName[substr(tab$FileName, 1, 1) != "/"])
        ## check if files exits
        checkFile <- file.exists(tab$FileName)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileName[!checkFile], collapse=", "))
        ## check if filename is duplicated
        if(any(checkFile <- duplicated(basename(tab$FileName))))
            stop("Duplicated filename ", paste(unique(basename(tab$FileName[checkFile])), collapse=", "))
        return(data.frame(name=tab$SampleName, filepath=I(tab$FileName), filetype=.fileType(tab$FileName)))
    }else{
        ## convert to absolut filename
        tab$FileName1[dirname(tab$FileName1) == "."] <-
            file.path(dirname(file), tab$FileName1[dirname(tab$FileName1) == "."])        
        tab$FileName1[substr(tab$FileName1, 1, 1) != "/"]  <-
            file.path(dirname(file), tab$FileName1[substr(tab$FileName1, 1, 1) != "/"])
        tab$FileName2[dirname(tab$FileName2) == "."] <-
            file.path(dirname(file), tab$FileName2[dirname(tab$FileName2) == "."])        
        tab$FileName2[substr(tab$FileName2, 1, 1) != "/"]  <-
            file.path(dirname(file), tab$FileName2[substr(tab$FileName2, 1, 1) != "/"])
        ## check if files exits
        checkFile <- file.exists(tab$FileName1)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileName1[!checkFile], collapse=", "))
        checkFile <- file.exists(tab$FileName2)
        if(any(!checkFile))
            stop("File not found: ", paste(tab$FileName2[!checkFile], collapse=", "))
        ## check if filename is duplicated
        if(any(checkFile <- duplicated(basename(tab$FileName1))))
            stop("Duplicated filename ", paste(unique(basename(tab$FileName1[checkFile])), collapse=", ")) #TODO for mate too
        return(data.frame(name=tab$SampleName, filepath1=I(tab$FileName1),
                          filepath2=I(tab$FileName2), filetype=.fileType(tab$FileName1)))
    }
}

.loadAlignments <- function(qproject){
    .progressReport("Read annotation file")
    ## emtpy alignment data.frame
    alignments <- as.data.frame(matrix(NA,
                                       nrow=nrow(qproject@env$samples),
                                       ncol=0,
                                       dimnames=dimnames(qproject@env$samples)[1]))
    ## load genomic alignment
    isBamfile <- qproject@env$samples$filetype == "bam"
    alignments$genome[!isBamfile] <- unlist(lapply(qproject@env$samples$filepath[!isBamfile], 
           function(filepath){
               if(is.null(qproject@env$bamfileDir))
                   lst <- list.files(dirname(filepath), 
                                     pattern=sprintf("%s.*\\.bam$", .baseFileName(filepath)), full.names=TRUE) #TODO exacter regex
               else
                   lst <- list.files(qproject@env$bamfileDir, 
                                     pattern=sprintf("%s.*\\.bam$", .baseFileName(filepath)), full.names=TRUE)
               ifelse(length(lst),
                      .getBamFile(lst, qproject@env$genome$name, qproject@env$alignmentParameter),
                      NA)
           }))
           #TODO what when more than one results
    ## load auxiliary alignment
    if(!is.null(qproject@env$annotations)){
        aux <- qproject@env$annotations$filepath[qproject@env$annotations$filetype == "fasta"]
        names(aux) <- qproject@env$annotations$feature[qproject@env$annotations$filetype == "fasta"]
        auxAlignment <- lapply(aux, function(a){
            unlist(lapply(qproject@env$samples$filepath, 
                      function(filepath){
                              if(is.null(qproject@env$bamfileDir))
                                  lst <- list.files(dirname(filepath), 
                                                    pattern=sprintf("%s.*\\.bam$", .baseFileName(filepath)), 
                                                    full.names=TRUE) #TODO exacter regex
                              else
                                  lst <- list.files(qproject@env$bamfileDir, 
                                                    pattern=sprintf("%s.*\\.bam$", .baseFileName(filepath)), 
                                                    full.names=TRUE)
                              ifelse(length(lst),
                                     .getBamFile(lst, a, qproject@env$alignmentParameter),
                                     NA)
                          }))
        })
        #names(auxAlignment) <- name(aux)
        alignments <- cbind.data.frame(alignments, 
                                    as.data.frame(auxAlignment, stringsAsFactors=FALSE), 
                                    stringsAsFactors=FALSE)
    }

    ## set external bamfiles
    if(any(idx <- qproject@env$samples$filetype == "bam")){
        alignments[idx,] <- ""
        alignments$genome[idx] <- qproject@env$samples$filepath[idx]
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

.loadIndex <- function(qproject, lib=NULL, lib.loc=NULL)
{
    .progressReport("Loading genome index")
    ## TODO check if genome and aligner is not NULL
    if(qproject@env$genome$bsgenome)
    {
        #Genome is BSgenome
        pkgname <- paste(qproject@env$aligner$pkgname,
                         qproject@env$genome$name,
                         sep=".")
        if(suppressWarnings(require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)))
        {
            index <- eval(parse(text=ls(sprintf("package:%s", pkgname))))
        } else {
            .progressReport("No index found. Create index now")
            ## create and install index
            srcPkgDir <- .createIndexPackage(qproject, lib.loc=lib.loc)
            .installIndexPackage(srcPkgDir, lib=lib)
            pkgname <- basename(srcPkgDir)
            require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)
            index <- get(ls(sprintf("package:%s", pkgname))[1])
        }
    } else {
        if(is.null(qproject@env$indexLocation))
            indexDir <- file.path(qproject@env$genome$dir, sprintf("%sIndex", qproject@env$aligner$pkgname))
        else
            indexDir <- file.path(qproject@env$indexLocation, sprintf("%sIndex", qproject@env$aligner$pkgname))
        if(!file.exists(indexDir)){
            .progressReport("No index found. Create index now")
            index <- .createGenomeIndex(qproject)
        } else {
            index <- read.table(file=file.path(indexDir,"index.tab"), sep="\t", header=TRUE)
            ## set index path
            index$path <- file.path(indexDir, index$shortname) # TODO or from index.tab file
        }
    }
    return(index)
}
