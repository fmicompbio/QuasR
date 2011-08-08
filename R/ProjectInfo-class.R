setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))

setClass("qProject",
         representation(id="character",
                        name="character",
                        samples="data.frameOrNULL",
                        alignments="data.frameOrNULL",
                        genome="listOrNULL",
                        aligner="listOrNULL",
                        index="listOrNULL",
                        annotations="data.frameOrNULL",
                        path="character",
                        paired="logical",
                        junction="logical",
                        stranded="logical"
                        ),
         prototype(samples=NULL,
                   alignments=NULL,
                   genome=NULL,
                   aligner=NULL,
                   index=NULL,
                   annotations=NULL,
                   paired=FALSE,
                   junction=FALSE,
                   stranded=FALSE
                   )
)

setMethod("initialize", "qProject", function(.Object, name, path, ...){
    if(missing(name))
        name <- "qProject"
    id <- paste(name, format(Sys.time(), "%Y%m%d_%H%M%S"), sep="_")
    if(missing(path) || path == "."){
        path <- getwd()
    }
    callNextMethod(.Object, name=name, id=id, path=path, ...)
})

qProject <- function(sampleFile="Sample.txt", genome=".", annotationFile="Annotation.txt", aligner="Rbowtie", projectName="qProject", path=".", bisulfiteCoversion=FALSE, lib.loc=NULL, ...)
{
    .progressReport("Gathering file path information", phase=-1)
    samples <- readSamples(sampleFile)
    annotations <- readAnnotations(annotationFile)
    genome <- checkGenome(genome, ...)
    aligner <- loadAligner(aligner, bisulfiteCoversion, ...)
    alignments <- as.data.frame(matrix(0,
                                       nrow=nrow(samples),
                                       ncol=0,
                                       dimnames=dimnames(samples)[1]))
    project <- new("qProject", projectName, path, samples=samples, annotations=annotations, genome=genome, aligner=aligner, alignments=alignments)
    .progressReport(sprintf("Successfully created project '%s'", projectName), phase=1)
    return(project)
}

setMethod("show","qProject", function(object){
    cat("QuasRProject\n")
    cat("Project: " , object@name, "\n", sep="")
    cat("Genome:  ", object@genome$name, " is BSgenome=", object@genome$bsgenome, "\n", sep="")
    cat("Aligner: ", object@aligner$pkgname, " Version ", object@aligner$pkgversion, "\n", sep="")
    cat("Samples:\n", paste(object@samples$name, object@samples$filepath, sep="\t", collapse="\n"), "\n", sep="")
    cat("Annotations:\n", paste(object@annotations$feature, object@annotations$filepath, sep="\t", collapse="\n"), "\n", sep="")
    ## TODO write alignments in a nice way
    cat("Alignments:\n")   
    lapply(names(object@alignments), 
           function(index){
               cat(index, paste(as.matrix(object@alignments[index]), collapse="\n"), sep="\n")
               })
})

saveqProject <- function(project, filename)
{
    if(missing(filename))
        filename <- file.path(project@path, paste(project@id, "rds", sep="."))
    ## TODO calculate checksum of files or save modification date
    saveRDS(project, file=filename)
    return(filename)
}

readqProject <- function(filename)
{
    project <- readRDS(file=filename)
    ## TODO check qProject with checksum, path ...
    return(project)
}

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

loadBSgenome <- function(pkgname, ...)
{
    if(!pkgname %in% installed.packages()[,'Package']){
        ## download BSgenome
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(pkgname)
    }
    ## load BSgenome
    require(pkgname, character.only=TRUE, quietly=TRUE, ...)
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

checkGenome <- function(genomeName, ...)
{
    .progressReport("Check genome name")
    ## check if fasta file or directory
    if(file.exists(genomeName))
        return(loadFastaGenome(genomeName))
    ## check for installed BSgenome
    if(genomeName %in% installed.packages()[,'Package'])
        return(list(name=genomeName, bsgenome=TRUE))
    ## check if there is a BSgenome available with this name
    require(BSgenome, quietly=TRUE, ...)
    if(genomeName %in% available.genomes()){
        warning("Genome '", genomeName, "' is not installed. It will be downloaded and installed during the alignment process.")
        return(list(name=genomeName, bsgenome=TRUE))
    } else
        stop("Genome '", genomeName, "' not found.\nChoose a 'fasta' file, a directory containing 'fasta' files or one of the following BSgenomes:\n", paste(available.genomes(), "\n", collapse=""), sep="")
}
