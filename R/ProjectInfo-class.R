setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
#setClassUnion("AlignerOrNULL", c("Aligner", "NULL"))
setClassUnion("listOrNULL", c("list", "BSgenome", "NULL"))

setClass("ProjectInfo", 
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

setMethod("initialize", "ProjectInfo", function(.Object, name, path, ...)
      {
          if(missing(name))
              name <- "projectinfo"
          id <- paste(name, format(Sys.time(), "%Y%m%d_%H%M%S"), sep="_")
          if(missing(path) || path == "."){
              path <- getwd()
          }
          callNextMethod(.Object, name=name, id=id, path=path, ...)
      })

ProjectInfo <- function(sampleFile="Sample.txt", genome=".", annotationFile="Annotation.txt", aligner="Rbowtie", projectName="projectinfo", path=".", bisulfiteCoversion=FALSE, lib=NULL, ...)
{
    .progressReport("Gathering file path information", phase=-1)
    samples <- readSamples(sampleFile)
    annotations <- readAnnotations(annotationFile)
    genome <- checkGenome(genome)
    aligner <- loadAligner(aligner, bisulfiteCoversion)
    alignments <- as.data.frame(matrix(0,
                                       nrow=nrow(samples),
                                       ncol=0,
                                       dimnames=dimnames(samples)[1]))                  
    project <- new("ProjectInfo", projectName, path, samples=samples, annotations=annotations, genome=genome, aligner=aligner, alignments=alignments)
    .progressReport(sprintf("Successfully created project '%s'", projectName), phase=1)
    return(project)
}

##setMethod("show","ProjectInfo",
##  function(object){
##    cat("Project:" ,object@name, "(id:", object@id, ")\n")
##    cat("Genome:\n", paste(object@genome, collapse="-"), "\n")
##    cat("Samples:\n", paste(object@samples, collapse="\n"), "\n")
##    cat("Reference:\n", paste(object@references, collapse="\n"), "\n")
##  }
##)

saveProjectInfo <- function(project, filename)
{
    if(missing(filename))
        filename <- file.path(project@path, paste(project@id, "rds", sep="."))
    ## TODO calculate checksum of files or save modification date
    saveRDS(project, file=filename)
    return(filename)
}

readProjectInfo <- function(filename)
{
    project <- readRDS(file=filename)
    ## TODO check projectInfo with checksum, path ...
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

loadBSgenome <- function(pkgname)
{
    if(!pkgname %in% installed.packages()[,'Package']){
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(pkgname)
    }
    require(pkgname, character.only=TRUE, quietly=TRUE)
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

checkGenome <- function(genomeName)
{
    .progressReport("Check genome name")
    ## check if fasta file or directory
    if(file.exists(genomeName))
        return(loadFastaGenome(genomeName))
    ## check for installed BSgenome    
    if(genomeName %in% installed.packages()[,'Package'])
        return(list(name=genomeName, bsgenome=TRUE))
    ## check if there is a BSgenome available with this name
    require(BSgenome, quietly=TRUE)
    if(genomeName %in% available.genomes()){
        warning("Genome '", genomeName, "' is not installed. It will be downloaded and installed during the alignment process.")
        return(list(name=genomeName, bsgenome=TRUE))
    } else
        stop("Genome '", genomeName, "' not found.\nChoose a 'fasta' file, a directory containing 'fasta' files or one of the following BSgenomes:\n", paste(available.genomes(), "\n", collapse="\n"), sep="")
}
