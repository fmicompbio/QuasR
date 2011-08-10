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
                        stranded="logical",
                        bisulfite="logical"
                        ),
         prototype(samples=NULL,
                   alignments=NULL,
                   genome=NULL,
                   aligner=NULL,
                   index=NULL,
                   annotations=NULL,
                   paired=FALSE,
                   junction=FALSE,
                   stranded=FALSE,
                   bisulfite=FALSE
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

qProject <- function(sampleFile="Sample.txt", genome=".", annotationFile="Annotation.txt", aligner="Rbowtie", projectName="qProject", path=".", lib.loc=NULL)
{
    .progressReport("Gathering file path information", phase=-1)
    samples <- readSamples(sampleFile)
    annotations <- readAnnotations(annotationFile)
    genome <- checkGenome(genome, lib.loc=lib.loc)
    aligner <- loadAligner(aligner, lib.loc=lib.loc)
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
    cat("Options: paired=", object@paired, 
        "\n         junction=", object@junction,
        "\n         stranded=", object@stranded,
        "\n         bisulfite=", object@bisulfite, "\n", sep="")
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

qSaveProject <- function(project, filename)
{
    if(missing(filename))
        filename <- file.path(project@path, paste(project@id, "rds", sep="."))
    ## TODO calculate checksum of files or save modification date
    saveRDS(project, file=filename)
    return(filename)
}

qReadProject <- function(filename)
{
    project <- readRDS(file=filename)
    ## TODO check qProject with checksum, path ...
    return(project)
}
