setClassUnion(".data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion(".listOrNULL", c("list", "NULL"))
setClassUnion(".characterOrNULL", c("character", "NULL"))

setClass("qProject",
         representation(id="character",
                        name="character",
                        samples=".data.frameOrNULL",
                        alignments=".data.frameOrNULL",
                        genome=".listOrNULL",
                        aligner=".listOrNULL",
                        index=".listOrNULL",
                        annotations=".data.frameOrNULL",
                        alignmentParameter="character",
                        qc=".listOrNULL",
                        path="character",
                        cacheDir="character",
                        indexLocation=".characterOrNULL",
                        paired="logical",
                        junction="logical",
                        bisulfite="logical",
                        maxHits="integer"
                        ),
         prototype(samples=NULL,
                   alignments=NULL,
                   genome=NULL,
                   aligner=NULL,
                   index=NULL,
                   annotations=NULL,
                   qc=NULL,
                   indexLocation=NULL,
                   paired=FALSE,
                   junction=FALSE,
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

qProject <- function(sampleFile="Sample.txt", genome=".",
                     annotationFile=NULL, aligner="Rbowtie",
                     projectName="qProject", path=".", paired=FALSE,
                     junction=FALSE, bisulfite=FALSE, lib.loc=NULL,
                     indexLocation=NULL, maxHits=99L,
                     cacheDir=NULL,
                     alignmentParameter)
{
    .progressReport("Gathering file path information", phase=-1)
    sampleFile <- tools::file_path_as_absolute(sampleFile)
    samples <- .readSamples(sampleFile, paired=paired)
    if(is.null(annotationFile))
        annotations <- NULL
    else{
        annotationFile <- tools::file_path_as_absolute(annotationFile)
        annotations <- .readAnnotations(annotationFile)
    }
    genome <- .checkGenome(genome, lib.loc=lib.loc)
    aligner <- .loadAligner(aligner, lib.loc=lib.loc)
    if(!is.null(indexLocation))
        indexLocation <- tools::file_path_as_absolute(indexLocation)
    if(is.null(cacheDir))
        cacheDir <- tempdir()
    else
        cacheDir <- tools::file_path_as_absolute(cacheDir)
    qproject <- new("qProject", projectName, path, samples=samples,
                   annotations=annotations, genome=genome, aligner=aligner,
                   paired=paired, junction=junction,
                   bisulfite=bisulfite, indexLocation=indexLocation, maxHits=maxHits, cacheDir=cacheDir)
    if(missing(alignmentParameter) || is.null(alignmentParameter) || alignmentParameter == "" )
        qproject@alignmentParameter <- .createAlignmentParameters(qproject)
    qproject@alignments <- .loadAlignments(qproject)
    .progressReport(sprintf("Successfully created project '%s'", projectName), phase=1)
    return(qproject)
}

setMethod("show","qProject", function(object){
    cat("QuasRProject\n")
    cat("Project: " , object@name, "\n", sep="")
    cat("Options: paired=", object@paired,
        "\n         junction=", object@junction,
        "\n         bisulfite=", object@bisulfite,
        "\n         maxHits=", object@maxHits, "\n", sep="")
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
    if(!is(project, "qProject"))
        stop("The variable project is not a qProject")
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

getGenomeInformation <- function(qProject, ...){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")

    if(is.null(qProject@genome$sequenceInfo)){
        if(qProject@genome$bsgenome)
            qProject@genome$sequenceInfo <- seqlengths(.loadBSgenome(qProject@genome$name, ...))
        else {
            faList <- open(FaFileList(file.path(qProject@genome$dir, qProject@genome$files)))
            qProject@genome$sequenceInfo <- seqlengths(IRanges::unlist(GRangesList(IRanges::lapply(faList, scanFaIndex))))
        }
    }
    return(qProject@genome$sequenceInfo)
}

getAlignments <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@alignments)
}

getSamples <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@samples)
}

paired <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@paired)
}

# path <- function(qProject){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
#     return(qProject@path)
# }
