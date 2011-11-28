setClass("qProject",
         representation(env="environment",
                        id="character",
                        name="character")
)

setMethod("initialize", "qProject", function(.Object, name, ...){
    env <- new.env(parent=emptyenv())
    if(missing(name))
        name <- "qProject"
    id <- paste(name, format(Sys.time(), "%Y%m%d_%H%M%S"), sep="_")
    callNextMethod(.Object, env=env, name=name, id=id, ...)
})

qProject <- function(sampleFile="Sample.txt", genome=".",
                     annotationFile=NULL, aligner="Rbowtie",
                     projectName="qProject", bamfileDir=NULL, paired=FALSE,
                     junction=FALSE, bisulfite=FALSE, lib.loc=NULL,
                     indexLocation=NULL, maxHits=100L,
                     cacheDir=NULL,
                     alignmentParameter)
{
    .progressReport("Gathering filepath information", phase=-1)
    qproject <- new("qProject", projectName)
    if(bamfileDir == ".")
         bamfileDir <- tools::file_path_as_absolute(".")
    assign("bamfileDir", bamfileDir, qproject@env)
    sampleFile <- tools::file_path_as_absolute(sampleFile)
    samples <- .readSamples(sampleFile, paired=paired)
    assign("samples", samples, qproject@env)
    if(is.null(annotationFile)){
        annotations <- NULL
    }else{
        annotationFile <- tools::file_path_as_absolute(annotationFile)
        annotations <- .readAnnotations(annotationFile)
    }
    assign("annotations", annotations, qproject@env)
    genome <- .checkGenome(genome, lib.loc=lib.loc)
    assign("genome", genome, qproject@env)
    aligner <- .loadAligner(aligner, lib.loc=lib.loc)
    assign("aligner", aligner, qproject@env)
    if(!is.null(indexLocation))
        indexLocation <- tools::file_path_as_absolute(indexLocation)
    assign("indexLocation", indexLocation, qproject@env)
    if(is.null(cacheDir))
        cacheDir <- tempdir()
    else
        cacheDir <- tools::file_path_as_absolute(cacheDir)
    assign("cacheDir", cacheDir, qproject@env)
    assign("paired", paired, qproject@env)
    assign("junction", junction, qproject@env)
    assign("bisulfite", bisulfite, qproject@env)
    assign("maxHits", maxHits, qproject@env) 
    if(missing(alignmentParameter) || is.null(alignmentParameter) || alignmentParameter == "" )
        alignmentParameter <- .createAlignmentParameters(qproject)
    assign("alignmentParameter", alignmentParameter, qproject@env)
    alignments <- .loadAlignments(qproject)
    assign("alignments", alignments, qproject@env)
    assign("index", NULL, qproject@env)
    
#     if(suppressWarnings(require(parallel, quietly=TRUE)) && is.null(getOption("quasr.cluster"))){ # inherits(getOption("quasr.cluster"), "cluster")
#         options(quasr.cluster=makeCluster(2)) # TODO number of cluster 
#         clusterCall(getOption("quasr.cluster"), function() library("QuasR"))
#     #stopCluster(cl)
#     }
    .progressReport(sprintf("Successfully created project '%s'", projectName), phase=1)
    return(qproject)
}

setMethod("show","qProject", function(object){
    cat("QuasRProject\n")
    cat("Project: " , object@name, "\n", sep="")
    cat("Options: paired=", object@env$paired,
        "\n         junction=", object@env$junction,
        "\n         bisulfite=", object@env$bisulfite,
        "\n         maxHits=", object@env$maxHits,
        "\n         alignmentParameter=", object@env$alignmentParameter, "\n", sep="")
    cat("Genome:  ", object@env$genome$name, " is BSgenome=", object@env$genome$bsgenome, "\n", sep="")
    cat("Aligner: ", object@env$aligner$pkgname, " Version ", object@env$aligner$pkgversion, "\n", sep="")
    cat("Samples:\n", paste(object@env$samples$name, object@env$samples$filepath, sep="\t", collapse="\n"), "\n", sep="")
    cat("Annotations:\n", paste(object@env$annotations$feature, object@env$annotations$filepath, sep="\t", collapse="\n"), "\n", sep="")
    ## TODO write alignments in a nice way
    cat("Alignments:\n")
    lapply(names(object@env$alignments),
           function(index){
               cat(index, paste(as.matrix(object@env$alignments[index]), collapse="\n"), sep="\n")
           })
})

qSaveProject <- function(project, filename)
{
    if(!is(project, "qProject"))
        stop("The variable project is not a qProject")
    if(missing(filename))
        filename <- file.path(project@env$bamfileDir, paste(project@id, "rds", sep="."))
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

    if(is.null(qProject@env$genome$sequenceInfo)){
        if(qProject@env$genome$bsgenome)
            qProject@env$genome$sequenceInfo <- seqlengths(.loadBSgenome(qProject@env$genome$name, ...))
        else {
            faList <- open(FaFileList(file.path(qProject@env$genome$dir, qProject@env$genome$files)))
            qProject@env$genome$sequenceInfo <- seqlengths(IRanges::unlist(GRangesList(IRanges::lapply(faList, scanFaIndex))))
        }
    }
    return(qProject@env$genome$sequenceInfo)
}

getAlignments <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@env$alignments)
}

getSamples <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@env$samples)
}

paired <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@env$paired)
}

bamfileDir <- function(qProject){
    if(!is(qProject, "qProject"))
        stop("The object '", class(qProject), "' is not a 'qProject' object.")
    return(qProject@env$bamfileDir)
}
