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

qProject <- function(sampleFile, genome,
                     auxiliaryFile=NULL, aligner="Rbowtie",
                     projectName="qProject", bamfileDir=NULL,
                     splicedAlignment=FALSE, bisulfite=FALSE, lib.loc=NULL,
                     indexLocation=NULL, maxHits=100L,
                     cacheDir=NULL,
                     alignmentParameter)
{
    if(missing(sampleFile))
        stop("Missing 'sampleFile' parameter.")
    if(missing(genome))
        stop("Missing 'genome' parameter.")

    .progressReport("Gathering filepath information", phase=-1)
    qproject <- new("qProject", projectName)
    if(!is.null(bamfileDir)){
        ## convert to absolut filename
        if(bamfileDir == ".")
            bamfileDir <- tools::file_path_as_absolute(".")
        if(substr(bamfileDir, 1, 1) != "/")
            bamfileDir <- file.path(getwd(), bamfileDir)
    }
    assign("bamfileDir", bamfileDir, qproject@env)
    sampleFile <- tools::file_path_as_absolute(sampleFile)
    samples <- .readSamples(sampleFile)
    assign("samples", samples, qproject@env)
    paired <- ifelse(length(samples) == 3, FALSE, TRUE)
    assign("paired", paired, qproject@env)    
    if(is.null(auxiliaryFile)){
        auxiliaries <- NULL
    }else{
        auxiliaryFile <- tools::file_path_as_absolute(auxiliaryFile)
        auxiliaries <- .readAuxiliaries(auxiliaryFile)
    }
    assign("auxiliaries", auxiliaries, qproject@env)
    genome <- .checkGenome(genome, lib.loc=lib.loc)
    assign("genome", genome, qproject@env)
    aligner <- .loadAligner(aligner, lib.loc=lib.loc)
    assign("aligner", aligner, qproject@env)
    if(!is.null(indexLocation))
        indexLocation <- tools::file_path_as_absolute(indexLocation)
    assign("indexLocation", indexLocation, qproject@env)
    if(is.null(cacheDir)){
#         cacheDir <- tempfile(pattern="quasr_")
        cacheDir <- tempdir()
        dir.create(path=cacheDir, showWarnings=FALSE)
    }else
        cacheDir <- tools::file_path_as_absolute(cacheDir)
    assign("cacheDir", cacheDir, qproject@env)
    if(splicedAlignment)
        stop("Spliced alignment mode is not implemented yet.")
    assign("splicedAlignment", splicedAlignment, qproject@env)
    if(bisulfite)
        stop("Bisulfite mode is not implemented yet.")
    assign("bisulfite", bisulfite, qproject@env)
    assign("maxHits", maxHits, qproject@env) 
    if(missing(alignmentParameter) || is.null(alignmentParameter) || alignmentParameter == "" )
        alignmentParameter <- .createAlignmentParameters(qproject)
    assign("alignmentParameter", alignmentParameter, qproject@env)
    alignments <- .loadAlignments(qproject)
    assign("alignments", alignments, qproject@env)
    assign("index", NULL, qproject@env)
    
    if(suppressWarnings(require(parallel, quietly=TRUE)) && is.null(getOption("quasr.cluster"))){ # inherits(getOption("quasr.cluster"), "cluster")
        if(is.null(getOption("quasr.clusterSize")))
            options(quasr.clusterSize=2)
        options(quasr.cluster=makeCluster(getOption("quasr.clusterSize")))
        clusterCall(getOption("quasr.cluster"), function() library("QuasR"))
        #stopCluster(cl)
    }
    .progressReport(sprintf("Successfully created project '%s'", projectName), phase=1)
    return(qproject)
}

setMethod("show","qProject", function(object){
    cat("Project: " , object@name, "\n", sep="")
    cat("Options: paired=", object@env$paired,
        "\n         splicedAlignment=", object@env$splicedAlignment,
        "\n         bisulfite=", object@env$bisulfite,
        "\n         maxHits=", object@env$maxHits,
        "\n         alignmentParameter=", object@env$alignmentParameter, "\n", sep="")
    cat("Genome:  ",
        if(object@env$genome$bsgenome) object@env$genome$name
        else .truncPath(object@env$genome$name, getOption("width")-26),
        " (BSgenome=", object@env$genome$bsgenome, ")\n", sep="")
    cat("Aligner: ", object@env$aligner$pkgname, " Version ", object@env$aligner$pkgversion, "\n", sep="")
    maxlen <- max(nchar(object@env$samples$name),0)
    cat("Samples:\n", paste(sprintf(" %*s  %s\n", maxlen, object@env$samples$name,
                                    .truncPath(object@env$samples$filepath, getOption("width")-maxlen-3)), collapse=""), sep="")
    maxlen <- max(nchar(object@env$auxiliaries$feature),0)
    cat("Auxiliaries:\n",
        if(is.null(object@env$auxiliaries)) " none\n"
        else paste(sprintf(" %*s  %s\n", maxlen, object@env$auxiliaries$feature,
                           .truncPath(object@env$auxiliaries$filepath, getOption("width")-maxlen-3)), collapse=""), sep="")
    cat("Alignments:\n")
    lapply(names(object@env$alignments),
           function(index){
               cat(paste(" ",index,sep=""),
                   paste("  ", .truncPath(object@env$alignments[index],getOption("width")), sep="", collapse="\n"),
                   sep="\n")
           })
    return(invisible(NULL))
})

qSaveProject <- function(qproject, filename)
{
    if(!is(qproject, "qProject"))
        stop("The variable project is not a qProject")
    if(missing(filename))
        filename <- file.path(paste(qproject@id, "rds", sep="."))
    ## TODO calculate checksum of files or save modification date
    saveRDS(qproject, file=filename)
    return(filename)
}

qReadProject <- function(filename)
{
    qproject <- readRDS(file=filename)
    ## TODO check qProject with checksum, path ...
    return(qproject)
}

# getGenomeInformation <- function(qProject, ...){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
# 
#     if(is.null(qProject@env$genome$sequenceInfo)){
#         if(qProject@env$genome$bsgenome)
#             qProject@env$genome$sequenceInfo <- seqlengths(.loadBSgenome(qProject@env$genome$name, ...))
#         else {
#             faList <- open(FaFileList(file.path(qProject@env$genome$dir, qProject@env$genome$files)))
#             qProject@env$genome$sequenceInfo <- seqlengths(IRanges::unlist(GRangesList(IRanges::lapply(faList, scanFaIndex))))
#         }
#     }
#     return(qProject@env$genome$sequenceInfo)
# }
# 
# getAlignments <- function(qProject){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
#     return(qProject@env$alignments)
# }
# 
# getSamples <- function(qProject){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
#     return(qProject@env$samples)
# }
# 
# paired <- function(qProject){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
#     return(qProject@env$paired)
# }
# 
# bamfileDir <- function(qProject){
#     if(!is(qProject, "qProject"))
#         stop("The object '", class(qProject), "' is not a 'qProject' object.")
#     return(qProject@env$bamfileDir)
# }
