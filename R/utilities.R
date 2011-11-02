.requirePkg <- function(pkgname, strict=TRUE, lib.loc=NULL)
{
    ## check package installed
    if(!require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc))
        stop("Can not find the '", pkgname, "' package.")
    ## check package version
    current <- installed.packages(lib.loc=lib.loc)[pkgname, 'Version']
    pkgs <- unlist(strsplit(installed.packages(lib.loc=lib.loc)['QuasR','Suggests'], ","))
    targetpkg <- grep(pkgname, pkgs, value=T)
    ## check if target is equal current version
    if(length(grep(paste("\\(== ", current, "\\)", sep=""), targetpkg)) == 0){
        if(strict)
            stop(gettextf("package '%s %s' was found, but '%s' is required by 'QuasR'",
                          pkgname, current, targetpkg))
        else
            warning(gettextf("package '%s' %s was found, but '%s' is required by 'QuasR'",
                             pkgname, current, targetpkg))
    }
}

.baseFileName <- function(filename)
{
    names <- strsplit(basename(filename),"\\.[^.]*$")
    names[filename == ""] <- "" # fix for empty filename
    return(unlist(names))
}

.fileExtension <- function(filename)
{
    pos <- regexpr("\\.[^.]*$", filename)
    pos[pos==-1] <- nchar(filename[pos==-1])
    extension <- substr(filename, pos+1, nchar(filename))
    return(extension)
}

.fileType <- function(filename)
{
    type <- as.character(.fileExtension(filename))
    type[ type %in% c("fa", "fna", "fasta") ] <- "fasta"
    type[ type %in% c("fq", "fastq") ] <- "fastq"
    return(type)
}

.createBamFilename <- function(path, sampleFilename, genomeName)
{
    if(is.null(path))
#         fn <- file.path(dirname(sampleFilename), sprintf("%s-%s.bam", .baseFileName(sampleFilename), genomeName))
        fn <- tempfile(tmpdir=dirname(sampleFilename),
                   pattern=sprintf("%s-%s-", .baseFileName(sampleFilename), genomeName), 
                   fileext=".bam")
    else 
#         fn <- file.path(path, sprintf("%s-%s.bam", .baseFileName(sampleFilename), genomeName))
        fn <- tempfile(tmpdir=path,
                   pattern=sprintf("%s-%s-", .baseFileName(sampleFilename), genomeName), 
                   fileext=".bam")
    return(fn)
}

.getBamFile <- function(listBamFilenames, genomeName, alignerParameter)
{
    # TODO qproject parameter  
    bfh <- scanBamHeader(listBamFilenames)
    bfhIdx <- unlist(lapply(bfh, function(x){
        idx <- grep("ID:QuasR", x$text)
        qTag <- x$text[[idx]]
        at <- unlist(strsplit(qTag[grep("at:", qTag)], ":"))[2]
        ap <- unlist(strsplit(qTag[grep("ap:", qTag)], ":"))[2]
        all(at == genomeName,  ap == paste(alignerParameter, collapse=","))
    }))
    if(sum(bfhIdx) > 1L)
        stop("More than one bamfile found '", paste(listBamFilenames[bfhIdx], collapse="', '"), "'.")
    if(sum(bfhIdx) < 1L)
        return(NA)
    else
        return(listBamFilenames[bfhIdx])
}

.progressReport <- function(msg, phase=0)
{
    qTag <- "[QuasR]"
    if(phase==0)
        cat("done\n")
    if(phase<=0)
        msg <- paste(msg, "...", sep="")
    if(phase>0)
        cat("done\n")
    cat(qTag, msg)
    if(phase>0)
        cat("\n")
}

.multiToSingleFasta <- function(inFiles, outFile){
  if(missing(outFile) || is.null(outFile))
    outFile <- sprintf("%s.fa", tempfile())
  append <- FALSE
  for(file in inFiles){
    seq <- read.DNAStringSet(filepath=file)
    write.XStringSet(seq, filepath=outFile, format="fasta", append=append)
    append <- TRUE
  }
  return(outFile)
}

.mapSeqnames <- function(genomeSeqnames, annotationSeqnames){
    if(!all(annotationSeqnames %in% genomeSeqnames)){
        ## try to map mitochondrion
        annotationSeqnames[ annotationSeqnames=="dmel_mitochondrion_genome" ] <- "M"
        annotationSeqnames[ annotationSeqnames=="MT" ] <- "M"
        ## add "chr" string
        annotationSeqnames <- paste("chr", annotationSeqnames, sep="")
        if(!all(annotationSeqnames %in% genomeSeqnames))
            stop("Can't match sequence names")
    }
    return(annotationSeqnames)
}
