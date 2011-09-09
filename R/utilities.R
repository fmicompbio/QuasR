.requirePkg <- function(pkgname, strict=TRUE, lib.loc=NULL)
{
    ## check package installed
    if(!require(pkgname, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc))
        stop("Can not find the '", pkgname, "' package.")
    ## check package version
    current <- installed.packages()[pkgname, 'Version']
    pkgs <- unlist(strsplit(installed.packages()['QuasR','Suggests'], ","))
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
    ext <- substr(filename, pos+1, nchar(filename))
    return(ext)
}

.fileType <- function(filename)
{
    type <- as.character(.fileExtension(filename))
    #type[ type %in% c("fa", "fna", "fasta") ] <- "fasta"
    return(type)
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
            error("Can't match sequence names")
    }
    return(annotationSeqnames)
}   
