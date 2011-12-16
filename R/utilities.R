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

.createBamFilename <- function(bamfileDir, sampleFilename, genomeName)
{
    if(is.null(bamfileDir))
#         fn <- file.path(dirname(sampleFilename), sprintf("%s-%s.bam", .baseFileName(sampleFilename), genomeName))
        fn <- tempfile(tmpdir=dirname(sampleFilename),
                   pattern=sprintf("%s-%s-", .baseFileName(sampleFilename), genomeName), 
                   fileext=".bam")
    else 
#         fn <- file.path(path, sprintf("%s-%s.bam", .baseFileName(sampleFilename), genomeName))
        fn <- tempfile(tmpdir=bamfileDir,
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
        if(length(idx) == 0)
            return(FALSE)
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

.progressReport <- function(msg, phase=0, qTag="[QuasR]", quiet=getOption("quasr.quiet"))
{
    if(!quiet){
        if(phase==0)
            cat("done\n")
        if(phase<=0)
            msg <- paste(msg, "...", sep="")
        if(phase>0)
            cat("done\n")
        if(msg != ""){
            cat(qTag, msg)
            if(phase>0)
                cat("\n")
        }
    }
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

.getMappingTagFromBam <- function(fnameBam){
    bfh <- scanBamHeader(fnameBam)[[1]]
    idx <- grep("ID:QuasR", bfh$text)
    if(length(idx) == 0)
        return(NULL)
    qTag <- bfh$text[[idx]]
    ## extrect subtags
    ml <- unlist(strsplit(qTag[grep("ml:", qTag)], ":"))[2]
    mc <- unlist(strsplit(qTag[grep("mc:", qTag)], ":"))[2]
    ## create r structure
    mc <- as.integer(strsplit(mc,";")[[1]])
    names(mc) <- strsplit(ml,";")[[1]]
    return(mc)
}

.getMappingStatsFromBam <- function(nm, fnameSeq, typeSeq, fnameBam) {
  # get mapping statistics from bam file (IH and XM tags) and total number of reads (n)
  if(typeSeq=="fasta")
    n <- length(fasta.info(fnameSeq, use.names=FALSE))
  else if(typeSeq=="fastq")
    n <- fastq.geometry(fnameSeq)[1]
  else
    stop(paste("Don't know how to count sequences in a '",typeSeq,"' file.",sep=""))

  ih <- scanBam(fnameBam, param=ScanBamParam(what=character(0), tag="IH"))[[1]]$tag$IH
  ihT <- table(ih)
  rm(ih)
  nhits <- as.numeric(names(ihT))
  ihTF <- numeric(max(nhits)+2) # from zero, 1:maxHits, >maxHits
  ihTF[nhits+1] <- ihT/nhits
  if("0" %in% names(ihT)) {
    xm <- scanBam(fnameBam, param=ScanBamParam(what=character(0), tag="XM", flag=scanBamFlag(isUnmappedQuery=TRUE)))[[1]]$tag$XM
    ihTF[max(nhits)+2] <- sum(xm>0)
    ihTF[1] <- ihT["0"] - ihTF[max(nhits)+2]
    rm(xm)
  }
  return(matrix(ihTF, nrow=1, dimnames=list(nm, c(as.character(0:max(nhits)), paste(">",max(nhits),sep="")))))
}
