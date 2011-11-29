qExportWig <- function(qproject, file=NULL, combined=TRUE, width=100L, shift=0L, maxHits=100L, normalize=TRUE, tracknames=NULL,
                       colors=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"))
{
  if(!is(qproject, "qProject"))
    stop("The object '", class(qproject), "' is not a 'qProject' object.")
  
  bamfiles <- qproject@env$alignments$genome
  n <- length(bamfiles)
  idx <- !is.na(qproject@env$alignments$genome)

  if(is.null(file))
    if(combined)
      stop("A value for 'file' has to be provided when creating a combined wig file.")
    else
      file <- paste(sub(".bam$","",bamfiles),".wig")
  else
    if(combined && length(file)!=1)
      stop("Only a single value for 'file' can be provided when creating a combined wig file.")
    else if(!combined && length(file) != n)
      stop("The length of 'file' has to match the number of samples when creating individual wig files.")

  width <- as.integer(width)
  if(any(is.na(width)) || any(width<1))
    stop("'width' has to be a vector of positive integer values.")
  if(length(width)!=1 && length(width)!=n)
    stop("'width' has to contain either a single value or exactly one value per sample.")
  if(length(width)==1)
    width <- rep(width,n)

  shift <- as.integer(shift)
  if(any(is.na(shift)))
    stop("'shift' has to be a vector of integer values.")
  if(length(shift)!=1 && length(shift)!=n)
    stop("'shift' has to contain either a single value or exactly one value per sample.")
  if(length(shift)==1)
    shift <- rep(shift,n)

  if(length(maxHits)>1)
    warning("'maxHits' contains more than one value; only the first one will be used.")
  maxHits <- as.integer(maxHits)[1]
  if(is.na(maxHits) || maxHits < 1)
    stop("'maxHits' has to contain a single positive integer value.")
  if(maxHits > qproject@env$maxHits)
    stop(paste("'maxHits' cannot be greater than the 'maxHits' value defined for 'qproject' (",qproject@env$maxHits,")",sep=""))
  
  if(is.null(tracknames)) {
    tracknames <- as.character(qproject@env$samples$name)
    for(x in as.character(unique(tracknames))) {
      i <- tracknames == x
      if(sum(i) > 1)
        tracknames[i] <- paste(tracknames[i], ".", 1:sum(i), sep="")
    }
  }
  
  fact <- rep(1,n)
  if(normalize) {
    ms <- qproject@env$qc$mappingStats$genome
    if(is.null(ms)) {
      if(proj@env$aligner=="Rbowtie") {
        .progressReport("Collecting mapping statistics from bam files", phase=-1)
        mapdataL <- lapply(seq.int(n), function(i)
                           .getMappingStatsFromBam(tracknames[i], qproject@env$samples$filepath[i], qproject@env$samples$filetype[i], bamfiles[i]) )
        mapdata <- do.call(rbind, mapdataL)
        qproject@env$qc$mappingStats$genome <- mapdata
        .progressReport("", phase=1)

      } else { # if not using bowtie, cannot guarantee to recover unmapped/overmapped counts from bam files
        stop("mapping statistics are needed to use 'normalize'.")
      }
    }
    N <- rowSums(ms[,2:(ncol(ms)-1)]) # fist/last columns contains un-/overmapped
    fact <- min(N) /N
  }

  if(length(colors) < n)
    colors <- colorRampPalette(colors)(n)

  .progressReport("Starting creation of wig file(s)", phase=-1)
  if(combined && file.exists(file))
    unlink(file)
  lapply(which(idx), function(i)
         bamfileToWig(bamFname=bamfiles[i], wigFname=ifelse(combined[1], file[1], file[i]), width=width[i], shift=shift[i],
                      maxHits=maxHits, normFactor=fact[i], trackname=tracknames[i], color=colors[i], append=combined, quiet=TRUE))
  .progressReport("", phase=1)

  return(invisible(file))
}

.getMappingStatsFromBam <- function(nm, fnameSeq, typeSeq, fnameBam) {
  # get mapping statistics from bam file (IH and XM tags) and total number of reads (n)
  n <- switch(typeSeq, fasta=length(fasta.info(fnameSeq, use.names=FALSE)), fastq=fastq.geometry(fnameSeq)[1], NA)
  if(is.na(n))
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
