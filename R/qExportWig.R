qExportWig <- function(qproject, file=NULL, combined=TRUE, width=100L, shift=0L, normalize=TRUE, tracknames=NULL,
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
    else if(length(file) != length(bamfiles))
      stop("The length of 'file' has to match the number of samples when creating individual wig files.")

  if(normalize && is.null(ms <- qproject@env$qc$mappingStats$genome))
    stop("mapping statistics are needed to use 'normalize'.")

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
  
  fact <- rep(1,n)
  if(normalize) {
    N <- rowSums(ms[,2:(ncol(ms)-1)]) # fist/last columns contains un-/overmapped
    fact <- min(N) /N
  }

  if(is.null(tracknames)) {
    tracknames <- as.character(qproject@env$samples$name)
    for(x in as.character(unique(tracknames))) {
      i <- tracknames == x
      if(sum(i) > 1)
        tracknames[i] <- paste(tracknames[i], ".", 1:sum(i), sep="")
    }
  }
  
  if(length(colors) < n)
    colors <- colorRampPalette(colors)(n)

  .progressReport("Starting creation of wig file(s)", phase=-1)
  if(combined && file.exists(file))
    unlink(file)
  lapply(which(idx), function(i)
         bamfileToWig(bamFname=bamfiles[i], wigFname=ifelse(combined[1], file[1], file[i]), width=width[i], shift=shift[i],
                      normFactor=fact[i], trackname=tracknames[i], color=colors[i], append=combined))
  .progressReport(phase=0)

  return(invisible(file))
}
