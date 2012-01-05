# functions that export alignments or alignment-derived data to external files

bamfileToWig <- function(bamFname, wigFname=NULL, width=100L, shift=0L, paired=FALSE, maxHits=100L,
                         normFactor=1.0, trackname=NULL, color="#1B9E77", append=FALSE, quiet=FALSE) {
  if(!file.exists(bamFname))
    stop(bamFname," not found.")

  if(is.null(wigFname))
    wigFname <- paste(sub(".bam$","",bamFname),".wig")

  if(is.null(trackname))
    trackname <- sub(".bam$","",basename(bamFname))

  if(!quiet && paired && shift!=0L)
    warning("'shift' will not be used for paired-end alignments.")

  rgbcolor <- paste(as.vector(col2rgb(color[1], alpha=FALSE)),collapse=",")

  if(!quiet)
    .progressReport(paste("converting",bamFname,"to wig"), phase=-1)
  res <- .Call(.bamfile_to_wig, bamFname, wigFname, as.integer(width), as.integer(shift), as.logical(paired),
               as.integer(maxHits), as.numeric(normFactor), trackname, rgbcolor, as.logical(append))
  if(!quiet)
    .progressReport(paste("successfully created",wigFname), phase=1)

  return(invisible(res))
}
