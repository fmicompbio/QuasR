# functions that export alignments or alignment-derived data to external files

bamfileToWig <- function(bamFname, wigFname=NULL, width=100L, shift=0L, normFactor=1.0, trackname=NULL, color="#1B9E77", append=FALSE) {
  if(!file.exists(bamFname))
    stop(bamFname," not found.")

  if(is.null(wigFname))
    wigFname <- paste(sub(".bam$","",bamFname),".wig")

  if(is.null(trackname))
    trackname <- sub(".bam$","",basename(bamFname))

  rgbcolor <- paste(as.vector(col2rgb(color[1], alpha=FALSE)),collapse=",")
  
  .progressReport(paste("converting",bamFname,"to wig"), phase=-1)
  res <- .Call(.bamfile_to_wig, bamFname, wigFname, as.integer(width), as.integer(shift),
               as.numeric(normFactor), trackname, rgbcolor, as.logical(append))
  .progressReport(paste("successfully created",wigFname), phase=1)

  return(invisible(res))
}
