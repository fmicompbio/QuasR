# functions that export alignments or alignment-derived data to external files

bamfileToWig <- function(bamFname, wigFname=NULL, width=100L, shift=0L, trackname=NULL) {
  if(!file.exists(bamFname))
    stop(bamFname," not found.")

  if(is.null(wigFname))
    wigFname <- paste(sub(".bam$","",bamFname),".wig")

  if(is.null(trackname))
    trackname <- sub(".bam$","",basename(bamFname))

  .progressReport(paste("converting",bamFname,"to wig"), phase=-1)
  res <- .Call(.bamfile_to_wig, bamFname, wigFname, as.integer(width), as.integer(shift), trackname)
  .progressReport(paste("successfully created",wigFname), phase=1)

  return(res)
}
