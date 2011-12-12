\name{qExportWig}
\alias{qExportWig}

\title{QuasR wig file export}

\description{
Create a fixed-step wig file from the alignments in the genomic bam
files of the \sQuote{QuasR} project.
}

\usage{
    qExportWig(qproject, file=NULL, combined=TRUE, width=100L, shift=0L, maxHits=100L,
    normalize=TRUE, tracknames=NULL, colors=c("#1B9E77",
    "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"))
}

\arguments{
  \item{qProject}{A \sQuote{qProject} object containing all information
    from the current project.}
  \item{file}{A character vector with the name(s) for the wig file(s) to
    be generated. Either \sQuote{NULL}, a single file name for
    \sQuote{combined=TRUE}, or a vector of the same length as the number
    of samples (bam files) in \sQuote{qproject} for
    \sQuote{combined=FALSE}. If \sQuote{NULL} (only possible if
    \sQuote{combinded=FALSE}), the wig file names are generated from the
    names of the genomic bam files with ".wig" replacing the ".bam" extension.}
  \item{combined}{If TRUE, a single multi-track wig file will be
    generated containing a track for each of the samples (genomic bam
    files) in \sQuote{qproject}.}
  \item{width}{Either a vector or a scalar value defining the window and
    step size for the wig file(s). If \sQuote{length(width)>1}, the
    length must match the number of samples in \sQuote{qproject}, and
    the i-th sample will be converted to wig using the value in
    \sQuote{width[i]}. \sQuote{width} will be coerced to integer().}
  \item{shift}{Either a vector or a scalar value defining the read shift (e.g. half of
    fragment length, see details). If \sQuote{length(shift)>1}, the
    length must match the number of samples in \sQuote{qproject}, and
    the i-th sample will be converted to wig using the value in
    \sQuote{shift[i]}. \sQuote{shift} will be coerced to integer().}
  \item{maxHits}{Only include alignments of reads with no more than
    \sQuote{maxHits} hits to the genome.}
  \item{normalize}{If TRUE, the output values in the wig file(s)
    will be normalized to improve comparability of multiple samples in
    \sQuote{qproject} (see details).}
  \item{tracknames}{A character vector with the names of the
    tracks to appear in the track header. If \sQuote{NULL}, the sample
    names in \sQuote{qproject} will be used.}
  \item{colors}{A character vector with R color names to be used for the tracks.}
}

\details{
  \sQuote{qExportWig()} is a wrapper function around
  \sQuote{bamfileToWig} and uses the genome bam files in
  \sQuote{qproject} as input to create a combined wig file
  (\sQuote{combined=TRUE}) or individual wig files
  (\sQuote{combined=FALSE}) with the number of alignments per window of
  \sQuote{width} base pairs. See
  (\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}) for the
  definition of the wig format.
  
  The genome is tiled with sequential windows of size \sQuote{width},
  and alignments in the bam file are assigned to these windows according
  to their 5'-end coordinate shifted by \sQuote{shift} towards the
  3'-end. The 5'-end is the leftmost coordinate for plus-strand
  alignments, and the rightmost coordinate for minus-strand alignments.

  For \sQuote{normalize=TRUE}, the number of alignments per window
  \sQuote{n} for the sample \sQuote{i} are linearly scaled to the
  minimum total number of alignments over all samples in
  \sQuote{qproject} according to:
  \sQuote{n_norm = n /N[i] *mean(N)}
  where \sQuote{n_norm} is the normalized number of alignments in the
  window, \sQuote{N} is a vector with the total number of alignments for
  each sample.

  The resulting wig file can be compressed and is suitable for uploading
  as a custom track to your favortie genome browser (e.g. UCSC or Ensembl).
}

\value{
  (invisible) The file name of the generated wig file(s).
}

\author{Michael Stadler}

\seealso{
    \code{\linkS4class{qProject}}, \code{\link{bamfileToWig}}
}

\examples{
\dontrun{
    qExportWig(qproject, file="combined_100nt.wig", combined=TRUE,
width=100L, shift=0L, normalize=TRUE)}
}

\keyword{utilites}
