\name{bamfileToWig}
\alias{bamfileToWig}

\title{Create wig file from bam file}

\description{
Create a fixed-step wig file from the alignments in a bam file.
}

\usage{
    bamfileToWig(bamFname, wigFname=NULL, width=100L, shift=0L,
    paired=FALSE, maxHits=100L, normFactor=1.0, trackname=NULL,
    color="#1B9E77", append=FALSE, quiet=FALSE)
}

\arguments{
  \item{bamFname}{A character(1) vector with the file name of the
    input bam file.}    
  \item{wigFname}{A character(1) vector with the file name of the
    output wig file. If \sQuote{NULL}, \sQuote{bamFname} with an added
    ".wig" will be used.}
  \item{width}{A scalar value defining the window size and step size
    for the wig file (will be coerced to integer()).}
  \item{shift}{A scalar value defining the read shift (e.g. half of
    fragment length, see details). \sQuote{shift} is ignored if
    \sQuote{paired=TRUE}.}
  \item{paired}{If \sQuote{TRUE}, only count paired alignments and
    assign them to the base in the middle (see details).}
  \item{maxHits}{Only include alignments of reads with no more than
    \sQuote{maxHits} hits to the genome. Warning: If
    \sQuote{bamfileToWig} is called directly (rather than by
    \sQuote{qExportWig}), there is no guarantee that all alignments up
    to \sQuote{maxHits} are contained in the bam file.}
  \item{normFactor}{A scalar value defining the normalization
    factor. The coverage for each window will be multiplied by \sQuote{normFactor}.}
  \item{trackname}{A character(1) vector with the name of the
    track to appear in the track header. If \sQuote{NULL}, the basename
    of the \sQuote{bamFname} with ".bam" removed will be used.}
  \item{color}{A character(1) vector with the color for the track.}
  \item{append}{If \sQuote{TRUE}, the generated track is appended to a
    the WIG file \sQuote{wigFname}. If \sQuote{wigFname} does not exist,
    it will be created.}
  \item{quiet}{If \sQuote{TRUE}, no progress reports will be printed to
    the console.}
}

\details{
  Using a BAM file as input, the function creates a WIG file
  (\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}) containing
  the coverage of alignments along the genome. Only alignments of reads
  with up to \sQuote{maxHits} hits will be included. Setting
  \sQuote{maxHits=1} will thus produce a wig file from uniquely mapped reads.
  
  The genome is tiled with sequential windows of size \sQuote{width},
  and alignments in the BAM file are assigned to these windows: If
  \sQuote{paired=FALSE} individual read alignments are assigned
  according to their 5'-end coordinate shifted by \sQuote{shift} towards
  the 3'-end (The 5'-end is the leftmost coordinate for plus-strand
  alignments, and the rightmost coordinate for minus-strand
  alignments). If \sQuote{paired=TRUE}, the paired-end alignments are
  assigned according to the base in the middle between the leftmost and
  rightmost coordinates of the aligned pair of reads. Each pair of reads
  is only counted once, and unpaired alignments are ignored.

  The resulting wig file can be compressed and is suitable for uploading
  as a custom track to your favortie genome browser (e.g. UCSC or Ensembl).
}

\value{
    (invisible) The file name of the generated wig file.
}

\author{Michael Stadler}

%\seealso{
%    \code{\linkS4class{qProject}}, \code{\link{qCount}}
%}

\examples{
\dontrun{
    bamfileToWig(bamFname="H3K27me3.bam",
                 wigFname="H3K27me3_w100.wig",
                 width=100L, shift=0,
                 trackname="H3K27me3")}  
}

\keyword{utilites}
