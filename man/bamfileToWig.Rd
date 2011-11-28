\name{bamfileToWig}
\alias{bamfileToWig}

\title{Create wig file from bam file}

\description{
Create a fixed-step wig file from the alignments in a bam file.
}

\usage{
    bamfileToWig(bamFname, wigFname=NULL, width=100L, shift=0L,
    normFactor=1.0, trackname=NULL, color="#1B9E77", append=FALSE)
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
    fragment length, see details).}
  \item{normFactor}{A scalar value defining the normalization
    factor. The coverage for each window will be multiplied by \sQuote{normFactor}.}
  \item{trackname}{A character(1) vector with the name of the
    track to appear in the track header. If \sQuote{NULL}, the basename
    of the \sQuote{bamFname} with ".bam" removed will be used.}
  \item{color}{A character(1) vector with the color for the track.}
  \item{append}{If \sQuote{TRUE}, the generated track is appended to a
    the WIG file \sQuote{wigFname}. If \sQuote{wigFname} does not exist,
    it will be created.}
}

\details{
  Using a BAM file as input, the function creates a WIG file
  (\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}) containing
  the coverage of alignments along the genome.
  
  The genome is tiled with sequential windows of size \sQuote{width},
  and alignments in the BAM file are assigned to these windows according
  to their 5'-end coordinate shifted by \sQuote{shift} towards the
  3'-end. The 5'-end is the leftmost coordinate for plus-strand
  alignments, and the rightmost coordinate for minus-strand alignments.

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
