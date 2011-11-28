\name{bamfileToWig}
\alias{bamfileToWig}

\title{QuasR Align}

\description{
Create a fixed-step wig file (\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}) from the alignments in a bam file.
}

\usage{
    bamfileToWig(bamFname, wigFname=NULL, width=100L, shift=0L, trackname=NULL)
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
    \item{trackname}{A character(1) vector with the name of the
    track to appear in the track header. If \sQuote{NULL}, the basename
    of the \sQuote{bamFname} with ".bam" removed will be used.}
}

\details{
 The target sequence (e.g. the genome) is tiled with sequential windows
 of size \sQuote{width}, and alignments in the bam file are assigned to
 these windows according to their leftmost coordinate +
 \sQuote{shift}. In most cases, this simple approach will produce
 very similar results as if each alignment would be shifted by
 \sQuote{shift} bases towards its 3'-end.

 The resulting wig file can be compressed and is suitable for uploading
 to your favortie genome browser (e.g. UCSC or Ensembl).
}

\value{
    The file name of the generated wig file.
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
