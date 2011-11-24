\name{bamfileToWig}
\alias{bamfileToWig}

\title{QuasR Align}

\description{
Create a fixed-step wig file (\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}) from the alignments in a bam file.
}

\usage{
    bamfileToWig(bamFname, wigFname=NULL, width=100L)
}

\arguments{
    \item{bamFname}{A character(1) vector with the file name of the
    input bam file.}    
    \item{wigFname}{A character(1) vector with the file name of the
    output wig file. If \sQuote{NULL}, \sQuote{bamFname} with an added
    ".wig" will be used.}
    \item{width}{A scalar value defining the window size and step size
    for the wig file (will be coerced to integer()).}
}

\details{
 TODO
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
    bamfileToWig(bamFname="H3K27me3.bam", wigFname="H3K27me3_w100.wig", width=100L)}  
}

\keyword{utilites}
