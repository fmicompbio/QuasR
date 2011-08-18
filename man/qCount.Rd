\name{qCount}
\alias{qCount}

\title{QuasR Count}

\description{
Retrieve the count of the alignments from a \sQuote{qProject} that are overlapped with the specified region. 
Function up to now not implemented.
}

\usage{
    qCount(qProject, gRange, stranded=FALSE, \dots)
}

\arguments{
    \item{qProject}{A \sQuote{qProject} object containing all meta-data needed for the alignment.}    
    \item{gRange}{TODO}
    \item{stranded}{logical(1).}
    \item{\dots}{logical(1).}
}

\details{
 TODO
}

\value{
    The modified \sQuote{GRanges} object containing the count overlapped with the specified region.
}

\author{Anita Lerch}

\seealso{
    \code{\linkS4class{qProject}}, \code{\link{qAlign}}
}

\examples{
\dontrun{
td <- tempdir()
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
aligner <- "Rbowtie"
project <- qProject(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

project <- qAlign(project)

project <- qCount(project, GRanges(seq="chr2L", IRanges(1000, width=500), strand="+"))
}  
}

\keyword{methods}
\keyword{manip}
