\name{qAlign}
\alias{qAlign}

\title{QuasR Align}

\description{
Align sequences from a \sQuote{qProject} hierarchical to a genome, exon-exon-junction and auxiliary sequence.
}

\usage{
    qAlign(qproject, lib=NULL, lib.loc=NULL)
}

\arguments{
    \item{qProject}{A \sQuote{qProject} object containing all meta-data needed for the alignment.}    
    \item{lib}{A character() vector giving the library directories where to install the packages. Recycled as needed. If missing, defaults to the first element of \code{.libPaths()}.}
    \item{lib.loc}{A character() vector describing the location of R library trees to search through, or \sQuote{NULL}.  The default value of \sQuote{NULL} corresponds to all libraries currently known. Non-existent library trees are silently ignored.}
}

\details{
 TODO
}

\value{
    A modified \sQuote{qProject} object containing the location of the \sQuote{BAM} files is returned.
}

\author{Anita Lerch}

\seealso{
    \code{\linkS4class{qProject}}, \code{\link{qCount}}
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
}  
}

\keyword{methods}
\keyword{manip}