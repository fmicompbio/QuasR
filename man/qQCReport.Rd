\name{qQCReport}
\alias{qQCReport}

\title{QuasR Quality Control Report}

\description{
Create the quality control report of the short reads and the alignments from a \sQuote{qProject}.
}

\usage{
    qQCReport(qproject, pdfFilename, ...)
}

\arguments{
    \item{qproject}{A \sQuote{qProject} object containing all meta-data needed for the alignment.}    
    \item{pdfFilename}{A character() vector of the output pdf filename. If missing the standard graphical device is used.}
    \item{...}{ lmat e.g. matrix(1:18, ncol=3) or rbind(1:2), breaks e.g. c(0,1,20,50)}
}

\details{}

\value{}

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

project <- qQCReport(project)
}  
}

\keyword{methods}
