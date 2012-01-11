\name{QuasR-package}
\alias{QuasR-package}
\alias{QuasR}
\docType{package}

\title{Quantify and Annotate Short Reads in R}

\description{
This package provides a pipeline for the quantification and analysis of Short Reads.
}

\details{
See packageDescription('QuasR') for package details.
}

\author{
Anita Lerch

Maintainer: Anita Lerch \email{anita.lerch@fmi.ch}
}

\seealso{
  \code{\linkS4class{qProject}}, \code{\link{qAlign}}, \code{\link{qCount}}
}

\examples{
\dontrun{
## Temporary directory, it will be deleted after closing the R session 
td <- tempdir()

## Example data
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")

## BSgenome name (see available.genomes()), a 'fasta' file 
## or a directory containing 'fasta' files
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"

## Create a QuasR project object which contains all information needed
## for the alignment. If the variable path is not defined, than output 
## will be written to the working directory.
project <- qProject(sampleFile, genome=genomeName, 
                    annotationFile=annotationFile, path=td)

## Execute the whole alignment including calculating the weights.
project <- qAlign(project)

## Get full report
qCount(project, "summary")

## Get count of a GenomicRange
gr <- GRanges(seq="chr2L", IRanges(1000, width=100000000), strand="+")
names(gr) <- "Example"
qCount(project, gr)

## Save the QuasR project object.
projectFilename <- qSaveProject(project)

## Restore the QuasR project object.
project <- qReadProject(projectFilename)
}
}

\keyword{package}

