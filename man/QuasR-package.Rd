\name{QuasR-package}
\alias{QuasR-package}
\alias{QuasR}
\docType{package}

\title{
Quantify and Annotate Short Reads in R
}

\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}

\details{
\tabular{ll}{
Package: \tab QuasR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-03-23\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}

\author{
Anita Lerch

Maintainer: Anita Lerch <anita.lerch@fmi.ch>
}

\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }

\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}

\examples{

    ## Temporary directory, will be deleted after closing R 
    td <- tempdir()
    
    ## Example data
    sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
    annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
    
    ## BSgenome name, a 'fasta' file or a directory containing 'fasta' files
    genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
    
    ## Create a QuasR project object which contains all information needed for the alignment.
    ## If the variable path is not defined, than output will be written to the working directory.
    project <- QuasRProject(sampleFile, genome=genomeName, annotationFile=annotationFile, path=td)

    ## Execute the whole alignment including calculating the weights
    project <- quasrAlign(project)

}
