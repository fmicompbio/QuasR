\name{qAlign}
\alias{qAlign}

\title{
    QuasR Align
}

\description{
    Align sequences from a \sQuote{qProject} hierarchical to a genome, exon-junction and auxiliary sequence.
}

\usage{

    qAlign(qProject, \dots)

}

\arguments{

    \item{qProject}{A \sQuote{qProject} object containing all meta-data needed for the alignment.}
    
    \item{...}{Further arguments to be passed to TODO.}
    
}

\details{
 TODO
}

\value{
    A modified \sQuote{qProject} object containing the location of the \sQuote{BAM} files is returned.
  
}

\author{
  Anita Lerch
}

\seealso{
  \code{\linkS4class{qProject}},
}

\examples{

  td <- tempdir()
  sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
  annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
  genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
  aligner <- "Rbowtie"
  project <- qProject(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

  project <- qAlign(project)
  
}

\keyword{methods}
\keyword{manip}