\name{qCount}
\alias{qCount}
\alias{qCount,qProject,GRanges-method}
\alias{qCount,qProject,character-method}
\alias{qCount,qProject,missing-method}

\title{QuasR Count}

\description{
Retrieve the count of the alignments from a \sQuote{qProject} that are overlapped with the specified region. 
Function up to now not implemented.
}

\usage{
    qCount(qproject, query, stranded=FALSE, collapseSamples=TRUE, collapseQuery=c("gene", "transcript", "exon"), overlap=c("any", "within", "startwithin", "endwithin"), maxHits=NULL, shift=0, minoverlap=1L)
    \S4method{qCount}{qProject,missing}(qproject, query, stranded, collapseSamples, collapseQuery, overlap, maxHits, shift, minoverlap)
    \S4method{qCount}{qProject,character}(qproject, query, stranded, collapseSamples, collapseQuery, overlap, maxHits, shift, minoverlap)
    \S4method{qCount}{qProject,GRanges}(qproject, query, stranded, collapseSamples, collapseQuery, overlap, maxHits, shift, minoverlap)
}

\arguments{
    \item{qproject}{A \sQuote{qProject} object containing all meta-data needed for the alignment.}    
    \item{query}{TODO}
    \item{stranded}{logical(1).}
    \item{collapseSamples}{logical(1).}
    \item{collapseQuery}{logical(1).}
    \item{overlap}{.}
    \item{maxHits}{.}
    \item{shift}{integer vector.}
    \item{minoverlap}{integer vector.}
}

\details{
 TODO
}

\value{
    A data.frame object containing the count overlapped with the specified region. The rownames are the same as the names of the GRange object and the colnames are the sample names if collapseSamples is \code{TRUE} otherwise the fasta/q filenames.
}

\author{Anita Lerch}

\seealso{
    \code{\linkS4class{qProject}}, \code{\link{qAlign}}
}

\examples{
\dontrun{
td <- tempdir()
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
aligner <- "Rbowtie"
project <- qProject(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

project <- qAlign(project)

project <- qCount(project, GRanges(seq="chr2L", IRanges(1000, width=500), strand="+"))
}  
}

\keyword{methods}
\keyword{manip}
