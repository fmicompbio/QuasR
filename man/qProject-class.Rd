\name{qProject}

\docType{class}
\alias{qProject-class}
\alias{qProject}
\alias{show,qProject-method}
\alias{qSaveProject,qProject-method}
\alias{qSaveProject}
\alias{qReadProject}

\title{A meta-data class for deep sequencing analysis}
  
\description{
Container class for meta-data needed during the deep sequencing analysis.
}

\usage{

## Constructor
qProject(sampleFile="Sample.txt",
     genome=".", 
     annotationFile="Annotation.txt", 
     aligner="Rbowtie", 
     projectName="qProject", 
     path=".", 
     paired=FALSE, 
     junction=FALSE, 
     bisulfite=FALSE,
     lib.loc=NULL)

## Accessors
\S4method{show}{qProject}(object)
\S4method{paired}{qProject}(object)

## Methods
\S4method{qSaveProject}{qProject}(object)
qReadProject(filename)

}

\arguments{

    \item{sampleFile}{A character() vector of the sample file path name.}
    \item{genome}{A character() vector of the \sQuote{BSgenome} name or of the directory path names fasta files path names or of a single fasta file path names.}
    \item{annotationFile}{A character() vector of the annotation file path name.}  
    \item{aligner}{A character() vector of the name of the alinger package used for the alignments.}
    \item{projectName}{A character() vector of the project name, used for saving the \sQuote{qProject}.}
    \item{path}{A character() vector of the directory path names, at which the resulting files are stored.}
    \item{lib.loc}{A character() vector describing the location of R library trees to search through, or \sQuote{NULL}.  The default value of \sQuote{NULL} corresponds to all libraries currently known.  Non-existent library trees are silently ignored.}
    \item{object}{A qProject() object.}
    \item{filename}{A character() vector of the \sQuote{qProject} rds-file.}
    \item{paired}{logical(1). If \code{TRUE} then the samples in the \sQuote{qProject} are handled as paired-end reads. Up to now not implemented.}
    \item{junction}{A logical(1). If \code{TRUE} then the samples are additionaly aligned to the Exon-exon-junction database. Up to now not implemented.}
    \item{bisulfite}{logical(1). If \code{TRUE} then the samples in the \sQuote{qProject} are handled as bisultite reads. Up to now not implemented.}
}

\section{Slots}{
    \describe{
    
        \item{name}{A character() vector of the project name.}        
        \item{id}{A character() vector of a unique identifier for the project, used for saving the \sQuote{qProject}.}        
        \item{path}{A character() vector of the directory path names, at which the resulting files are stored.}      
        \item{samples}{A data.frame() containing the path to of the sequence filenames and the sample name.}
        \item{annotations}{A data.frame() containing the path to of the annotation filenames and the feature names.}
        \item{genome}{The BSgenome object or a list containing information of a genomic fasta sequence.}
        \item{aligner}{The Aligner object.}
        \item{index}{A list containing the information of the genom index.}
        \item{alignments}{A data.frame() containing the path to of the alignment filenames.}
        \item{paired}{logical(1). Samples in \sQuote{qProject} are paired-end reads.}
        \item{junction}{logical(1). Exon-exon-junction database is created for alignment.}
        \item{bisulfite}{logical(1). Samples in \sQuote{qProject} are bisultite reads.}
    }
}

\section{Functions and methods}{

    Constructor:
    \describe{
        \item{qProject}{Returns a \sQuote{qProject} object.}
    }
    
    Accessors:
    \describe{
        \item{show}{Compactly display the object.}
    }
    
    Methods:
    \describe{
        \item{qSaveProject}{Save the \sQuote{qProject} object, returning the file name of the saved R object.}    
        \item{qReadProject}{Restore the \sQuote{qProject} object and returns it.}    
    }

}

\author{Anita Lerch}

\seealso{
  \code{\link{qAlign}}, \code{\link{qCount}}
}

\examples{
\dontrun{
td <- tempdir()
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
aligner <- "Rbowtie"

project <- qProject(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

show(project)

project <- qAlign(project)

show(project)

projectFilename <- qSaveProject(project)
project <- qReadProject(projectFilename)
}
}

\keyword{classes}
