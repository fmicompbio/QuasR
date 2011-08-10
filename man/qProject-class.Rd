\name{qProject}

\docType{class}
\alias{qProject-class}
\alias{qProject}
\alias{show,qProject-method}
\alias{qSaveProject,qProject-method}
\alias{qSaveProject}
\alias{qReadProject}

\title{A class for meta-data for deep sequencing analysis}
  
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
     lib.loc=NULL)

## Accessors


## Methods
\S4method{show}{qProject}(object)
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
}

\section{Slots}{
    \describe{
    
        \item{name}{A character() vector of the project name.}        
        \item{id}{A character() vector of a unique identifier for the project, used for saving the \sQuote{qProject}.}        
        \item{path}{A character() vector of the directory path names, at which the resulting files are stored.}      
        \item{samples}{A data.frame() containing the path to of the sequence filenames and the sample name.}
        \item{annotations}{A data.frame() containing the path to of the annotation filenames and the feature name. TODO rewrite}
        \item{genome}{The BSgenome object or a list containing information of a genomic fasta sequence. TODO rewrite}
        \item{aligner}{The Aligner object.}
        \item{index}{A list containing the information of the genom index.}
        \item{alignments}{A data.frame() containing the path to of the alignment filenames.}
        \item{paired}{A logical() TODO.}
        \item{junction}{A logical() TODO.}
        \item{stranded}{A logical() TODO.}
        \item{bisulfite}{A logical() TODO.}
    }
}

\section{Functions and methods}{

    Constructor:
    \describe{
        \item{qProject}{Returns a \sQuote{qProject} object.}
    }
    
    Accessors:
    %\describe{\item{}{}}
    
    Methods:
    \describe{
        \item{show}{Compactly display the object.}    
        \item{qSaveProject}{Save the \sQuote{qProject} object, returning the file name of the saved R object.}    
        \item{qReadProject}{Restore the \sQuote{qProject} object and returns it.}    
    }

}

\author{Anita Lerch}

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
