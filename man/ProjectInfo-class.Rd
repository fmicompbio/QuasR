\name{ProjectInfo}

\docType{class}
\alias{ProjectInfo-class}
\alias{ProjectInfo}

% accessors

% methods
\alias{alignHierarchical,ProjectInfo-method}
\alias{saveProjectInfo,ProjectInfo-method}
\alias{readProjectInfo,ProjectInfo-method}

\title{A class for meta-data for deep sequencing analysis.}
  
\description{
  Container class for meta-data needed during the deep sequencing analysis.
}

\usage{

## Constructor
ProjectInfo(sampleFile="Sample.txt",
     genome=".", 
     annotationFile="Annotation.txt", 
     aligner="Rbowtie", 
     projectName="projectinfo", 
     path=".", 
     lib=NULL)

## Accessors


## Methods
\S4method{show}{ProjectInfo}(object)
\S4method{alignHierarchical}{ProjectInfo}(object)
\S4method{saveProjectInfo}{ProjectInfo}(object)
readProjectInfo(filename)

}

\arguments{

  \item{sampleFile}{A character() vector of the sample file path name.}

  \item{genome}{A character() vector of the \sQuote{BSgenome} name or 
    of the directory path names fasta files path names or
    of a single fasta file path names.}

  \item{annotationFile}{A character() vector of the annotation file path name.}
  
  \item{Aligner}{A character() vector of the name of the alinger package used for the alignments.}

  \item{projectName}{A character() vector of the project name, used for saving the ProjectInfo.}

  \item{path}{A character() vector of the directory path names, at which the resulting files are stored.}

  \item{lib}{A character() vector describing the location of R library trees
          to search through, or ‘NULL’.  The default value of ‘NULL’
          corresponds to all libraries currently known.  Non-existent
          library trees are silently ignored.}

}

\section{Slots}{
  \describe{

  \item{name:}{A character() vector of the project name.}
  
  \item{id:}{A character() vector of a unique identifier for the project, used for saving the ProjectInfo.}
  
  \item{path:}{A character() vector of the directory path names, at which the resulting files are stored.}  

  \item{samples:}{A data.frame() containing the path to of the sequence filenames and the sample name.}

  \item{annotations:}{A data.frame() containing the path to of the annotation filenames and the feature name.} ##TODO rewrite
    
  \item{genome:}{The BSgenome object or a list containing information of a genomic fasta sequence.} ##TODO rewrite
 
  \item{aligner:}{The Aligner object.}

  \item{index:}{A list containing the information of the genom index.}

  \item{alignments:}{A data.frame() containing the path to of the alignment filenames.}   
  }
}

\section{Functions and methods}{

  Constructor:
  \describe{

	\item{ProjectInfo:}{Returns a \code{ProjectInfo} object.}

  }

  Accessors:
  \describe{
  }

  Methods:
  \describe{

	\item{show}{Compactly display the object.}
	
	\item{saveProjectInfo}{Save the \sQuote{ProjectInfo} object, returning the file name of the saved R object.}
	
	\item{readProjectInfo}{Restore the \sQuote{ProjectInfo} object, returning the \sQuote{ProjectInfo} object.}

  }

}

\author{Anita Lerch}

\examples{

  td <- tempdir()
  sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
  annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
  genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
  aligner <- "Rbowtie"

  project <- ProjectInfo(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

  show(project)

  project <- alignHierarchical(project)

  show(project)
  
  projectFilename <- saveProjectInfo(project)
  project <- readProjectInfo(projectFilename)

}

\keyword{classes}
