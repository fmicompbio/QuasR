\name{Aligner}

\docType{class}
\alias{Aligner-class}
\alias{Aligner}

% accessors

% methods
\alias{execute,Aligner-method}

\title{Sequence aligner wrapper class.}
  
\description{
  A wrapper class for sequence aligner packages.
}

\usage{

## Constructor
Aligner(aligner=character())

## Accessors


## Methods
\S4method{show}{Aligner}(object)
\S4method{execute}{Aligner}(object, callstr)
\S4method{execute}{character}(object, callstr)
help{Aligner}(object) 
}

\arguments{

  \item{object}{A character() vector of the aligner package name.}
  
  \item{callstr}{A character() vector containing the command line execution command, see the \code{help} for more information.}

}

\section{Slots}{
  \describe{

  \item{pkgname:}{A character() vector of the aligner package name.}
  
  \item{pkgversion:}{A character() vector of the aligner package version.}
  
  }
}

\section{Functions and methods}{

  Constructor:
  \describe{

	\item{Aligner:}{Returns a \code{Aligner} object.}

  }

  Accessors:
  \describe{
  }

  Methods:
  \describe{

	\item{show}{Compactly display the object.}
	
	\item{execute}{Execute the command line execution command of the aligner, returning a character vector giving the output of the
     command.}
	
	\item{help}{Displays the help page of the aligner.}

  }

}

\author{Anita Lerch}

\examples{

  ## Example perfomed with the bowtie aligner R wrapper package 
  aligner <- Aligner("Rbowtie")
  show(aligner)
  help(aligner)
  
  td <- tempdir()
  
  ## Building a bowtie index
  refs <- dir(system.file(package="Rbowtie", "samples", "refs"), full=TRUE)
  indexdir <- file.path(td, "index")
  execute(alinger, sprintf("bowtie-build %s %s", refs, indexdir)

  ## Alignments
  reads <- system.file(package="Rbowtie", "samples", "reads", "reads.fastq")
  samfile <- tempfile()
  execute(alinger, sprintf("bowtie %s -f %s -S %s", refs, reads, samfile)  

}

\keyword{classes}
