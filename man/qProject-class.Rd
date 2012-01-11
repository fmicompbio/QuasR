\name{qProject}

\docType{class}
\alias{qProject-class}
\alias{qProject}
\alias{show,qProject-method}
\alias{qSaveProject,qProject-method}
\alias{qSaveProject}
\alias{qReadProject}

\alias{getGenome,qProject-method}
\alias{getGenome}
\alias{getSamples,qProject-method}
\alias{getSamples}
\alias{getAlignments,qProject-method}
\alias{getAlignments}
\alias{paired,qProject-method}
\alias{paired}


\title{A meta-data class for deep sequencing analysis}
  
\description{
Container class for meta-data needed during the deep sequencing analysis.
}

\usage{

## Constructor
qProject(sampleFile,
     genome, 
     auxiliaryFile=NULL, 
     aligner="Rbowtie", 
     projectName="qProject", 
     bamfileDir=NULL, 
     splicedAlignment=FALSE, 
     bisulfite=FALSE,
     lib.loc=NULL,
     indexLocation=NULL,
     maxHits=100L,
     cacheDir=NULL,
     alignmentParameter)


## Accessors
\S4method{show}{qProject}(object)
%\S4method{getGenome}{qProject}(object)
%\S4method{getSamples}{qProject}(object)
%\S4method{getAlignments}{qProject}(object)
%\S4method{paired}{qProject}(object)

## Methods
\S4method{qSaveProject}{qProject}(object)
qReadProject(filename)

}

\arguments{
    \item{sampleFile}{A character() vector of the sample file path name.}
    \item{genome}{A character() vector of the \sQuote{BSgenome} name or of the directory path names fasta files path names or of a single fasta file path names.}
    \item{auxiliaryFile}{A character() vector of the auxiliary file path name.}  
    \item{aligner}{A character() vector of the name of the alinger package used for the alignments.}
    \item{projectName}{A character() vector of the project name, used for saving the \sQuote{qProject}.}
    \item{bamfileDir}{A character() vector of the directory path names, at which the resulting bam files will be stored. If it is NULL, then the bam files are stored at the same location as the fasta/q files.}
    \item{indexLocation}{A character() vector of the directory path names, at which the index files are stored.}    
    \item{lib.loc}{A character() vector describing the location of R library trees to search through, or \sQuote{NULL}.  The default value of \sQuote{NULL} corresponds to all libraries currently known.  Non-existent library trees are silently ignored.}
    \item{object}{A qProject() object.}
    \item{filename}{A character() vector of the \sQuote{qProject} rds-file.}
    \item{splicedAlignment}{A logical(1). If \code{TRUE} then the samples are additionaly aligned to the Exon-exon-junction database. Up to now not implemented.}
    \item{bisulfite}{logical(1). If \code{TRUE} then the samples in the \sQuote{qProject} are handled as bisultite reads. Up to now not implemented.}
    \item{cacheDir}{A character() vector of the directory path names, where the temporary files are saved.}
    \item{maxHits}{A integer(1) vector of the maximal allowed hits in the alignment file.}
    \item{alignmentParameter}{A character(1) vector of the alignment parameter (e.g. q,f,v,k,a,m,M,best,strata) without the index, sample and output location/parameter. The Parameter \code{maxHits} will be ignored. For example: \dQuote{bowtie -v 1 -a -f}}
}

\section{Slots}{
    \describe{
        \item{name}{A character() vector of the project name.}
        \item{id}{A character() vector of a unique identifier for the project, used for saving the \sQuote{qProject}.}
        \item{env}{An environment containing the variable used by QuasR, e.g. samples, auxiliaries, genome, aligner, alignment.}
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
        \item{getGenome}{Get the genome object.}
        \item{getSamples}{Get the samples object.}
        \item{getAlignments}{Get the alignments object.}
        \item{paired}{Return 'TRUE' if the project is a paired-end reads project.}
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

sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"

project <- qProject(sampleFile, genome=genomeName, auxiliaryFile=auxiliaryFile, aligner="Rbowtie", bamfileDir=tempdir())

show(project)

project <- qAlign(project)

show(project)

regions <- GRanges(seqnames="chr2L", ranges=IRanges(1, width=1000000))
count <- qCount(project, regions)

projectFilename <- qSaveProject(project)
project <- qReadProject(projectFilename)
}
}

\keyword{classes}
