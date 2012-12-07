\name{qProject-class}
\docType{class}

% Class:
\alias{class:qProject}
\alias{qProject-class}
\alias{qProject}

% Constructors:
%\alias{qProject}
%\alias{updateObject,qProject-method}

% Coercion:
%\alias{as.data.frame,qProject-method}

% Accessors:
\alias{length,qProject-method}
%\alias{sampleNames}
%\alias{sampleNames<-}
%\alias{sampleNames,qProject-method}
%\alias{sampleNames<-,qProject-method}
\alias{genome,qProject-method}
\alias{auxiliaries}
\alias{auxiliaries,qProject-method}
\alias{alignments}
\alias{alignments,qProject-method}
\alias{[,qProject,ANY,ANY-method}

% "show" method:
\alias{show,qProject-method}

\title{qProject objects}

\description{
  The qProject class is a container for the meta-data (e.g. sample
  names, paths and names of sequence and alignment files) associated
  with a high-throughput sequencing experiment analyzed with \code{QuasR}.
}

\details{
  The qProject class is returned by \link{qAlign} and stores all
  information on a high-throughput sequencing experiment analyzed with
  \code{QuasR}. qProject objects can be conveniently passed to
  \sQuote{q}-functions (function name starting with the letter
  \sQuote{q}). The information is stored in the following slots:
  \describe{
    \item{\code{reads}}{a 'data.frame' with sequence read files.}
    \item{\code{reads_md5subsum}}{a 'data.frame' with fingerprints for
      sequence read files.}
    \item{\code{alignments}}{a 'data.frame' with alignment files.}
    \item{\code{samplesFormat}}{a 'character(1)' specifying the format
      of input files.}
    \item{\code{genome}}{a 'character(1)' specifying the reference genome.}
    \item{\code{genomeFormat}}{a 'character(1)' specifying the format of
      the reference genome.}
    \item{\code{aux}}{a 'data.frame' with auxiliary reference sequence files.}
    \item{\code{auxAlignments}}{a 'data.frame' with alignment files for
      auxiliary reference sequence files.}
    \item{\code{aligner}}{a 'character(1)' specifying the aligner.}
    \item{\code{maxHits}}{a 'numeric(1)' specifying the maximum number
      of alignments per sequence.}
    \item{\code{paired}}{a 'character(1)' specifying the paired-type;
      one of "no", "fr", "rf", "ff".}
    \item{\code{splicedAlignment}}{a 'logical(1)'; \code{TRUE} when
      performing spliced-alignments.}
    \item{\code{snpFile}}{a 'character(1)' with a file name containing
      SNP information.}
    \item{\code{bisulfite}}{a 'character(1)' defining the bisulfite
  type; one of "no", "dir", "undir".}
    \item{\code{alignmentParameter}}{a 'character(1)' with aligner
      command line parameters.}
    \item{\code{projectName}}{a 'character(1)' with the project name.}
    \item{\code{alignmentsDir}}{a 'character(1)' with the directory to
      be used to store alignment files.}
    \item{\code{lib.loc}}{a 'character(1)' with the library directory to
      use for installing of alignment index packages.}
    \item{\code{cacheDir}}{a 'character(1)' with a directory to use for
      temporary files.}
    \item{\code{alnModeID}}{a 'character(1)' used internally to indicate
      the alignment mode.}
  }
}

%\section{Constructor}{
%  \describe{
%    \item{}{
%      \code{qProject()}:
%      Creates a qProject object.
%      \describe{
%        \item{\code{reads}}{...}
%      }
%    }
%  }
%}

%\section{Coercion}{
%  In the code snippets below, \code{x} is a qProject object.
%
%  \describe{
%    \item{}{
%      \code{as.data.frame(x)}:
%      Creates a data.frame with the information stored in the qProject object.
%    }
%  }
%}

\section{Accessors}{
  In the following code snippets, \code{x} is a qProject object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of input files.
    }
    %    \item{}{
    %      \code{sampleNames(x)}, \code{samplesNames(x) <- value}:
    %      Gets or sets the sample names associated to input files.
    %      \code{value} can be a character vector or a factor.
    %    }
    \item{}{
      \code{genome(x)}:
      Gets the reference genome as a 'character(1)'. The type of genome
      is stored as an attribute in
      \code{attr(genome(x),"genomeFormat")}: "BSgenome" indicates that
      \code{genome(x)} refers to the name of a BSgenome package, "file"
      indicates that it contains the path and filename of a genome in
      FASTA format.
    }
    \item{}{
      \code{auxiliaries(x)}:
      Gets a \code{data.frame} with auxiliary target sequences, with one
      row per auxiliary target, and columns "FileName" and "AuxName".
    }
    \item{}{
      \code{alignments(x)}:
      Gets a list with two elements "genome" and
      "aux". \code{alignments(x)$genome} contains a \code{data.frame}
      with \code{length(x)} rows and the columns "FileName" (containing
      the path to bam files with genomic alignments) and
      "SampleName". \code{alignments(x)$aux} contains a
      \code{data.frame} with one row per auxiliary target sequence (with
      auxiliary names as row names), and \code{length(x)} columns.
    }
  }
}

\section{Subsetting}{
  In the following code snippets, \code{x} is a qProject object.

  \describe{
    \item{}{
      \code{x[i]}:
      Get \code{qProject} object instance with \code{i} input files,
      where \code{i} can be an NA-free logical, numeric, or character vector.
    }
  }
}

\author{
Anita Lerch, Dimos Gaidatzis and Michael Stadler
}

\seealso{
  \link{qAlign}
}

\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create alignments
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"
auxFile <- "extdata/auxiliaries.txt"

proj <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile)
proj

# alignment statistics using a qProject
alignmentStats(proj)

# alignment statistics using bam files
alignmentStats(alignments(proj)$genome$FileName)
alignmentStats(unlist(alignments(proj)$aux))
}
