\name{qQCReport}
\alias{qQCReport}

\title{QuasR Quality Control Report}

\description{
Generate quality control plots for a \sQuote{qProject} object. The
available plot vary depending on whether base qualities are available
(fastq rather than fasta files) and whether alignments have been made already.
}

\usage{
    qQCReport(qproject, pdfFilename=NULL, chunkSize=1e6, \dots)
}

\arguments{
    \item{qproject}{A \sQuote{qProject} object instance.}    
    \item{pdfFilename}{The path and name of a pdf file to store the
      report. If \sQuote{NULL}, the quality control plots will be
      generated in separate plotting windows on the standard graphical
      device.}
    \item{chunkSize}{The number of sequences or alignments that will be
      sampled from each file to collect quality statistics.}
    \item{\dots}{additional arguments that will be passed to the functions
      generating the individual quality control plots, see Details.}
}

\details{This function generates quality control plots for all samples
  in a \sQuote{qProject} object, allowing assessment of the qualty of a
  sequencing experiment. \sQuote{qQCReport} uses functionality from the
  \pkg{ShortRead} package to collect quality data, and visualizes the
  results similarly as the FastQC quality control tool from Simon
  Andrews.

  Some types of plots will only be made if the input sequences contain
  information on base qualities (i.e. are in fastq format), or if the
  input sequences have already been aligned using \sQuote{qAlign}.

  The following types of plots are available:
  \describe{
    \item{\emph{Quality score boxplot}}{shows the distribution of
      base quality values as a box plot for each position in the input
      sequence. The background color (green, orange or red)
      indicates ranges of high, intermediate and low qualtities. The
      plot is available for fastq and external bam files.}
    \item{\emph{Nucleotide frequency}}{plot shows the frequency of A, C,
      G, T and N bases by position in the input sequence.}
    \item{\emph{Duplication level}}{plot shows for each sample the
      fraction of reads observed at different duplication levels
      (e.g. once, two-times, three-times, etc.). In addition, the most
      frequent sequences are listed.}
    \item{\emph{Mapping statistics}}{shows the number and fraction of
      reads that were unmappable, or that were mapped with different
      numbers of hits to the reference genome are shown. This plot is
      available for projects with alignments created by the
      \sQuote{qAlign} function.}
    \item{\emph{Mismatch frequency}}{shows the frequency and position
      (relative to the short read) of mismatches in the alignments
      against the reference genome. This plot is available for projects
      with alignments created by the \sQuote{qAlign} function.}
  }

  Additional arguments can be passed to the internal functions that
  generate the individual quality control plots using \sQuote{\dots{}}:
  \describe{
    \item{\code{lmat}:}{a matrix (e.g. \code{matrix(1:18, ncol=3)}) used
      by an internal call to the \sQuote{layout} function to specify the
      ositioning of multiple panels on a device page. Individual panels
      correspond to different samples.}
    \item{\code{breaks}:}{a numerical vector
      (e.g. \code{c(0,1,5,10,50,100)}) defining the bins used by
      the \sQuote{Mapping statistics} and \sQuote{Duplication level} plots.}
  }
}

\value{No return value; the function is called for its side effect of
  generating quality control plots.}

\references{FastQC quality control tool at \url{http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/}}

\author{Anita Lerch, Michael Stadler}

\seealso{
  \code{\linkS4class{qProject}}, \code{\link{qAlign}},
  \code{\link[=ShortReadBase-package]{ShortRead}} package
}

\examples{
\dontrun{
td <- tempdir()
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
annotationFile <- system.file(package="QuasR", "extdata", "annotations.txt")
genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
aligner <- "Rbowtie"
project <- qProject(sampleFile, genome=genomeName, annotationFile=annotationFile, aligner=aligner, path=td)

project <- qAlign(project)

project <- qQCReport(project)
}  
}

\keyword{methods}
