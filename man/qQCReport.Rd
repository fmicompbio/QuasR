\name{qQCReport}
\alias{qQCReport}

\title{QuasR Quality Control Report}

\description{
Generate quality control plots for a \code{qProject} object or a vector
of fasta/fastq/bam files. The available plot vary depending on the types
of available input (fasta, fastq or bam files).
}

\usage{
    qQCReport(input, pdfFilename=NULL, chunkSize=1e6L, clObj=NULL, \dots)
}

\arguments{
    \item{input}{A vector of files or a \code{qProject} object as
      returned by \code{qAlign}}    
    \item{pdfFilename}{The path and name of a pdf file to store the
      report. If \code{NULL}, the quality control plots will be
      generated in separate plotting windows on the standard graphical
      device.}
    \item{chunkSize}{The number of sequences, sequence pairs (for
      paired-end data) or alignments that will be sampled from each data
      file to collect quality statistics}
    \item{clObj}{a cluster object to be used for parallel processing of
      multiple input files.}
    \item{\dots}{additional arguments that will be passed to the functions
      generating the individual quality control plots, see \sQuote{Details}.}
}

\details{This function generates quality control plots for all input
  files or the sequence and alignment files contained in a
  \code{qProject} object, allowing assessment of the quality of a
  sequencing experiment. \code{qQCReport} uses functionality from the
  \pkg{ShortRead} package to collect quality data, and visualizes the
  results similarly as the \sQuote{FastQC} quality control tool from
  Simon Andrews (see \sQuote{References} below). It is recommended to
  create PDF reports (\code{pdfFilename} argument), for which the
  plot layouts have been optimised.
  
  Some plots will only be generated if the necessary information is
  available (e.g. base qualities in fastq sequence files).

  The currently available plot types are:
  \describe{
    \item{\emph{Quality score boxplot}}{shows the distribution of
      base quality values as a box plot for each position in the input
      sequence. The background color (green, orange or red)
      indicates ranges of high, intermediate and low qualtities. The
      plot is available for fastq only (bam files may contain base
      quality information, which is however not used here because reads
      contained in the bam file, e.g. aligned reads, may not be a
      representative subsample of all sequenced reads).}
    \item{\emph{Nucleotide frequency}}{plot shows the frequency of A, C,
      G, T and N bases by position in the read. The plot is always
      available.}
    \item{\emph{Duplication level}}{plot shows for each sample the
      fraction of reads observed at different duplication levels
      (e.g. once, two-times, three-times, etc.). In addition, the most
      frequent sequences are listed. The plot is available for fasta or
      fastq files, but not for bam files, again because contained
      reads may not be representative for the experiment.}
    \item{\emph{Mapping statistics}}{shows fractions of reads that were
      (un)mappable to the reference genome. This plot is available for
      bam input.}
    \item{\emph{Library complexity}}{shows fractions of unique
      read(-pair) alignment positions, as a measure of the complexity in
      the sequencing library. Please note that this measure is not
      independent from the total number of reads in a library, and is best
      compared between libraries of similar sizes. This plot is
      available for bam input.}
    \item{\emph{Mismatch frequency}}{shows the frequency and position
      (relative to the read sequence) of mismatches in the alignments
      against the reference genome. The plot is available for bam input.}
    \item{\emph{Mismatch types}}{shows the frequency of
      read bases that caused mismatches in the alignments to the
      reference genome, separately for each genome base. This plot is
      available for bam input.}
    \item{\emph{Fragment size}}{shows the distribution of fragment sizes
      inferred from aligned read pairs. This plot is available for
      paired-end bam input.}
  }

  One approach to assess the quality of a sample is to compare its
  control plots to the ones from other samples and search for relative
  differences. Special quality measures are expected for certain types
  of experiments: A genomic re-sequencing sample with an
  overrepresentation of T bases may be suspicious, while such a
  nucleotide bias is normal for a directed bisulfite-sequencing sample.

  Additional arguments can be passed to the internal functions that
  generate the individual quality control plots using \code{\dots{}}:
  \describe{
    \item{\code{lmat}:}{a matrix (e.g. \code{matrix(1:12, ncol=2)}) used
      by an internal call to the \code{layout} function to specify the
      positioning of multiple plot panels on a device page. Individual panels
      correspond to different samples.}
    \item{\code{breaks}:}{a numerical vector
      (e.g. \code{c(1:10)}) defining the bins used by
      the \sQuote{Duplication level} plot.}
  }
}
  
\value{The function is called for its side effect of generating quality
  control plots. It invisibly returns a list with components that
  contain the data used to generate each of the QC plots. Available
  components are (depending on input data, see \sQuote{Details} above):
  \itemize{
    \item{\emph{qualByCycle}:}{ quality score boxplot}
    \item{\emph{nuclByCycle}:}{ nucleotide frequency plot}
    \item{\emph{duplicated}:}{ duplication level plot}
    \item{\emph{mappings}:}{ mapping statistics barplot}
    \item{\emph{uniqueness}:}{ library complexity barplot}
    \item{\emph{errorsByCycle}:}{ mismatch frequency plot}
    \item{\emph{mismatchTypes}:}{ mismatch type plot}
    \item{\emph{fragDistribution}:}{ fragment size distribution plot}
  }
}

\references{FastQC quality control tool at \url{http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/}}

\author{
Anita Lerch, Dimos Gaidatzis and Michael Stadler
}

\seealso{
  \code{\linkS4class{qProject}}, \code{\link{qAlign}},
  \code{\link[=ShortReadBase-package]{ShortRead}} package
}

\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create alignments
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"

proj <- qAlign(sampleFile, genomeFile)

# create quality control report
qQCReport(proj, pdfFilename="qc_report.pdf")
}

\keyword{methods}
