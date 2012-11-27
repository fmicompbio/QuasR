\name{qExportWig}
\alias{qExportWig}

\title{QuasR wig file export}

\description{
Create a fixed-step wig file from the alignments in the genomic bam
files of the \sQuote{QuasR} project.
}

\usage{
qExportWig(proj, file=NULL, collapseBySample=TRUE, binsize=100L,
           shift=0L, strand=c("*","+","-"), scaling=TRUE,
           tracknames=NULL, log2p1=FALSE,
           colors=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                    "#66A61E", "#E6AB02", "#A6761D", "#666666"))
}

\arguments{
  \item{proj}{A \code{qProject} object as returned by \code{qAlign}.}
  \item{file}{A character vector with the name(s) for the wig file(s) to
    be generated. Either \code{NULL} or a vector of the same length as
    the number of bam files (for \code{collapseBySample=FALSE}) or the
    number of unique sample names (for \code{collapseBySample=TRUE}) in
    \code{proj}. If \code{NULL}, the wig file names are generated from
    the names of the genomic bam files or unique sample names with an
    added \dQuote{.wig.gz} extension.}
  \item{collapseBySample}{If \code{TRUE}, genomic bam files with identical
    sample name will be combined (summed) into a single track.}
  \item{binsize}{a numerical value defining the bin and step size for the wig
    file(s). \code{binsize} will be coerced to \code{integer()}.}
  \item{shift}{Either a vector or a scalar value defining the read shift (e.g.
    half of fragment length, see \sQuote{Details}). If \code{length(shift)>1},
    the length must match the number of bam files in \sQuote{proj}, and
    the i-th sample will be converted to wig using the value in
    \code{shift[i]}. \code{shift} will be coerced to \code{integer()}. For
    paired-end alignments, \code{shift} will be ignored, and a warning
    will be issued if it is set to a non-zero value (see \sQuote{Details}).}
  \item{strand}{Only count alignments of \code{strand}. The default (\dQuote{*})
    will count all alignments.}
  \item{scaling}{If TRUE, the output values in the wig file(s)
    will be linearly scaled by the total number of aligned reads per sample
    to improve comparability (see \sQuote{Details}).}
  \item{tracknames}{A character vector with the names of the
    tracks to appear in the track header. If \code{NULL}, the sample
    names in \code{proj} will be used.}
  \item{log2p1}{If \code{TRUE}, the number of alignments \code{x} per bin will
    be transformed using the formula \code{log2(x+1)}.}
  \item{colors}{A character vector with R color names to be used for the tracks.}
}

\details{
  \code{qExportWig()} uses the genome bam files in \code{proj} as input
  to create wig files with the number of alignments (pairs) per window of
  \code{binsize} nucleotides. By default (\code{collapseBySample=TRUE}),
  one file per unique sample will be created. If \code{collapseBySample=FALSE},
  one file per genomic bam file will be created.
  See \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html} for the
  definition of the wig format.
  
  The genome is tiled with sequential windows of length \code{binsize},
  and alignments in the bam file are assigned to these windows: Single
  read alignments are assigned according to their 5'-end coordinate
  shifted by \code{shift} towards the 3'-end (assuming that the 5'-end
  is the leftmost coordinate for plus-strand alignments, and the rightmost
  coordinate for minus-strand alignments). Paired-end alignments are
  assigned according to the base in the middle between the leftmost and
  rightmost coordinates of the aligned pair of reads. Each pair of reads
  is only counted once, and not properly paired alignments are ignored.

  For \code{scaling=TRUE}, the number of alignments per bin \sQuote{n}
  for the sample \sQuote{i} are linearly scaled to the minimum total
  number of alignments over all samples in \code{proj} according to:
  \sQuote{n_s = n /N[i] *min(N)} where \sQuote{n_s} is the scaled number
  of alignments in the bin and \sQuote{N} is a vector with the total
  number of alignments for each sample.

  If \code{file} ends with \sQuote{.gz}, the resulting wig file will be
  compressed using gzip and is suitable for uploading as a custom track
  to your favortie genome browser (e.g. UCSC or Ensembl).
}

\value{
  (invisible) The file name of the generated wig file(s).
}
\author{
Anita Lerch, Dimos Gaidatzis and Michael Stadler
}

\seealso{
  \code{\linkS4class{qProject}}, \code{\link{qAlign}}
}

\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create alignments
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"
proj <- qAlign(sampleFile, genomeFile)

# export wiggle file
qExportWig(proj, binsize=100L, shift=0L, scaling=TRUE)
}

\keyword{utilites}
