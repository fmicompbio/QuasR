\name{alignmentStats}
\alias{alignmentStats}
\title{Get statistics on alignments}
\description{
  Get statistics on alignments from bam file or \code{qProjec} object.
}
\usage{
alignmentStats(x, collapseBySample=TRUE)
}
\arguments{
  \item{x}{the source of alignment bam files, one of:
    \itemize{
      \item a \code{character} vector with bam files
      \item a \code{qProject} object
    }
  }
  \item{collapseBySample}{If \code{TRUE} and \code{x} is a \code{qProject}
    object, sum counts for bam files with identical sample names.}
}
\details{
  Internally, \code{alignmentStats} queries the bam index files similar
  to 'idxstats' from samtools.

  If \code{x} is a \code{qProject} object, the auxiliary bam files will
  not contain any unmapped reads, and the corresponding unmapped counts
  are calculated by substracting auxiliary mapped counts from the total
  reads. The latter correspond to the unmapped counts from the
  corresponding genome bam files.
}
\value{
  A \code{matrix} with one row per bam file and three columns
  ("seqlength", "mapped" and "unmapped").
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Anita Lerch, Dimos Gaidatzis and Michael Stadler
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[QuasR]{qProject}},
}
\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# using a qProject
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"
auxFile <- "extdata/auxiliaries.txt"

proj <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile)
alignmentStats(proj)

# using bam files
alignmentStats(alignments(proj)$genome$FileName)
alignmentStats(unlist(alignments(proj)$aux))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
