\name{QuasR-package}
\alias{QuasR-package}
\alias{QuasR}
\docType{package}

\title{Quantify and Annotate Short Reads in R}

\description{
This package provides a pipeline for the quantification and analysis of Short Reads.
}

\details{
See packageDescription('QuasR') for package details.
}

\author{
Anita Lerch, Dimos Gaidatzis and Michael Stadler
}

%\references{TODO}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
  \code{\link{qAlign}}, \code{\link{qCount}}, \code{\link{qProfile}}, \code{\link{qMeth}}, \code{\link{qQCReport}}
}

\examples{
library(GenomicRanges)

# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)


# create alignments
sampleFile <- "extdata/samples_chip_single.txt"
genomeFile <- "extdata/hg19sub.fa"
auxFile <- "extdata/auxiliaries.txt"

proj <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile)


# get alignment statistics
alignmentStats(proj)


# quantify alignments
gr <- GRanges(seq=paste("chr",1:3,sep=""), IRanges(start=1000*1:3, width=500))
names(gr) <- paste(seqnames(gr), start(gr), sep=":")
gr

qCount(proj, gr)
}
