\name{preprocessReads}
\alias{preprocessReads}

\title{Preprocess Short Reads}

\description{
Remove adaptors, truncate bases and filter short reads from 'fasta' or 'fastq' files.
}

\usage{
    preprocessReads(filename, filenameMate=NULL,
                    outputFilename=NULL, outputFilenameMate=NULL,
                    Lpattern="", Rpattern="",
                    max.Lmismatch=rep(0:2, c(6,3,100)), max.Rmismatch=rep(0:2, c(6,3,100)),
                    with.Lindels=FALSE, with.Rindels=FALSE,
                    truncateStartBases=NULL, truncateEndBases=NULL, 
                    minLength=14L, nBases=2L, complexity=NULL,
                    nrec=1000000L
                    )
}

\arguments{
    \item{filename}{A \sQuote{character} vector containing the fasta/q filenames}
    \item{filenameMate}{A \sQuote{character} vector containing the mates of the fasta/q filenames}
    \item{outputFilename}{A \sQuote{character} vector containing the output fasta/q filenames}
    \item{outputFilenameMate}{A \sQuote{character} vector containing the output mates of the fasta/q filenames}
    \item{Lpattern}{A character(1) containing the left adaptor pattern.}
    \item{Rpattern}{A character(1) containing the right adaptor pattern.}
    \item{max.Lmismatch}{, for more information xxx}
    \item{max.Rmismatch}{, for more information xxx}
    \item{with.Lindels}{, for more information xxx}
    \item{with.Rindels}{, for more information xxx}
    \item{truncateStartBases}{A integer(1) vector containing the number of bases which should be truncatet at the begining of the sequence.}
    \item{truncateEndBases}{A integer(1) vector containing the number of bases which should be truncatet at the end of the sequence.}
    \item{minLength}{A integer(1) vector representing the minimal allowed length of the sequences.}
    \item{nBases}{A integer(1) vector representing the maximal allowed Ns in the sequences.}
    \item{complexity}{A numeric(1) vector containig the treshold of the minimal allowd entropie of the sequences.}
    \item{nrec}{Single integer. The maximum of number of records to read in per chunk. TODO check what happen with: Negative values are ignored.}
}

\details{
 TODO
}

\value{
   TODO
}

\author{Anita Lerch}

}

\examples{
\dontrun{
#TODO

td <- tempdir()
sampleFile <- system.file(package="QuasR", "extdata", "samples.txt")
project <- preprocessReads(c("",""))
}  
}

\keyword{methods}
\keyword{manip}
