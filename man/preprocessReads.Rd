\name{preprocessReads}
\alias{preprocessReads}
\title{Preprocess Short Reads}
\description{
Truncate sequences, remove parts matching to adapters and filter out low
quality or low complexity sequences from (compressed) 'fasta' or 'fastq' files.
}
\usage{
preprocessReads(filename, outputFilename=NULL,
             filenameMate=NULL, outputFilenameMate=NULL,
             truncateStartBases=NULL, truncateEndBases=NULL, 
             Lpattern="", Rpattern="",
             max.Lmismatch=rep(0:2, c(6,3,100)), max.Rmismatch=rep(0:2, c(6,3,100)),
             with.Lindels=FALSE, with.Rindels=FALSE,
             minLength=14L, nBases=2L, complexity=NULL,
             nrec=1000000L, clObj=NULL)
}
\arguments{
    \item{filename}{the name(s) of the input sequence file(s).}
    \item{outputFilename}{the name(s) of the output sequence file(s).}
    \item{filenameMate}{for paired-end experiments, the name(s) of the
      input sequence file(s) containing the second read (mate) of each pair.}
    \item{outputFilenameMate}{for paired-end experiments, the name(s) of the
      output sequence file(s) containing the second read (mate) of each pair.}
    \item{truncateStartBases}{integer(1): the number of bases to be
      truncated (removed) from the begining of each sequence.}
    \item{truncateEndBases}{integer(1): the number of bases to be
      truncatet (removed) from the end of each sequence.}
    \item{Lpattern}{character(1): the left (5'-end) adapter sequence.}
    \item{Rpattern}{character(1): the right (3'-end) adapter sequence.}
    \item{max.Lmismatch}{mismatch tolerance when searching for matches
      of \code{Lpattern} (see \sQuote{Details}).}
    \item{max.Rmismatch}{mismatch tolerance when searching for matches
      of \code{Rpattern} (see \sQuote{Details}).}
    \item{with.Lindels}{if \code{TRUE}, indels are allowed in the
      alignments of the suffixes of \code{Lpattern} with the subject,
      at its beginning (see \sQuote{Details}).}
    \item{with.Rindels}{same as \code{with.Lindels} but for alignments
      of the prefixes of \code{Rpattern} with the subject, at its end
      (see \sQuote{Details}).}
    \item{minLength}{integer(1): the minimal allowed sequence length.}
    \item{nBases}{integer(1): the maximal number of Ns allowed per sequence.}
    \item{complexity}{\code{NULL} (default) or numeric(1): If not \code{NULL},
      the minimal sequence complexity, as a fraction of the average complexity
      in the human genome (~3.9bits). For example, \code{complexity = 0.5} will
      filter out sequences that do not have at least half the complexity of the
      human genome. See \sQuote{Details} on how the complexity is calculated.}
    \item{nrec}{integer(1): the number of sequence records to read at a
      time.}
    \item{clObj}{a cluster object to be used for parallel processing of
      multiple files (see \sQuote{Details}).}
}
\details{
Sequence files can be in fasta or fastq format, and can be compressed by
either gzip, bzip2 or xz (extensions .gz, .bz2 or .xz). Multiple files
can be processed by a single call to \code{preprocessReads}; in that
case all sequence file vectors must have identical lengths.

\code{nrec} can be used to limit the memory usage when processing
large input files. \code{preprocessReads} iteratively loads chunks of
\code{nrec} sequences from the input until all data been processed.

Sequence pairs from paired-end experiments can be processed by
specifying pairs of input and output files (\code{filenameMate} and
\code{outputFilenameMate} arguments). In that case, it is assumed that
pairs appear in the same order in the two input files, and only pairs
in which both reads pass all filtering criteria are written to the
output files, maintaining the consistent ordering.

If output files are compressed, the processed sequences are first
written to temporary files (created in the same directory as the final
output file), and the output files are generated at the end by compressing
the temporary files.

For the trimming of left and/or right flanking sequences (adapters) from
sequence reads, the \code{\link[Biostrings]{trimLRPatterns}} function
from package \pkg{Biostrings} is used, and the arguments \code{Lpattern},
\code{Rpattern}, \code{max.Lmismatch}, \code{max.Rmismatch},
\code{with.Lindels} and \code{with.Rindels} are used in the call to
\code{trimLRPatterns}. \code{Lfixed} and \code{Rfixed} arguments
of \code{trimLRPatterns} are set to \code{TRUE}, thus only fixed
patterns (without IUPAC codes for ambigous bases) can be
used. Currently, trimming of adapters is only supported for single read
experiments.

Sequence complexity (\eqn{H}) is calculated based on the dinucleotide
composition using the formula (Shannon entropy): \deqn{H = -\sum_i {f_i \log_2 f_i}},
where \eqn{f_i} is the fraction of dinucleotide \eqn{i} from all
dinucleotides in the sequence. Sequence reads that fulfill the condition
\eqn{H/H_r \ge c} are retained (not filtered out), where \eqn{H_r =
  3.908} is the reference complexity in bits obtained from the human
genome, and \eqn{c} is the value given to the argument \code{complexity}.

If an object that inherits from class \code{cluster} is provided to
the \code{clObj} argument, for example an object returned by
\code{\link[parallel]{makeCluster}} from package \pkg{parallel},
multiple files will be processed in parallel using
\code{\link[parallel]{clusterMap}} from package \pkg{parallel}.
}
\value{
A matrix with summary statistics on the processed
sequences, containing:
\itemize{
  \item One column per input file (or pair of input files for paired-end
  experiments).
  \item The number of sequences or sequence pairs in rows:
  \itemize{
    \item \code{totalSequences} - the total number in the input
    \item \code{matchTo5pAdaptor} - matching to \code{Lpattern}
    \item \code{matchTo3pAdaptor} - matching to \code{Rpattern}
    \item \code{tooShort} - shorter than \code{minLength}
    \item \code{tooManyN} - more than \code{nBases} Ns
    \item \code{lowComplexity} - relative complexity below
    \code{complexity}
    \item \code{totalPassed} - the number of sequences/sequence pairs
    that pass all filtering criteria and were written to the output
    file(s).
  }
}
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
  \code{\link[Biostrings]{trimLRPatterns}} from pacakge \pkg{Biostrings},
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
}
\examples{
# sample files
infiles <- system.file(package="QuasR", "extdata",
                       c("rna_1_1.fq.bz2","rna_1_2.fq.bz2"))
outfiles <- paste(tempfile(pattern=c("output_1_","output_2_")),".fastq",sep="")

# single read example
preprocessReads(infiles, outfiles, nBases=0, complexity=0.6)
unlink(outfiles)

# paired-end example
preprocessReads(filename=infiles[1],
                outputFilename=outfiles[1],
                filenameMate=infiles[2],
                outputFilenameMate=outfiles[2],
                nBases=0, complexity=0.6)
unlink(outfiles)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
