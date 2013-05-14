\name{qMeth}
\alias{qMeth}
\title{Quantify DNA methylation}
\description{
Quantify methylation of cytosines from bisulfite sequencing data.
}
\usage{
qMeth(proj, query=NULL, reportLevel=c("C","alignment"),
      mode=c("CpGcomb","CpG","allC","var"),
      collapseBySample=TRUE, collapseByQueryRegion=FALSE,
      asGRanges=TRUE, mask=NULL, reference="genome",
      keepZero=TRUE, clObj=NULL)
}
\arguments{
  \item{proj}{a \code{qProject} object from a bisulfite sequencing experiment}
  \item{query}{a \code{GRanges} object with the regions to be
    quantified. If \code{NULL}, all available target sequences (e.g. the
    whole genome) will be analyzed. Available target sequences are
    extracted from the header of the first bam file.}
  \item{reportLevel}{report results combined for C's
    (\code{reportLevel}=\dQuote{C}}, the default) or individually for single
    alignments (\code{reportLevel}=\dQuote{alignment}). The latter imposes
    further restrictions on some arguments (see \sQuote{Details}).
  \item{mode}{cytosine quantification mode, one of:
    \itemize{
      \item \code{CpGcomb} : only C's in CpG context (strands combined)
      \item \code{CpG} : only C's in CpG context (strands separate)
      \item \code{allC} : all C's (strands separate)
      \item \code{var} : variant detection (all C's, strands separate)
    }
    \code{CpGcomb} is the default.}
  \item{collapseBySample}{if \code{TRUE}, combine (sum) counts from
    bamfiles with the same sample name.}
  \item{collapseByQueryRegion}{if \code{TRUE}, combine (sum) counts for
    all cytosines contained in the same query region.}
  \item{asGRanges}{if \code{TRUE}, return results as a \code{GRanges} object;
    if \code{FALSE}, the results are returned as a \code{data.frame}.}
  \item{mask}{an optional \code{GRanges} object with genomic regions to
    be masked, i.e. excluded from the analysis (e.g. unmappable regions).}
  \item{reference}{source of bam files; can be either \dQuote{genome}
    (then the alignments against the genome are used) or the name of an
    auxiliary target sequence (then alignments against this target
    sequence will be used). The auxiliary name must correspond to the
    name contained in the auxiliary file refered by the
    \code{auxiliaryFile} argument of \code{\link[QuasR]{qAlign}}.}
  \item{keepZero}{if \code{FALSE}, only cytosines covered by at least
    one alignment will be returned; \code{keepZero} must be \code{TRUE}
    if multiple samples have the same sample name and
    \code{collapseBySample} is \code{TRUE}.}
  \item{clObj}{a cluster object to be used for parallel processing of
    multiple files (see \sQuote{Details}).}
}
\details{
  \code{qMeth} cam be used on a \code{qProject} object from a bisulfite
  sequencing experiment (sequencing of bisulfite-converted DNA), such as
  the one returned by \code{\link[QuasR]{qAlign}} when its parameter \code{bisulfite} is
  set to a different value than \dQuote{no}.

  \code{qMeth} quantifies DNA methylation by counting total and
  methylated events for individual cytosines, using the alignments that
  have been generated in converted (three-letter) sequence space for
  example by \code{\link[QuasR]{qAlign}}. A methylated event corresponds to a C/C
  match in the alignment, an unmethylated event to a T/C mismatch (or G/G
  matches and A/G mismatches on the opposite strand).

  Both directed (\code{bisulfite}=\dQuote{dir}) and undirected
  (\code{bisulfite}=\dQuote{undir}) experimental protocols are supported
  by \code{\link[QuasR]{qAlign}} and \code{qMeth}.

  By default, results are returned per C nucleotide. If
  \code{reportLevel}=\dQuote{alignment}, results are reported separately
  for individual alignments. In that case, \code{query} has to be a
  \code{GRanges} object with exactly one region, \code{mode} has to be
  either \dQuote{CpG} or \dQuote{allC}, the arguments
  \code{collapseByQueryRegion}, \code{asGRanges}, \code{mask} and
  \code{keepZero} have no effect and allele-specific projects are
  treated in the same way as normal (non-allele specific) projects.
  
  Using the parameter \code{mode}, quantification can be limited to
  cytosines in CpG context, and counts obtained for the two cytosines on
  opposite strands within a single CpG can be combined (summed).

  The quantification of methylation for all cytosines in the query region(s)
  (\code{mode}=\dQuote{allC}) should be done with care, especially for
  large query regions, as the return value may require a large amount of
  memory.

  If \code{mode} is set to \dQuote{var}, \code{qMeth} only counts reads
  from the strand opposite of the cytosine and reports total and
  matching alignments. For a position identical to the reference
  sequence, only matches (and very few sequencing errors) are
  expected, independent on the methylation state of the cytosine. A
  reduced fraction of alignments matching the reference are indicative
  of sequence variations in the sequenced sample.

  If an object that inherits from class \code{cluster} is provided to
  the \code{clObj} argument, for example an object returned by
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel},
  the quantification task is split into multiple chunks and processed in
  parallel using \code{\link[parallel]{clusterApplyLB}} from package
  \pkg{parallel}. Not all tasks will be efficiently parallelized: For
  example, a single query region and a single (group of) bam files will
  not be split into multiple chunks.
}
\value{
  For \code{reportLevel}=\dQuote{C}, a \code{GRanges} object if
  \code{asGRanges}=\code{TRUE}, otherwise a \code{data.frame}.

  Each row contains the coordinates of individual cytosines for
  \code{collapseByQueryRegion}=\code{FALSE} or query regions
  for \code{collapseByQueryRegion}=\code{TRUE}.

  In addition to the coordinates columns (or \code{seqnames},
  \code{ranges} and \code{strand} slots for \code{GRanges} objects),
  each row contains per bam file:

  Two values (total and methylated events, with suffixes _T and _M), or
  if the \code{qProject} object was created including a SNP table,
  six values (total and methylated events for Reference, Unknown and
  Alternative genotypes, with suffixed _TR, _TU, _TA, _MR, _MU and _MA).
  In the latter case, C's or CpG's that overlap with SNPs in the table
  are removed.

  If \code{collapseBySample}=\code{TRUE}, groups of bam files with
  identical sample name are combined (summed) and will be represented by
  a single set of total and methylated count columns.

  If \code{mode}=\dQuote{var}, the _T and _M columns correspond to total
  and matching alignments overlapping the guanine paired to the
  cytosine.

  For \code{reportLevel}=\dQuote{alignment}, a \code{list} with one
  element per bam file or sample (depending on
  \code{collapseBySample}. Each list element is another list with the
  elements:
  \itemize{
    \item \code{aid}: character vector with unique alignment identifiers
    \item \code{Cid}: integer vector with genomic coordinate of C base
    \item \code{strand}: character vector with the strand of the C base
    \item \code{meth}: integer vector with methylation state for
    alignment and C defined by \code{aid} and \code{Cid}. The values are
    1 for methylated or 0 for unmethylated states.
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
  \code{\link[QuasR]{qAlign}},
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
}
\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create alignments
sampleFile <- "extdata/samples_bis_single.txt"
genomeFile <- "extdata/hg19sub.fa"
proj <- qAlign(sampleFile, genomeFile, bisulfite="dir")
proj

# calculate methylation states
meth <- qMeth(proj, mode="CpGcomb")
meth
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
