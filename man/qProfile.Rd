\name{qProfile}
\alias{qProfile}
\title{Quantify alignments by relative position}
\description{
Quantify alignments from sequencing data, relative to their position in query regions.
}
\usage{
qProfile(proj, 
         query,
         upstream=1000,
         downstream=upstream,
         selectReadPosition=c("start", "end"),
         shift=0L,
         orientation=c("any", "same", "opposite"),
         useRead=c("any", "first", "last"),
         auxiliaryName=NULL,
         mask=NULL,
         collapseBySample=TRUE,
         includeSpliced=TRUE,
         includeSecondary=TRUE,
         mapqMin=0L,
         mapqMax=255L,
         absIsizeMin=NULL,
         absIsizeMax=NULL,
         maxInsertSize=500L,
         clObj=NULL)
}
\arguments{
  \item{proj}{a \code{\linkS4class{qProject}} object representing a
    sequencing experiment as returned by \code{\link[QuasR]{qAlign}}}
  \item{query}{an object of type \code{\link[GenomicRanges]{GRanges}}
    with the regions to be profiled. All regions in \code{query} will be
    anchored at their biological start position (\code{start(query)} for
    regions on strand \dQuote{+} or \dQuote{*}, \code{end(query)} for
    regions on strand \dQuote{-}). This position will become position
    zero in the return value.}
  \item{upstream}{An \dQuote{integer} vector of length one or the same
    length as \code{query} indicating the number of bases upstream of the
    anchor position to include in the profile.}
  \item{downstream}{An \dQuote{integer} vector of length one or the same
    length as \code{query} indicating the number of bases downstream of the
    anchor position to include in the profile.}
  \item{selectReadPosition}{defines the part of the alignment that has
    to be contained within a query region to produce an overlap (see
    Details), and that is used to calculate the relative position within
    the query region. Possible values are:
    \itemize{
      \item \code{start} (default): start of the alignment
      \item \code{end}: end of the alignment
    }
  }
  \item{shift}{controls the shifting alignments towards their 3'-end before
    quantification. \code{shift} can be one of:
    \itemize{
      \item an \dQuote{integer} vector of the same length as the
      number of alignment files
      \item a single \dQuote{integer} value
      \item the character string \code{"halfInsert"} (only available for
      paired-end experiments)
    }
    The default of \code{0} will not shift any alignments.}
  \item{orientation}{sets the required orientation of the alignments relative to the query
    region in order to be counted, one of:
    \itemize{
      \item \code{any} (default): count alignment on the same and opposite strand
      \item \code{same} : count only alignment on the same strand
      \item \code{opposite} : count only alignment on the opposite strand
    }
  }
  \item{useRead}{For paired-end experiments, selects the read mate whose
    alignments should be counted, one of:
    \itemize{
      \item \code{any} (default): count all alignments
      \item \code{first} : count only alignments from the first read
      \item \code{last} : count only alignments from the last read
    }
  }
  \item{auxiliaryName}{which bam files to use in an experiments with
    auxiliary alignments (see Details).}
  \item{mask}{If not \code{NULL}, a \code{\link[GenomicRanges]{GRanges}}
    object with reference regions to be masked, i.e. excluded from the
    quantification, such as unmappable or highly repetitive regions (see
    Details).}
  \item{collapseBySample}{if \code{TRUE} (the default), sum alignment
    counts from bam files with the same sample name.}
  \item{includeSpliced}{if \code{TRUE} (the default), include spliced
    alignments when counting. A spliced alignment is defined as an
    alignment with a gap in the read of at least 60 bases.}
  \item{includeSecondary}{if \code{TRUE} (the default), include alignments
    with the secondary bit (0x0100) set in the \code{FLAG} when counting.}
  \item{mapqMin}{minimal mapping quality of alignments to be included when
    counting (mapping quality must be greater than or equal to
    \code{mapqMin}). Valid values are between 0 and 255. The default (0)
    will include all alignments.}
  \item{mapqMax}{maximal mapping quality of alignments to be included when
    counting (mapping quality must be less than or equal to \code{mapqMax}).
    Valid values are between 0 and 255. The default (255) will include
    all alignments.}
  \item{absIsizeMin}{For paired-end experiments, minimal absolute insert
    size (TLEN field in SAM Spec v1.4) of alignments to be included when
    counting. Valid values are greater than 0 or \code{NULL} (default),
    which will not apply any minimum insert size filtering.}
  \item{absIsizeMax}{For paired-end experiments, maximal absolute insert
    size (TLEN field in SAM Spec v1.4) of alignments to be included when
    counting. Valid values are greater than 0 or \code{NULL} (default),
    which will not apply any maximum insert size filtering.}
  \item{maxInsertSize}{Maximal fragment size of the paired-end experiment. 
    This parameter is used if \code{shift="halfInsert"} and will
    ensure that query regions are made wide enough to emcompass all
    alignment pairs whose mid falls into the query region. The default
    value is \code{500} bases.}
  \item{clObj}{a cluster object to be used for parallel processing (see
    \sQuote{Details}).}
}
\details{
  \code{qProfile} is used to count alignments in each sample from a
  \code{qProject} object, relative to their position in query regions.
  
  Most arguments are identical to the ones of \code{\link[QuasR]{qCount}}.

  The \code{query} argument is a \code{\link[GenomicRanges]{GRanges}}
  object that defines the regions for the profile. All regions in
  \code{query} will be aligned to one another at their anchor position,
  which corresponds to their biological start position (\code{start(query)}
  for regions on strand \dQuote{+} or \dQuote{*}, \code{end(query)} for
  regions on strand \dQuote{-}).
  
  This anchor position will be extended (with regard to strand) by
  the number of bases specified by \code{upstream} and \code{downstream}.
  In the return value, the anchor position will be at position zero.
  
  Regions with identical names in \code{names{query}} will be summed, and
  profiles will be padded with zeros to accomodate the length of all 
  profiles (\code{max(upstream)+max(downstream)+1}).
}
\value{
  A \code{list} of matrices with \code{length(unique(names(query)))} rows
  with profile names, and \code{max(upstream)+max(downstream)+1} columns
  indicating relative position.
  
  The first list element is called \dQuote{coverage} and contains, for each
  profile and relative position, the number of overlapping regions
  that contributed to the profile.
  
  Subsequent list elements contain the alignment counts for individual
  sequence files (\code{collapseBySample=FALSE}) or samples
  (\code{collapseBySample=TRUE}) in \code{proj}.

  For projects with allele-specific quantification, i.e. if a file with
  single nucleotide polymorphisms was supplied to the \code{snpFile}
  argument of \code{\link[QuasR]{qAlign}}, there will be three rows
  instead of one row with counts per unique region name, with numbers
  of alignments for Reference, Unknown and Alternative genotypes
  (suffixed _R, _U and _A).
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
  \code{\link[QuasR]{qCount}},
  \code{\link[QuasR]{qAlign}},
  \code{\linkS4class{qProject}},
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
}
\examples{
# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create alignments (single-end experiment)
genomeFile <- "extdata/hg19sub.fa"
sampleFile <- "extdata/samples_chip_single.txt"
proj <- qAlign(sampleFile, genomeFile)

# load transcript start site coordinates
library(rtracklayer)
annotationFile <- "extdata/hg19sub_annotation.gtf"
tssRegions <- import.gff(annotationFile, format="gtf",
                         asRangedData=FALSE,
                         feature.type="start_codon")

# obtain a combined TSS profile
pr1 <- qProfile(proj, tssRegions)
lapply(pr1, dim)
lapply(pr1, "[", , 1:5)

prComb <- do.call("+", lapply(pr1[-1], function(x) x/pr1[[1]]))
barplot(prComb, xlab="Position", ylab="Mean no. of alignments")

# obtain TSS profiles for individual regions
names(tssRegions) <- mcols(tssRegions)$transcript_id
pr2 <- qProfile(proj, tssRegions)
lapply(pr2, dim)
lapply(pr2, "[", 1:3, 1:5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
