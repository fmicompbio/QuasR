\name{qCount}
\alias{qCount}
\title{Quantify alignments}
\description{
Quantify alignments from sequencing data.
}
\usage{
qCount(proj, 
       query,
       reportLevel=c(NULL, "gene", "exon","promoter","junction"),
       selectReadPosition=c("start", "end"),
       shift=0L,
       orientation=c("any", "same", "opposite"),
       useRead=c("any", "first", "last"),
       auxiliaryName=NULL,
       mask=NULL,
       collapseBySample=TRUE,
       includeSpliced=TRUE,
       mapqMin=0L,
       mapqMax=255L,
       maxInsertSize=500L,
       clObj=NULL)
}
\arguments{
  \item{proj}{a \code{\linkS4class{qProject}} object representing a
    sequencing experiment as returned by \code{\link[QuasR]{qAlign}}}
  \item{query}{an object of type \code{\link[GenomicRanges]{GRanges}}, 
    \code{\link[GenomicRanges]{GRangesList}} or 
    \code{\link[GenomicFeatures]{TxDb}} with the regions to be
    quantified. The type of \code{query} will determine the mode of
    quantification (see \sQuote{Details}). For
    \code{reportLevel="junction"}, \code{query} is ignored and can also
    be \code{NULL}.}
  \item{reportLevel}{level of quantification (\code{query} of type
    \code{TxDb} or \code{NULL}), one of
    \itemize{
      \item \code{gene} (default): one value per gene
      \item \code{exon}: one value per exon
      \item \code{promoter}: one value per promoter
      \item \code{junction}: one count per detected exon-exon junction
      (\code{query} will be ignored in this case)
    }
  }
  \item{selectReadPosition}{defines the part of the alignment that has
    to be contained within a query region to produce an overlap (see
    Details). Possible values are:
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
  \item{mapqMin}{minimal mapping quality of alignments to be included when
    counting (mapping quality must be greater than or equal to
    \code{mapqMin}). Valid values are between 0 and 255. The default (0)
    will include all alignments.}
  \item{mapqMax}{maximal mapping quality of alignments to be included when
    counting (mapping quality must be less than or equal to \code{mapqMax}).
    Valid values are between 0 and 255. The default (255) will include
    all alignments.}
  \item{maxInsertSize}{Maximal fragment size of the paired-end experiment. 
    This parameter is used if \code{shift="halfInsert"} and will
    ensure that query regions are made wide enough to emcompass all
    alignment pairs whose mid falls into the query region. The default
    value is \code{500} bases.}
  \item{clObj}{a cluster object to be used for parallel processing (see
    \sQuote{Details}).}
}
\details{
  \code{qCount} is used to count alignments in each sample from a
  \code{qProject} object. The features to be quantified, together with
  the mode of quantification, are specified by the \code{query}
  argument, which is one of:
  \itemize{
    \item \code{\link[GenomicRanges]{GRanges}}: Overlapping alignments
    are counted separately for each coordinate region. If multiple
    regions have identical names, their counts will be summed, counting
    each alignment only once even if it overlaps more than one of these
    regions. Alignments may be counted more than once if they overlap
    multiple regions that have different names.
    This mode is for example used to quantify ChIP-seq alignments in
    promoter regions, or gene expression levels in an RNA-seq experiment
    (using a \code{query} with exon regions named by gene).
    \item \code{\link[GenomicRanges]{GRangesList}}: Alignments are
    counted and summed for each list element in \code{query} if they
    overlap with any of the regions contained in the list element. The
    order of the list elements defines a hierarchy for quantification:
    Alignment will only be counted for the first element (the one with
    the lowest index in \code{query}) that they overlap, but not for any
    potential further list elements containing overlapping regions.
    This mode can be used to hierarchically and uniquely count (assign)
    each alignment to a one of several groups of regions (the elements
    in \code{query}), for example to estimate the fractions of different
    classes of RNA in an RNA-seq experiment (rRNA, tRNA, snRNA, snoRNA,
    mRNA, etc.)
    \item \code{\link[GenomicFeatures]{TxDb}}: Used to extract
    regions from annotation and report alignment counts depending on the
    value of \code{reportLevel}. If \code{reportLevel="exon"},
    alignments overlapping each exon in \code{query} are counted.
    If \code{reportLevel="gene"}, alignment counts for all exons of a
    gene will be summed, counting each alignment only once even if it
    overlaps multiple annotated exons of a gene. These are useful to
    calculate exon or gene expression levels in RNA-seq experiments
    based on the annotation in a \code{TranscriptDB} object. If
    \code{reportLevel="promoter"}, the \code{promoters} function from package
    \pkg{GenomicFeatures} is used with default arguments to extract
    promoter regions around transcript start sites, e.g. to quantify
    alignments inf a ChIP-seq experiment.
    \item any of the above or \code{NULL} for
    \code{reportLevel="junction"}: The \code{query} argument is ignored
    if \code{reportLevel} is set to \code{"junction"}, and \code{qCount}
    will count the number of alignments supporting each exon-exon
    junction detected in any of the samples in \code{proj}. The
    arguments \code{selectReadPosition}, \code{shift},
    \code{orientation}, \code{useRead} and \code{mask} will have no
    effect in this quantification mode.
  }

  The additional arguments allow to fine-tune the quantification:

  \code{selectReadPosition} defines the part of the alignment that has
  to be contained within a query region for an overlap. The values
  \code{start} (default) and \code{end} refer to the biological start
  (5'-end) and end (3'-end) of the alignment. For example, the
  \code{start} of an alignment on the plus strand is its leftmost
  (lowest) base, and the \code{end} of an alignment on the minus strand
  is also the leftmost base.

  \code{shift} allows on-the-fly shifting of alignments towards their
  3'-end prior to overlap determination and counting. This can be
  helpful to increase resolution of ChIP-seq experiments by moving
  alignments by half the immuno-precipitated fragment size towards the
  middle of fragments. \code{shift} is either an \dQuote{integer} vector
  with one value per alignment file in \code{proj}, or a single
  \dQuote{integer} value, in which case all alignment files will be
  shifted by the same value. For paired-end experiments, it can be
  alternatively set to "halfInsert", which will estimate the true
  fragment size from the distance between aligned read pairs and shift
  the alignments accordingly.
  
  \code{orientation} controls the interpretation of alignment strand
  when counting, relative to the strand of the query region. \code{any}
  will count all overlapping alignments, irrespective of the alignment
  strand (e.g. used in an unstranded RNA-seq experiment). \code{same}
  will only count the alignments on the same strand as the query region
  (e.g. in a stranded RNA-seq experiment), and \code{opposite} will only
  count the alignments on the opposite strand from the query region
  (e.g. to quantify anti-sense transcription in a stranded RNA-seq
  experiment).
  
  \code{includeSpliced} can be used to include or exclude spliced
  alignments. \code{mapqMin} and \code{mapqMax} allow to select alignments
  based on their mapping qualities. \code{mapqMin} and \code{mapqMax} can
  take integer values between 0 and 255 and equal to
  \eqn{-10 log10 Pr(mapping position is wrong)}, rounded to the nearest
  integer. A value 255 indicates that the mapping quality is not available.

  In paired-end experiments, \code{useRead} allows to quantify either
  all alignments (\code{useRead="any"}), or only the first
  (\code{useRead="first"}) or last (\code{useRead="last"}) read from a
  read pair or read group. Note that for \code{useRead="any"} (the
  default), an alignment pair that is fully contained within a query
  region will contribute two counts to the value of that region.

  \code{auxiliaryName} selects the reference sequence for which
  alignments should be quantified. \code{NULL} (the default) will
  select alignments against the genome. If set to a character string
  that matches one of the auxiliary target names (as specified in
  the \code{auxiliaryFile} argument of \code{\link[QuasR]{qAlign}}), the
  corresponding alignments will be counted.

  \code{mask} can be used to specify a
  \code{\link[GenomicRanges]{GRanges}} object with regions in the
  reference sequence to be excluded from quantification. The regions
  will be considered unstranded (\code{strand="*"}). Alignments that
  overlap with a region in \code{mask} will not be counted. Masking may
  reduce the effective width of query regions reported by \code{qCount},
  even down to zero for regions that are fully contained in \code{mask}.

  %collapseBySample
  %maxInsertSize
  
  If \code{clObj} is set to an object that inherits from class
  \code{cluster}, for example an object returned by
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}, the
  quantification task is split into multiple chunks and processed in
  parallel using \code{\link[parallel]{clusterMap}}. Currently, not all
  tasks will be efficiently parallelized: For example, a single query
  region and a single (group of) bam files will not be split into
  multiple chunks.
}
\value{
  A \code{matrix} with effective query regions width in the first
  column, and alignment counts in subsequent columns, or a
  \code{GRanges} object if \code{reportLevel="junction"}.

  The effective query region width returned as first column in the
  matrix is calculated by the number of unique, non-masked bases in the
  reference sequence that contributed to the count of this query
  name. An effective width of zero indicates that the region was fully
  masked and will have zero counts in all samples.
  
  The alignment counts in the matrix are contained from column two
  onwards. For projects with allele-specific quantification, i.e. if a
  file with single nucleotide polymorphisms was supplied to the
  \code{snpFile} argument of \code{\link[QuasR]{qAlign}}, there will be
  three columns per bam file (number of alignments for Reference,
  Unknown and Alternative genotypes, with suffixed _R, _U and
  _A). Otherwise there is a single columns per bam file.
  
  If \code{collapseBySample}=\code{TRUE}, groups of bam files with identical 
  sample name are combined by summing their alignment counts.

  For \code{reportLevel="junction"}, the return value is a
  \code{GRanges} object. The start and end coordinates correspond to the
  first and last base in each detected intron. Plus- and minus-strand
  alignments are quantified separately, so that in an unstranded RNA-seq
  experiment, the same intron may be represented twice; once for each
  strand. The counts for each sample are contained in the \code{mcols}
  of the \code{GRanges} object.
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
  \code{\linkS4class{qProject}},
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
}
\examples{
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

# copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)


# load genome sequence
genomeFile <- "extdata/hg19sub.fa"
gseq <- readDNAStringSet(genomeFile)
chrRegions <- GRanges(names(gseq), IRanges(start=1,width=width(gseq),names=names(gseq)))


# create alignments (paired-end experiment)
sampleFile <- "extdata/samples_rna_paired.txt"
proj <- qAlign(sampleFile, genomeFile, splicedAlignment=TRUE)


# count reads using a "GRanges" query
qCount(proj, query=chrRegions)
qCount(proj, query=chrRegions, useRead="first")


# hierarchical counting using a "GRangesList" query
library(rtracklayer)
annotationFile <- "extdata/hg19sub_annotation.gtf"
gtfRegions <- import.gff(annotationFile, format="gtf", asRangedData=FALSE,
                    feature.type="exon")
names(gtfRegions) <- mcols(gtfRegions)$source
gtfRegionList <- split(gtfRegions, names(gtfRegions))
names(gtfRegionList)

res3 <- qCount(proj, gtfRegionList)
res3


# gene expression levels using a "TxDb" query
library("GenomicFeatures")
genomeRegion <- scanFaIndex(genomeFile)
chrominfo <- data.frame(chrom=as.character(seqnames(genomeRegion)),
                        length=end(genomeRegion),
                        is_circular=rep(FALSE, length(genomeRegion)))
txdb <- makeTranscriptDbFromGFF(annotationFile, 
                                format="gtf", 
                                exonRankAttributeName="exon_number",  
                                gffGeneIdAttributeName="gene_name", 
                                chrominfo=chrominfo,
                                dataSource="Ensembl modified",
                                species="Homo sapiens")

res4 <- qCount(proj, txdb, reportLevel="gene")
res4

# exon-exon junctions
res5 <- qCount(proj, NULL, reportLevel="junction")
res5
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
