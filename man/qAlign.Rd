\name{qAlign}
\alias{qAlign}
\title{Align reads}
\description{
Create read alignments against reference genome and optional auxiliary
targets if not yet existing. If necessary, also build target indices for
the aligner.}
\usage{
qAlign(sampleFile, 
       genome, 
       auxiliaryFile=NULL, 
       aligner="Rbowtie", 
       maxHits=1, 
       paired=NULL, 
       splicedAlignment=FALSE, 
       snpFile=NULL, 
       bisulfite="no", 
       alignmentParameter=NULL, 
       projectName="qProject", 
       alignmentsDir=NULL, 
       lib.loc=NULL, 
       cacheDir=NULL, 
       clObj=NULL)
}
\arguments{
  \item{sampleFile}{the name of a text file listing input sequence files
    and sample names (see \sQuote{Details}).}
  \item{genome}{the reference genome for primary alignments, one of:
    \itemize{
      \item a string referring to a \dQuote{BSgenome} package
      (e.g. \dQuote{"BSgenome.Hsapiens.UCSC.hg19"}), which will be
      downloaded automatically from Bioconductor if not present
      \item the name of a fasta sequence file containing one or
      several sequences (chromosomes) to be used as a reference. The aligner
      index will be created when neccessary and stored in a default
      location (see \sQuote{Details}).
    }
  }
  \item{auxiliaryFile}{the name of a text file listing sequences to be
    used as additional targets for alignment of reads not mapping to the
    reference genome (see \sQuote{Details}).}
  \item{aligner}{selects the aligner program to be used for aligning the
    reads. Currently, only \dQuote{Rbowtie} is supported, which is an R
    wrapper package for \sQuote{bowtie} and \sQuote{SpliceMap} (see
    \code{\link[Rbowtie]{Rbowtie}} package).}
  \item{maxHits}{sets the maximal number of allowed mapping positions
    per read (default: 1). If a read produces more than \code{maxHits}
    alignments, no alignments will be reported for it. In case of a
    multi-mapping read, a single alignment is randomly selected}
  \item{paired}{defines the type of paired-end library and can be set to
    one of \code{no} (single read experiment, default), \code{fr} (fw/rev),
    \code{ff} (fw/fw) or \code{rf} (rev/fw).}
  \item{splicedAlignment}{if \code{TRUE}, reads will be aligned by
    SpliceMap to produce spliced alignments (without using a database of
    known exon-exon junctions). Using splicedAlignment=TRUE will
    increase alignment times roughly by a factor of ten. The option can
    only be used for reads with a minimal length of 50nt; SpliceMap
    ignores reads that are shorter. Such short reads will not be
    contained in the BAM file, neither as mapped or unmapped reads.}
  \item{snpFile}{the name of a text file listing single nucleotide
    polymorphisms to be used for allele-specific alignment and
    quantification (see \sQuote{Details}).}
  \item{bisulfite}{for bisulfite-converted samples (Bis-seq), the type
    of bisulfite library (\dQuote{dir} for directional libraries,
    \dQuote{undir} for undirectional libraries).}
  \item{alignmentParameter}{a optional string containing command line
    parameters to be used for the aligner, to overrule the default
    alignment parameters used by \code{QuasR}. Please use with caution;
    some alignment parameters may brake assumptions made by
    \code{QuasR}. Default parameters are listed in \sQuote{Details}.}
  \item{projectName}{an optional name for the alignment project.}
  \item{alignmentsDir}{the directory to be used for storing alignments
    (bam files). If set to \code{NULL} (default), bam files will be
    generated at the location of the input sequence files.}
  \item{lib.loc}{can be used to change the default library path of
    R. The library path is used by \code{QuasR} to store aligner index
    packages created from \code{BSgenome} reference genomes, or to
    install newly downloaded \code{BSgenome} packages.}
  \item{cacheDir}{specifies the location to store (potentially huge)
    temporary files. If set to \code{NULL} (default), the temporary
    directory of the current R session as returned by \code{tempdir()}
    will be used.}
  \item{clObj}{a cluster object, created by the package
    \code{\link[parallel]{parallel}}, to enable parallel processing and
    speed up the alignment process.}
}
\details{
  Before generating new alignments, \code{qAlign} looks for previously
  generated alignments as well as for an aligner index. If no aligner
  index exists, it will be automatically created and stored in the same
  directory as the provided fasta file, or as an R package in the case
  of a BSgenome reference. The name of this R package will be the same
  as the BSgenome package name, with an additional suffix from the
  aligner (e.g. \code{BSgenome.Hsapiens.UCSC.hg19.Rbowtie}).
  
  \code{sampleFile} is a tab-delimited text file listing all the input
  sequences to be included in a given analysis. The file has either two
  (single-end) or three columns (paired-end). The first row contains the
  column names, and additional rows contain relative or absolute path
  and name of input sequence file(s), as well as the according sample
  name. Three input file formats are supported (fastq, fasta and
  bam). All input files in one \code{sampleFile} need to be in the same
  format, and are recognized by their extension (.fq, .fastq, .fa,
  .fasta, .fna, .bam), in raw or compressed form (e.g. .fastq.gz). If
  bam files are provided, then no alignments are generated by
  \code{qAlign}, and the alignments contained in the bam files will be
  used instead. 

  The column names in \code{sampleFile} have to match to the ones in the
  examples below, for a single-read experiment:
  \tabular{ll}{
    FileName        \tab SampleName \cr
    chip_1_1.fq.bz2 \tab Sample1    \cr
    chip_2_1.fq.bz2 \tab Sample2
  }
  and for a paired-end experiment:     
  \tabular{lll}{
    FileName1      \tab FileName2      \tab SampleName \cr
    rna_1_1.fq.bz2 \tab rna_1_2.fq.bz2 \tab Sample1    \cr
    rna_2_1.fq.bz2 \tab rna_2_2.fq.bz2 \tab Sample2
  }

  The \dQuote{SampleName} column is the human-readable name for each
  sample that will be used as sample labels. Multiple sequence files may
  be associated to the same sample name, which instructs \code{QuasR} to
  combine those files.

  \code{auxiliaryFile} is a tab-delimited text file listing one or
  several additional target sequence files in fasta format. Reads that
  do not map against the reference genome will be aligned against each
  of these target sequence files. The first row contains the column
  names which have to match to the ones in the example below:
  \tabular{ll}{
    FileName       \tab AuxName \cr
    NC_001422.1.fa \tab phiX174
  }

  \code{snpFile} is a tab-delimited text file without a header and
  contains four columns with chromosome name, position, reference allele
  and alternative allele, as in the example below:
  \tabular{llll}{
    chr1 \tab  8596 \tab G \tab A \cr
    chr1 \tab 18443 \tab G \tab A \cr
    chr1 \tab 18981 \tab C \tab T \cr
    chr1 \tab 19341 \tab G \tab A
  }

  The reference and alternaive alleles will be injected into the
  reference genome, resulting in two separate genomes. All reads will be
  aligned separately to both of these genomes, and the alignments will
  be combined, only retaining the best alignment for each read. In the
  final alignment, each read will be marked with a tag that classifies
  it into reference (\code{R}), alternative (\code{A}) or unknown
  (\code{U}), if the reads maps equally well to both genomes.

  If \code{bisulfite} is set to \dQuote{dir} or \dQuote{undir}, reads
  will be C-to-T converted and aligned to a similarly converted genome.

  If \code{alignmentParameter} is \code{NULL} (recommended),
  \code{qAlign} will select default parameters that are suitable for the
  experiment type. Please note that for bisulfite or allele-specific
  experiments, each read is aligned multiple times, and resulting
  alignments need to be combined. This requires special settings for the
  alignment parameters that are not recommended to be changed. For
  \sQuote{simple} experiments (neither bisulfite, allele-specific, nor
  spliced), alignments are generated using the parameters \dQuote{-m
  maxHits --best --strata}. This will align reads with up to
  \dQuote{maxHits} best hits in the genome and selects one of them randomly.
}
\value{
  A \code{\link[QuasR]{qProject}} object.
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
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel},
  \code{\link[Rbowtie]{Rbowtie}} package
}
\examples{
\dontrun{
    # see qCount, qMeth and qProfile manual pages for examples
    example(qCount)
    example(qMeth)
    example(qProfile)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{misc}% __ONLY ONE__ keyword per line
