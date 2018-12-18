requireNamespace("GenomicRanges")
requireNamespace("parallel")

# copy sample data
file.copy(system.file("extdata", package = "QuasR"), ".", recursive = TRUE)

# create cluster object
clObj <- parallel::makeCluster(2L)

# load QuasR on cluster nodes
parallel::clusterEvalQ(cl = clObj, expr = library(QuasR))

# create annotation and GRanges
# ... for hg19sub
gtfGr <- createGtfGr()
txdb  <- createTxDb()
auxGr <- Rsamtools::scanFaIndex(auxGenomeFile)

# ... for synthetic genome with chrV
qGenome <- GenomicRanges::GRanges("chrV", IRanges(start = 1, width = 800), strand = "+")
q01.20  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(20), width = 1))
q01.99  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(99), width = 1))
q21.40  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(20) + 20, width = 1))
q41.60  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(20) + 40, width = 1))
q61.80  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(20) + 60, width = 1))
q81.99  <- GenomicRanges::GRanges("chrV", IRanges(start = seq.int(19) + 80, width = 1))
qTiles  <- createTiles()

# create qProject instances
# ... for hg19sub
genomeFile    <- file.path("extdata", "hg19sub.fa")
snpFile       <- file.path("extdata", "hg19sub_snp.txt")
auxGenomeFile <- file.path("extdata", "NC_001422.1.fa")
auxFile       <- file.path("extdata", "auxiliaries.txt")

sChipSingle   <- file.path("extdata", "samples_chip_single.txt")
sRnaSingle    <- file.path("extdata", "samples_rna_single.txt")
sRnaPaired    <- file.path("extdata", "samples_rna_paired.txt")
sBisSingle    <- file.path("extdata", "samples_bis_single.txt")
#sMirnaSingle  <- file.path("extdata", "samples_mirna.txt")

pChipSingle       <- qAlign(sChipSingle,  genomeFile, clObj = clObj)
pChipSingleAux    <- qAlign(sChipSingle,  genomeFile, auxiliaryFile = auxFile, clObj = clObj)
pChipSingleSnps   <- qAlign(sChipSingle,  genomeFile, snpFile = snpFile, clObj = clObj)
pRnaPaired        <- qAlign(sRnaPaired,   genomeFile, clObj = clObj)
pRnaSingleSpliced <- qAlign(sRnaSingle,   genomeFile, splicedAlignment = TRUE, aligner = "Rbowtie", clObj = clObj)
pRnaPairedSpliced <- qAlign(sRnaPaired,   genomeFile, splicedAlignment = TRUE, aligner = "Rbowtie", clObj = clObj)
pBis              <- qAlign(sBisSingle,   genomeFile, bisulfite = "dir", clObj = clObj)
pBisUndir         <- qAlign(sBisSingle,   genomeFile, bisulfite = "undir", clObj = clObj)
#pBisSnps          <- qAlign(sBisSingle,   genomeFile, bisulfite = "dir", snpFile = snpFile, clObj = clObj)
pPhiX             <- qAlign(file.path("extdata", "phiX_paired_withSecondary_sampleFile.txt"),
                            file.path("extdata", "NC_001422.1.fa"), paired = "fr")

# ... for synthetic genome with chrV
pSingle        <- createProjectSingle(allelic = FALSE)
pSingleAllelic <- createProjectSingle(allelic = TRUE)
pPaired        <- createProjectPaired()
