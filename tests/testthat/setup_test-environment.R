# copy sample data
file.copy(system.file("extdata", package = "QuasR"), ".", recursive = TRUE)

# temporary R library
rlibdir     <- tempfile(pattern = "Rlib", tmpdir = "extdata")
dir.create(rlibdir)
# ... add rlibdir to R_LIBS for "R CMD INSTALL" and cluster nodes to find it 
# ... add rlibdir to R_LIBS for "R CMD INSTALL" and cluster nodes to find it 
oldRlibs <- Sys.getenv("R_LIBS")
Sys.setenv(R_LIBS = paste(tools::file_path_as_absolute(rlibdir), oldRlibs,
                          sep = .Platform$path.sep))

# create cluster object
clObj <- parallel::makeCluster(2L)

# load QuasR on cluster nodes
parallel::clusterEvalQ(cl = clObj, expr = library(QuasR))

# install BSgenome.HSapiens.QuasR.hg19sub into temporary library
bsgPkg <- system.file("extdata", "BSgenome.HSapiens.QuasR.hg19sub_0.1.0.tar.gz",
                      package = "QuasR")
utils::install.packages(pkgs = bsgPkg, lib = rlibdir, repos = NULL,
                        type = "source", INSTALL_opts = "--no-test-load")

# copy one of the ChIP fastq files and give it a long name, to test the 
# truncation rules
# also generate a corresponding sample file
file.copy(from = file.path("extdata", "chip_1_1.fq.bz2"), 
          to = file.path("extdata", 
                         "chip_1_1_with_long_file_name_that_will_be_truncated.fq.bz2"))
tmp <- read.delim(file.path("extdata", "samples_chip_single.txt"))
tmp$FileName <- sub("chip_1_1", 
                    "chip_1_1_with_long_file_name_that_will_be_truncated", 
                    tmp$FileName)
write.table(tmp, file = file.path("extdata", "samples_chip_single_longfname.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
sChipSingleLongFname <- file.path("extdata", "samples_chip_single_longfname.txt")

# create qProject instances
# ... for hg19sub
genomePkg     <- "BSgenome.HSapiens.QuasR.hg19sub"
genomeFile    <- file.path("extdata", "hg19sub.fa")
snpFile       <- file.path("extdata", "hg19sub_snp.txt")
auxGenomeFile <- file.path("extdata", "NC_001422.1.fa")
auxFile       <- file.path("extdata", "auxiliaries.txt")
gtfFile       <- file.path("extdata", "hg19sub_annotation.gtf")
txdbFile      <- file.path("extdata", "hg19sub_annotation.sqlite")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtfFile, format = "gtf")
AnnotationDbi::saveDb(txdb, file = txdbFile)

sChipSingle   <- file.path("extdata", "samples_chip_single.txt")
sRnaSingle    <- file.path("extdata", "samples_rna_single.txt")
sRnaPaired    <- file.path("extdata", "samples_rna_paired.txt")
sBisSingle    <- file.path("extdata", "samples_bis_single.txt")
#sMirnaSingle  <- file.path("extdata", "samples_mirna.txt")

# ... ... as a BSgenome
pChipSingle       <- qAlign(sChipSingle,  genomePkg,  clObj = clObj, lib.loc = rlibdir, cacheDir = tempdir())
pChipSingleLongFname <- qAlign(sChipSingleLongFname,  genomePkg,  clObj = clObj, lib.loc = rlibdir, cacheDir = tempdir())
pBis              <- qAlign(sBisSingle,   genomePkg,  bisulfite = "dir", clObj = clObj, lib.loc = rlibdir)

# ... ... as a fasta genome
pChipSingleAux    <- qAlign(sChipSingle,  genomeFile, auxiliaryFile = auxFile, clObj = clObj)
pChipSingleSnps   <- qAlign(sChipSingle,  genomeFile, snpFile = snpFile, clObj = clObj)
pRnaPaired        <- qAlign(sRnaPaired,   genomeFile, clObj = clObj)
pRnaSingleSpliced <- qAlign(sRnaSingle,   genomeFile, splicedAlignment = TRUE, aligner = "Rbowtie", clObj = clObj)
pRnaPairedSpliced <- qAlign(sRnaPaired,   genomeFile, splicedAlignment = TRUE, aligner = "Rbowtie", clObj = clObj)

pRnaSingleSplicedHisat2 <- qAlign(sRnaSingle, genomeFile, alignmentsDir = "extdata", 
                                  splicedAlignment = TRUE, aligner = "Rhisat2", clObj = clObj)
pRnaPairedUnsplicedHisat2 <- qAlign(sRnaPaired, genomeFile, splicedAlignment = FALSE, 
                                    aligner = "Rhisat2", clObj = clObj)
pRnaPairedSplicedHisat2Gtf <- qAlign(sRnaPaired, genomeFile, splicedAlignment = TRUE,
                                     aligner = "Rhisat2", clObj = clObj, geneAnnotation = gtfFile)
pRnaPairedSplicedHisat2TxDb <- qAlign(sRnaPaired, genomeFile, splicedAlignment = TRUE, 
                                      aligner = "Rhisat2", clObj = clObj, geneAnnotation = txdbFile)

pBis              <- qAlign(sBisSingle,   genomeFile, bisulfite = "dir", clObj = clObj)

pBisUndir         <- qAlign(sBisSingle,   genomeFile, bisulfite = "undir", clObj = clObj)
pBisSnps          <- qAlign(sBisSingle,   genomeFile, bisulfite = "dir", snpFile = snpFile, clObj = clObj)
pPhiX             <- qAlign(file.path("extdata", "phiX_paired_withSecondary_sampleFile.txt"),
                            file.path("extdata", "NC_001422.1.fa"), paired = "fr")

# ... for synthetic genome with chrV
pSingle        <- createProjectSingle(allelic = FALSE)
pSingleAllelic <- createProjectSingle(allelic = TRUE)
pPaired        <- createProjectPaired()

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

