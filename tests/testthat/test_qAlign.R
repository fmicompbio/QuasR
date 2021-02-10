
context("alignmentStats")

test_that("alignmentStats correctly digests its arguments", {
  expect_error(alignmentStats(1L))
  expect_error(alignmentStats("nonexistent.bam"))
})

test_that("alignmentStats works as expected", {
  resSoll <- matrix(c(95000,95000,
                      2339,3609,
                      258,505), nrow = 2, ncol = 3)
  res <- alignmentStats(pChipSingle)
  expect_identical(resSoll, unname(res))
  res <- alignmentStats(alignments(pChipSingle)$genome$FileName)
  expect_identical(resSoll, unname(res))

  resSoll <- matrix(c(95000,95000,5386,5386,
                      2339,3609,251,493,
                      258,505,7,12), nrow = 4, ncol = 3)
  res <- alignmentStats(pChipSingleAux)
  expect_identical(resSoll, unname(res))
})


context("qAlign")

test_that("qAlign digests arguments correctly", {
  # test redundant sample files
  tmpf1 <- tempfile(fileext = ".txt")
  tmpf2 <- tempfile(fileext = ".txt")
  tab1 <- pChipSingle@reads
  tab2 <- pChipSingle@alignments
  write.table(tab1[c(1,1), ], tmpf1, sep = "\t", row.names = FALSE)
  write.table(tab2[c(1,1), ], tmpf2, sep = "\t", row.names = FALSE)
  expect_error(qAlign(tmpf1, genomeFile))
  expect_error(qAlign(tmpf2, genomeFile))

  # other arguments
  expect_error(qAlign("nonexistent_file"))
  expect_error(qAlign(genome = genomeFile))
  expect_error(qAlign(sChipSingle, genome = genomeFile[c(1,1)]))
  expect_error(qAlign(sampleFile = sChipSingle))
  expect_error(qAlign(sChipSingle, genomeFile, splicedAlignment = TRUE, checkOnly = TRUE))
  expect_error(qAlign(sChipSingle, genomeFile, aligner = "unsupportedAligner"))
  expect_error(qAlign(sChipSingle, genomeFile, lib.loc = "notADirectory"))
  expect_error(qAlign(sChipSingle, genomeFile, alignmentsDir = "notADirectory"))
  expect_error(qAlign(sChipSingle, genomeFile, aligner = "Rhisat2", bisulfite = "dir"))
  expect_error(qAlign(sChipSingle, genomeFile, aligner = "Rbowtie", splicedAlignment = TRUE, bisulfite = "dir"))
  expect_error(qAlign(sRnaPaired, genomeFile, aligner = "Rbowtie", splicedAlignment = TRUE, paired = "rf"))
  expect_error(qAlign(sBisSingle, genomeFile, bisulfite = "invalid"))
  expect_error(qAlign(sRnaPaired, genomeFile, bisulfite = "dir", paired = "rf"))
  expect_error(qAlign(sRnaPaired, genomeFile, auxiliaryFile = "nonExisting"))
  expect_error(qAlign(sRnaPaired, genomeFile, aligner = "Rhisat2", geneAnnotation = "nonExisting"))
  expect_error(qAlign(sRnaPaired, genomeFile, paired = "invalid"))
  expect_error(qAlign(sampleFile = "nonExisting", genomeFile))
  expect_error(qAlign(sChipSingle, genomeFile, paired = "fr"))
  
  # inconsistent bam header and BSgenome reference
  # ... create bam file with inconsistent header
  pBad <- pChipSingle
  goodbamfile <- pChipSingle@alignments$FileName[1]
  goodsamfile <- Rsamtools::asSam(goodbamfile, tempfile())
  goodsamlines <- readLines(goodsamfile)
  badsamlines <- c(goodsamlines[1:4], "@SQ\tSN:chr4\tLN:9999", goodsamlines[5:length(goodsamlines)])
  badsamfile <- tempfile(fileext = ".sam")
  writeLines(badsamlines, badsamfile)
  badbamfile <- Rsamtools::asBam(badsamfile, file.path(dirname(goodbamfile), paste0("bad_", sub(".bam$", "", basename(goodbamfile)))))
  # ... incorporate it into a sample file
  goodfnames <- c(pChipSingle@reads$FileName[1], paste0(pChipSingle@alignments$FileName[1], ".txt"))
  badfnames <- file.path(dirname(goodfnames), paste0("bad_", basename(goodfnames)))
  file.copy(goodfnames, badfnames)
  txtlines <- readLines(badfnames[2])
  txtlines[1] <- sub("chip_1_1", "bad_chip_1_1", txtlines[1])
  writeLines(txtlines, badfnames[2])
  badsamplefile <- tempfile(tmpdir = dirname(goodfnames[1]), fileext = ".txt")
  samplelines <- readLines(sChipSingle)
  samplelines[2] <- paste0("bad_", samplelines[2])
  writeLines(samplelines, badsamplefile)
  expect_warning(qAlign(badsamplefile, genomePkg, clObj = clObj, lib.loc = rlibdir, cacheDir = tempdir()))
})

test_that("qAlign correctly works for single reads", {
  aln <- GenomicAlignments::readGAlignments(pChipSingle@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 1073L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(590, 708, 1041)))))
  expect_identical(sum(as.numeric(start(aln))), 36493854)
  expect_identical(sum(as.numeric(end(aln))), 36575719)
})

test_that("qAlign correctly works for single reads with auxiliaries", {
  aln <- GenomicAlignments::readGAlignments(pChipSingleAux@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 1073L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(590, 708, 1041)))))
  expect_identical(sum(as.numeric(start(aln))), 36493854)
  expect_identical(sum(as.numeric(end(aln))), 36575719)

  aln <- GenomicAlignments::readGAlignments(pChipSingleAux@auxAlignments["phiX174","Sample2"], use.names = TRUE)
  expect_length(runValue(strand(aln)), 238L)
  expect_identical(seqnames(aln), Rle(factor(rep("phiX174", 493))))
  expect_identical(sum(as.numeric(start(aln))), 1339781)
  expect_identical(sum(as.numeric(end(aln))), 1357036)
})

test_that("qAlign correctly works for paired-end reads", {
  aln <- GenomicAlignments::readGAlignments(pRnaPaired@alignments$FileName[2], use.names = TRUE)
  expect_length(runValue(strand(aln)), 394L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(898, 64, 1690)))))
  expect_identical(sum(as.numeric(start(aln))), 54091733)
  expect_identical(sum(as.numeric(end(aln))), 54221681)
})

test_that("qAlign correctly works for single reads (spliced, Rbowtie)", {
  aln <- GenomicAlignments::readGAlignments(pRnaSingleSpliced@alignments$FileName[3], use.names = TRUE)
  expect_length(runValue(strand(aln)), 435L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(732, 1112, 1156)))))
  expect_identical(sum(as.numeric(start(aln))), 42121859)
  expect_identical(sum(as.numeric(end(aln))), 43006110)
})

test_that("qAlign correctly works for paired-end reads (spliced, Rbowtie)", {
  aln <- GenomicAlignments::readGAlignments(pRnaPairedSpliced@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 1620L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(850, 2924, 2228)))))
  expect_identical(sum(as.numeric(start(aln))), 69243139)
  expect_identical(sum(as.numeric(end(aln))), 71021433)
})

test_that("qAlign correctly works for single reads (spliced, without splice site file, Rhisat2)", {
  aln <- GenomicAlignments::readGAlignments(pRnaSingleSplicedHisat2@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 601L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(420, 1457, 1104)))))
  expect_identical(sum(as.numeric(start(aln))), 34122868)
  expect_identical(sum(as.numeric(end(aln))), 35028897)
})

test_that("qAlign correctly works for paired reads (unspliced, without splice site file, Rhisat2)", {
  aln <- GenomicAlignments::readGAlignments(pRnaPairedUnsplicedHisat2@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 868L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(766, 2055, 2089)))))
  expect_identical(sum(as.numeric(start(aln))), 62093887)
  expect_identical(sum(as.numeric(end(aln))), 62334263)
})

test_that("qAlign correctly works for paired reads (spliced, with splice site file from gtf, Rhisat2)", {
  aln <- GenomicAlignments::readGAlignments(pRnaPairedSplicedHisat2Gtf@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 1000L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(839, 2913, 2207)))))
  expect_identical(sum(as.numeric(start(aln))), 68278227)
  expect_identical(sum(as.numeric(end(aln))), 70288159)
})

test_that("qAlign correctly works for paired reads (spliced, with splice site file from TxDb, Rhisat2)", {
  aln <- GenomicAlignments::readGAlignments(pRnaPairedSplicedHisat2TxDb@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 1000L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(839, 2913, 2207)))))
  expect_identical(sum(as.numeric(start(aln))), 68278227)
  expect_identical(sum(as.numeric(end(aln))), 70288159)
})

test_that("qAlign correctly works in allelic mode", {
  aln <- GenomicAlignments::readGAlignments(pChipSingleSnps@alignments$FileName[1], use.names = TRUE,
                                            param = Rsamtools::ScanBamParam(tag = "XV"))
  expect_length(runValue(strand(aln)), 931L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(503, 634, 905)))))
  expect_identical(sum(as.numeric(start(aln))), 31658312)
  expect_identical(sum(as.numeric(end(aln))), 31729782)
  expect_equal(as.vector(table(mcols(aln)$XV)[c("A","R","U")]), c(13, 192, 1837))
})

test_that("qAlign correctly works in directional bisulfite mode", {
  aln <- GenomicAlignments::readGAlignments(pBis@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 8154L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(5973, 3997, 15528)))))
  expect_identical(sum(as.numeric(start(aln))), 494072411)
  expect_identical(sum(as.numeric(end(aln))), 495957570)
})

test_that("qAlign correctly works in undirectional bisulfite mode", {
  aln <- GenomicAlignments::readGAlignments(pBisUndir@alignments$FileName[1], use.names = TRUE)
  expect_length(runValue(strand(aln)), 8160L)
  expect_identical(seqnames(aln), Rle(factor(rep(paste0("chr",1:3), c(5973, 3997, 15526)))))
  expect_identical(sum(as.numeric(start(aln))), 493991845)
  expect_identical(sum(as.numeric(end(aln))), 495876924)
})


context("qAlign helper functions")

test_that("pathAsAbsoluteRedirected works as expected", {
    expect_null(QuasR:::pathAsAbsoluteRedirected("test.txt", "../"))
})

test_that("resolveCacheDir works as expected", {
    expect_identical("extdata", QuasR:::resolveCacheDir("extdata"))
    expect_identical(tempdir(), QuasR:::resolveCacheDir(NA))
})

test_that("determineSamplesFormat correctly digests its arguments", {
    expect_error(QuasR:::determineSamplesFormat(c("file1.notSupported")))
    expect_error(QuasR:::determineSamplesFormat(c("file1.fasta", "file2.fastq")))
})
