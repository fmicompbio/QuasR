# cover parts that are not already covered indirectly by other tests

context("compiled code")

test_that("extractUnmappedReads works as expected", {
  fun    <- function(...) .Call(QuasR:::extractUnmappedReads, ...)
  bamf1  <- pChipSingle@alignments$FileName[1]
  bamf2  <- pRnaPaired@alignments$FileName[1]
  fafile <- tempfile(fileext = rep(".fa",2), tmpdir = "extdata")
  fqfile <- tempfile(fileext = rep(".fq",2), tmpdir = "extdata")
  
  # arguments
  expect_error(fun(1L,         fafile[1],             FALSE, FALSE))
  expect_error(fun(bamf1,             1L,             FALSE, FALSE))
  expect_error(fun(bamf1,      fafile[1],                1L, FALSE))
  expect_error(fun(bamf1,      fafile[1],             FALSE,    1L))
  expect_error(fun("no-file",  fafile[1],             FALSE, FALSE))
  expect_error(fun(genomeFile, fafile[1],             FALSE, FALSE))
  expect_error(fun(bamf1,      "not-there/fafile", FALSE, FALSE))
  
  # results
  f <- fun(bamf1, fafile[1], FALSE, FALSE); expect_length(readLines(f), 516L)
  f <- fun(bamf1, fafile[1], FALSE, TRUE);  expect_length(readLines(f), 516L)
  f <- fun(bamf1, fqfile[1], TRUE,  FALSE); expect_length(readLines(f), 1032L)
  f <- fun(bamf1, fqfile[1], TRUE,  TRUE);  expect_length(readLines(f), 1032L)
  f <- fun(bamf2, fafile,    FALSE, FALSE); expect_length(readLines(f[1]), 3746L)
  f <- fun(bamf2, fafile,    FALSE, TRUE);  expect_length(readLines(f[2]), 3746L)
  f <- fun(bamf2, fqfile,    TRUE,  FALSE); expect_length(readLines(f[1]), 7492L)
  f <- fun(bamf2, fqfile,    TRUE,  TRUE);  expect_length(readLines(f[2]), 7492L)
})

test_that("catBam works as expected", {
  requireNamespace("Rsamtools")
  
  fun    <- function(...) .Call(QuasR:::catBam, ...)
  bamf1  <- pSingle@alignments$FileName[1]
  bfile  <- tempfile(fileext = ".bam", tmpdir = "extdata")
  
  # arguments
  expect_error(fun(   1L, bfile))
  expect_error(fun(bamf1,    1L))
  expect_error(fun(bamf1, "not-there/bfile"))

  # results
  expect_identical(fun(bamf1, bfile), 0L)
  Rsamtools::indexBam(bfile)
  expect_identical(alignmentStats(bamf1)[1,],
                   alignmentStats(bfile)[1,])
})

test_that("splitSamChr works as expected", {
  fun    <- function(...) .Call(QuasR:::splitSamChr, ...)
  bamf1  <- pSingle@alignments$FileName[1]
  samf1  <- sub(".bam$", ".sam", bamf1)
  outdir <- tempfile(pattern = "sam", tmpdir = "extdata")
  dir.create(outdir)
  
  # arguments
  expect_error(fun(   1L, outdir))
  expect_error(fun(samf1,     1L))

  # results
  chrs <- fun(samf1, outdir)
  expect_identical(chrs, c("chrV", "splitChrSam_unaligned"))
  expect_identical(list.files(outdir), paste0(chrs, ".sam"))
  expect_length(readLines(file.path(outdir, "chrV.sam")), 514L)
  expect_length(readLines(file.path(outdir, "splitChrSam_unaligned.sam")), 2L)
  
  unlink(outdir, recursive = TRUE, force = TRUE)
})

test_that("countJunctions works as expected", {
  fun    <- function(...) .Call(QuasR:::countJunctions, ...)
  bamf1  <- pSingleAllelic@alignments$FileName[1]
  samf1  <- sub(".bam$", ".sam", bamf1)
  
  # arguments
  expect_error(fun(   1L, 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, "", 0L, 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, "", 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L,    "", FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L,    "", FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE,    "", 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, -1, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, 0L,   -1))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, 2L,   1L))
  expect_error(fun("err", 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(samf1, 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L))
  
  # results
  r1 <- fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L)
  r2 <- fun(bamf1, 0L, 0L, 1000L, TRUE,  FALSE, 0L, 255L)
  expect_is(r1, "integer")
  expect_length(r1, 512L)
  expect_is(r2, "list")
  expect_length(r2, 4L)
  expect_true(all(lengths(r2) == 512L))
})

test_that("bamfileToWig works as expected", {
  fun    <- function(...) .Call(QuasR:::bamfileToWig, ...)
  bamf1  <- pSingle@alignments$FileName[1]
  wig1   <- tempfile(fileext = ".wig",    tmpdir = "extdata")
  wig2   <- tempfile(fileext = ".wig.gz", tmpdir = "extdata")
  wig3   <- tempfile(fileext = ".wig",    tmpdir = "extdata")
  wig4   <- tempfile(fileext = ".wig.gz", tmpdir = "extdata")
  
  # arguments
  expect_error(fun(   1L, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1,   1L, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1,    "", 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE,  "", 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, "", "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L,  1L, 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*",  "", "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0,          1L, FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname",    1L,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   1L,          FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0",    1L, FALSE, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE,    1L, 0L, 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, "", 255L, 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L,   "", 0L, 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, "", 1000L, 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L,    "", 16L))
  expect_error(fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L,  ""))
  expect_error(fun("err", wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
                   "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 1000L, 16L))

  
  # results
  r1 <- fun(bamf1, wig1, FALSE, 10L, 0L, "*", 1.0, "trackname", FALSE,
            "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 800L, 16L)
  r2 <- fun(bamf1, wig3, FALSE, 20L, 0L, "*", 1.0, "trackname", TRUE,
            "255,212,0", FALSE, FALSE, 0L, 255L, 0L, 800L, 16L)
  r3 <- fun(bamf1, wig2, FALSE, 15L, 0L, "*", 1.0, "trackname", FALSE,
            "255,212,0", TRUE,  FALSE, 0L, 255L, 0L, 800L, 16L)
  r4 <- fun(bamf1, wig4, FALSE, 30L, 0L, "*", 1.0, "trackname", TRUE,
            "255,212,0", TRUE,  FALSE, 0L, 255L, 0L, 800L, 16L)
  l1 <- readLines(wig1)
  l2 <- readLines(wig2)
  l3 <- readLines(wig3)
  l4 <- readLines(wig4)
  expect_length(l1, 82L)
  expect_length(l2, 55L)
  expect_length(l3, 42L)
  expect_length(l4, 28L)
  expect_identical(sum(as.numeric(l1[-c(1, 2)])), sum(as.numeric(l2[-c(1, 2)])))
  expect_equal(sum(as.numeric(l3[-c(1, 2)])), 113.45)
  expect_equal(sum(as.numeric(l4[-c(1, 2)])), 85.73)
})

test_that("mergeReorderSam works as expected", {
  fun    <- function(...) .Call(QuasR:::mergeReorderSam, ...)
  bamf1  <- pPaired@alignments$FileName[1]
  samf1  <- sub(".bam$", ".sam", bamf1)
  samf2  <- tempfile(fileext = ".sam", tmpdir = "extdata")
  
  # arguments
  expect_error(fun(   1L, samf2,     0L, 1L))
  expect_error(fun(samf1,    1L,     0L, 1L))
  expect_error(fun(samf1, samf2,     "", 1L))
  expect_error(fun(samf1, samf2,     0L, ""))
  expect_error(fun(samf1, samf2,     4L, 1L))
  expect_error(fun(samf1, samf2,     1L, 1L))
  expect_error(fun(samf1, samf2,     2L, 1L))
  expect_error(fun(samf1, "err/err", 0L, 1L))

  # results
  expect_identical(fun(rep(samf1, 2), samf2, 0L, 1L), 1L)
  expect_length(readLines(samf2), 44L)
  expect_identical(fun(rep(samf1, 2), samf2, 2L, 1L), 1L)
  expect_length(readLines(samf2), 24L)
})

test_that("removeUnmappedFromSamAndConvertToBam works as expected", {
  fun    <- function(...) .Call(QuasR:::removeUnmappedFromSamAndConvertToBam, ...)
  bamf1  <- pPaired@alignments$FileName[1]
  samf1  <- sub(".bam$", ".sam", bamf1)
  bamf2  <- tempfile(fileext = ".bam", tmpdir = "extdata")
  
  # arguments
  expect_error(fun(   1L, bamf2))
  expect_error(fun(samf1,    1L))
  expect_error(fun("err", bamf2))
  expect_error(fun(samf1, "err/err"))

  # results
  expect_identical(fun(samf1, bamf2), bamf2)
})


test_that("filterHisat2 works as expected", {
  # create example sam file
  fun    <- function(...) .Call(QuasR:::filterHisat2, ...)
  samf1 <- tempfile(fileext = ".sam", tmpdir = "extdata")
  samf2 <- tempfile(fileext = ".sam", tmpdir = "extdata")
  writeLines(text = c("@HD\tVN:1.4",
                      "@SQ\tSN:chr1\tLN:40000",
                      "@SQ\tSN:chr2\tLN:10000",
                      "@SQ\tSN:chr3\tLN:45000",
                      "seq1\t163\tchr3\t12355\t255\t38M560N12M\t=\t13074\t769\tCAGCCCTTGAACGGAGAATAGAGTACATTGAAGCTCGGGTGACAAAAGGT\tggggggggggggggggggggggggggggghggggfggegccceacLSSTS\tNH:i:1",
                      "seq1\t83\tchr3\t13074\t255\t50M\t=\t12355\t-769\tAGAGCAACAGGGCTTATTCTTGTTTTTCTTTTTTCAAAAGTGTGGCCTTT\tgggggggggggggggggggggggfffffffffgggggggggggggggggg\tNH:i:1",
                      "seq2\t163\tchr2\t2927\t255\t14M241N36M\t=\t3424\t547\tCGGCGCTCGGCAAGTTCTCCCAGGAGAAAGCCATGTTCAGTTCGAGCGCC\tgggggggggghgggfggggggdggggggggggggggggggefeggdgggg\tNH:i:2",
                      "seq2\t83\tchr2\t3424\t255\t50M\t=\t2927\t-547\tGATGAACTCGGACCTCAAGGCTCAGCTCAGGGAGCTGAATATTACGGCAG\tgggggggggggggfggggggggggggfgggfggggggggggggggggghg\tNH:i:2",
                      "seq3\t163\tchr2\t5346\t255\t19M2335N31M\t=\t7796\t2500\tAGCAAAAGCGTCCCAGGAGCCGTACTCTGACAGCTGTGCACGATGCCATC\tgggggggggggggggggfggggggggggggggggdggggggggggegggg\tNH:i:3",
                      "seq3\t83\tchr2\t7796\t255\t50M\t=\t5346\t-2500\tCCGGCTCATAAAGGTTCATTTGGACAAAGCACAGCAGAACAATGTGGAAC\tgggggggggggggggggggggggggggggggggggggggggggggggggg\tNH:i:3",
                      "seq4\t99\tchr3\t2413\t255\t50M\t=\t12343\t9980\tCGGGAGATTCACCAGGACTGGGCTGACCAGGAGTACATTGAGATAATCAC\teggegddceedggggeZcec^]`^`ffffcTTSUTTSSTTfffaf]``]^\tNH:i:4",
                      "seq4\t147\tchr3\t12343\t255\t50M\t=\t2413\t-9980\tACGAGAAATTGACAGCCCTTGAACGGAGAATAGAGTACATTGAAGCTCGG\tacaaccaeccTTTSS_]_]`__P__SSTTTgdgggghgbgTTTTSggggg\tNH:i:4",
                      "seq5\t419\tchr1\t23628\t255\t50M\t=\t24226\t648\tCTCGTCCACTTTGAGTTCCTCGTTGAGCCTGATGGCGTCGGCAACCTCCT\tggfegfeeffTSSTSd`bbcc^c```^cc[aeefa_TT[ZTTTTTZ_[ZZ\tNH:i:1",
                      "seq5\t339\tchr1\t24226\t255\t50M\t=\t23628\t-648\tTCCACGGCGCGGAAGTGTGTCTTGCTCTCCTCCATGGCCTCCTGGAAGTG\tgggggggfgfW`]``ggggggfdgd^^c``[d_cbgdggggdgeggggb^\tNH:i:1",
                      "seq6\t77\t*\t0\t0\t*\t*\t0\t0\tGCCGCCCTCGCAGTGCTGCCAGAGAAGGGAGGGCATCCCCTGAGCCGCCG\tggggggggggggggggggggggggggggggggfffggggggggggdgggg",
                      "seq6\t141\t*\t0\t0\t*\t*\t0\t0\tGGCCCCAAGGCCCCCGTCCCGCAGCCACGCTTGTGGTCGCTGCGTCCCGG\tgggggggggggggggggggggggggggdgcgggedbe`^bEbV__ed^bc"),
             con = samf1)

  # arguments
  expect_error(fun(   1L, samf2,     1L))
  expect_error(fun(samf1,    1L,     1L))
  expect_error(fun(samf1, samf2,    "1L"))
  expect_error(fun("err", samf2,     1L))

  # results
  r0 <- fun(samf1, samf2, 0L); n0 <- length(readLines(samf2))
  r1 <- fun(samf1, samf2, 1L); n1 <- length(readLines(samf2))
  r2 <- fun(samf1, samf2, 2L); n2 <- length(readLines(samf2))
  r3 <- fun(samf1, samf2, 3L); n3 <- length(readLines(samf2))
  r4 <- fun(samf1, samf2, 4L); n4 <- length(readLines(samf2))

  expect_identical(n0, 14L)
  expect_identical(n1, 14L)
  expect_identical(n2, 14L)
  expect_identical(n3, 14L)
  expect_identical(n4, 14L)
  
  expect_identical(r0, c(n_secondary=2L, n_overmapped=8L))
  expect_identical(r1, c(n_secondary=2L, n_overmapped=6L))
  expect_identical(r2, c(n_secondary=2L, n_overmapped=4L))
  expect_identical(r3, c(n_secondary=2L, n_overmapped=2L))
  expect_identical(r4, c(n_secondary=2L, n_overmapped=0L)) 
})
