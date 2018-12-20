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
})

test_that("countJunctions works as expected", {
  fun    <- function(...) .Call(QuasR:::countJunctions, ...)
  bamf1  <- pChipSingleSnps@alignments$FileName[1]
  
  # arguments
  expect_error(fun(   1L, 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, "", 0L, 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, "", 1000L, FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L,    "", FALSE, FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L,    "", FALSE, 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE,    "", 0L, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, -1, 255L))
  expect_error(fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, 0L,   -1))
  expect_error(fun("err", 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L))

  # results
  r1 <- fun(bamf1, 0L, 0L, 1000L, FALSE, FALSE, 0L, 255L)
  r2 <- fun(bamf1, 0L, 0L, 1000L, TRUE,  FALSE, 0L, 255L)
  expect_is(r1, "integer")
  expect_length(r1, 0L)
  expect_is(r2, "list")
  expect_length(r2, 4L)
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
  expect_identical(sum(as.numeric(l1[-(1:2)])), sum(as.numeric(l2[-(1:2)])))
  expect_equal(sum(as.numeric(l3[-(1:2)])), 113.45)
  expect_equal(sum(as.numeric(l4[-(1:2)])), 85.73)
})
