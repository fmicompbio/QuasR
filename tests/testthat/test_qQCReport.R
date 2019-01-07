
context("qQCReport")

test_that("qQCReport works as expected", {
  # arguments
  expect_error(qQCReport("error"))
  expect_error(qQCReport(c("f1.bam", "f2.fastq")))
  
  # results
  tmppdf <- tempfile(fileext = rep(".pdf", 6), tmpdir = "extdata")

  resqc  <- qQCReport(pChipSingle, pdfFilename = tmppdf[1])
  resqc2 <- qQCReport(pChipSingle, pdfFilename = tmppdf[2], clObj = clObj)
  resqc3 <- qQCReport(pRnaPaired,  pdfFilename = tmppdf[3], useSampleNames = FALSE)
  resqc4 <- qQCReport(alignments(pChipSingle)$genome$FileName[1],
                      pdfFilename = tmppdf[4])
  resqc5 <- qQCReport(pChipSingle@reads$FileName[1],
                      pdfFilename = tmppdf[5])

  expect_identical(resqc, resqc2)
  expect_is(resqc, "list")
  expect_length(resqc,  8L)
  expect_length(resqc3, 9L)
  expect_length(resqc4, 5L)
  expect_length(resqc5, 4L)
  expect_is(resqc3$fragDistribution, "matrix")
  expect_is(resqc$nuclByCycle, "list")
  expect_is(resqc$nuclByCycle[[1]], "matrix")
  expect_is(resqc$nuclByCycle[[2]], "matrix")
  expect_identical(dim(resqc$nuclByCycle[[1]]), c(5L, 36L))
  expect_identical(dim(resqc$nuclByCycle[[2]]), c(5L, 36L))

  unlink(tmppdf)
})

test_that("qQCReport helper functions work as expected", {
  requireNamespace("Rsamtools")
  
  samempty <- tempfile(fileext = ".sam", tmpdir = "extdata")
  writeLines(c("@HD\tVN:1.0\tSO:coordinate","@SQ\tSN:chrV\tLN:99"), con = samempty)
  bamempty <- Rsamtools::asBam(samempty)
  
  bamf1 <- pChipSingle@alignments$FileName[1]
  
  tmppdf <- tempfile(fileext = ".pdf", tmpdir = "extdata")

  # calcQaInformation
  expect_warning(r <- QuasR:::calcQaInformation("err.fa.gz", "err", "fasta", 1e6))
  expect_null(r)
  
  # calcMmInformation
  l <- QuasR:::calcMmInformation(bamempty, genomeFile, 1e6)
  expect_is(l, "list")
  expect_length(l, 4L)
  expect_true(all(is.na(l[[3]])))

  l <- QuasR:::calcMmInformation(bamf1, genomeFile, 10)
  expect_length(l, 4L)

  # truncStringToPlotWidth
  s <- c("012345",
         "01234567890",
         "01234567890012345")
  pdf(tmppdf, height = 5, width = 5)
  plot(0:1, 0:1)
  expect_identical(QuasR:::truncStringToPlotWidth(s[1], 0.1), s[1])
  expect_identical(QuasR:::truncStringToPlotWidth(s[2], 0.1), "012345...")
  expect_identical(QuasR:::truncStringToPlotWidth(s[3], 0.1), "012345...")
  dev.off()
  unlink(tmppdf)
})
