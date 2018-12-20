
context("qQCReport")

test_that("qQCReport correctly digests its arguments", {
  expect_error(qQCReport("error"))
  expect_error(qQCReport(c("f1.bam", "f2.fastq")))
})

test_that("qQCReport works as expected", {
  tmppdf <- tempfile(fileext = rep(".pdf", 6), tmpdir = "extdata")

  resqc  <- qQCReport(pChipSingle, pdfFilename = tmppdf[1])
  resqc2 <- qQCReport(pChipSingle, pdfFilename = tmppdf[2], clObj = clObj)
  resqc3 <- qQCReport(pRnaPaired,  pdfFilename = tmppdf[3], useSampleNames = FALSE)
  resqc4 <- qQCReport(alignments(pChipSingle)$genome$FileName[1],
                      pdfFilename = tmppdf[4])
  resqc5 <- qQCReport(pChipSingle@reads$FileName[1],
                      pdfFilename = tmppdf[5])

  pdf(tmppdf[6], height = 5, width = 5)
  QuasR:::plotDuplicated(resqc$raw$qa, lmat = rbind(1:2))
  dev.off()

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
