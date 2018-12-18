
context("qQCReport")

test_that("qQCReport works as expected", {
  tmppdf <- tempfile(fileext = rep(".pdf", 2), tmpdir = "extdata")
  resqc  <- qQCReport(pChipSingle, pdfFilename = tmppdf[1])
  resqc2 <- qQCReport(pChipSingle, pdfFilename = tmppdf[2], clObj = clObj)
  expect_identical(resqc, resqc2)
  expect_is(resqc, "list")
  expect_length(resqc, 8L)
  expect_is(resqc$nuclByCycle, "list")
  expect_is(resqc$nuclByCycle[[1]], "matrix")
  expect_is(resqc$nuclByCycle[[2]], "matrix")
  expect_identical(dim(resqc$nuclByCycle[[1]]), c(5L, 36L))
  expect_identical(dim(resqc$nuclByCycle[[2]]), c(5L, 36L))
  unlink(tmppdf)
})
