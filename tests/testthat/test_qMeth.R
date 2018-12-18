
context("qMeth")

test_that("qMeth correctly digests its arguments", {
  expect_error(qMeth("error"))
  expect_error(qMeth(pChipSingle))
  expect_error(qMeth(pBis, NULL, reportLevel = "alignment"))
  expect_error(qMeth(pBis, "error"))
  expect_error(qMeth(pBis, keepZero = "error"))
  expect_error(qMeth(pBis, asGRanges = "error"))
  expect_error(qMeth(pBis[c(1,1)], keepZero = FALSE, collapseBySample = FALSE))
  expect_error(qMeth(pBis, mode = "var", collapseByQueryRegion = TRUE))
  expect_error(qMeth(pBis, reportLevel = "alignment", mode = "var"))
  expect_error(qMeth(pBis, gtfGr, reportLevel = "alignment"))
  expect_error(qMeth(pBis, mask = gtfGr))
})

test_that("qMeth works as expected in (un-)directional mode", {
  requireNamespace("GenomicRanges")
  meth <- qMeth(pBisUndir, clObj = clObj)
  expect_identical(colSums(as(GenomicRanges::mcols(meth), "matrix")),
                   c("Sample1_T" = 37467, "Sample1_M" = 31907))
  expect_equal(GenomicRanges::mcols(meth)[201:210, 2],
               c(0, 3, 5, 7, 8, 6, 5, 0, 0, 2))

  meth <- qMeth(pBis)
  expect_identical(colSums(as(GenomicRanges::mcols(meth), "matrix")),
                   c("Sample1_T" = 37467, "Sample1_M" = 31907))
  expect_equal(GenomicRanges::mcols(meth)[201:210, 2],
               c(0, 3, 5, 7, 8, 6, 5, 0, 0, 2))
  
  meth <- qMeth(pBisSnps)
  expect_length(meth, 3013L)
  expect_identical(ncol(GenomicRanges::mcols(meth)), 6L)
})

test_that("qMeth correctly works with different 'mode' and 'reportLevel' arguments", {
  requireNamespace("GenomicRanges")
  
  meth <- qMeth(pBis, mode = "allC")
  expect_length(meth, 50969L)

  meth <- qMeth(pBis, mode = "CpG", collapseBySample = FALSE)
  expect_length(meth, 6220L)

  meth <- qMeth(pBis, mode = "CpGcomb")
  expect_length(meth, 3110L)

  meth <- qMeth(pBis, GenomicRanges::GRanges("chr1", IRanges(start = c(100, 300), width = 50)), collapseByQueryRegion = TRUE)
  expect_length(meth, 2L)
  
  meth <- qMeth(pBis, mode = "var")
  expect_length(meth, 6220L)
  
  gr <- GenomicRanges::GRanges("chr1", IRanges(start = 100, width = 50))
  meth <- qMeth(pBis, gr, mode = "CpG", reportLevel = "alignment")
  expect_is(meth, "list")
  expect_length(meth[[1]], 4L)
  expect_equal(meth[[1]]$meth, c(1,1,1,0,1,1,0,1,1,0))

  meth2 <- qMeth(pBis, gr, mode = "CpG", reportLevel = "alignment", collapseByQueryRegion = TRUE)
  expect_identical(meth, meth2)
})
