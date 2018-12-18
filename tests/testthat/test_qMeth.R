
context("qMeth")

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
})

test_that("qMeth correctly works with different 'mode' arguments", {
  meth <- qMeth(pBis, mode = "allC")
  expect_length(meth, 50969L)

  meth <- qMeth(pBis, mode = "CpG")
  expect_length(meth, 6220L)

  meth <- qMeth(pBis, mode = "CpGcomb")
  expect_length(meth, 3110L)

  meth <- qMeth(pBis, mode = "var")
  expect_length(meth, 6220L)
})
