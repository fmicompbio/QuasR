context("utilites")

test_that("md5subsum works as expected", {
  expect_identical_md5subsum <- function(fname, expval) {
    fname <- system.file("extdata", fname, package = "QuasR")
    expect_identical(QuasR:::md5subsum(fname), expval)
  }

    expect_warning(QuasR:::md5subsum("nonexistent"))
  expect_identical_md5subsum("hg19sub.fa",      "cf9c426a33b0e99261f310c67f0df0b6")
  expect_identical_md5subsum("bis_1_1.fa.bz2",  "fb9bd7b28edc59b41757833c68ae94ff")
  set.seed(95874)
  expect_identical_md5subsum("bis_1_1.fa.bz2",  "fb9bd7b28edc59b41757833c68ae94ff")
  set.seed(948620)
  expect_identical_md5subsum("chip_1_1.fq.bz2", "653105d10a200f5663ceb174027e4eb9")
  set.seed(95)
  expect_identical_md5subsum("rna_1_1.fq.bz2",  "89282741f79188e8b8bfb517f606c035")
  expect_identical_md5subsum("NC_001422.1.fa",  "9d91fb2b59c4134ab1dc3249ed81fbce")
})

test_that("displayNames correctly digests its arguments", {
  expect_error(displayNames("error"))
})

test_that("freeDiskSpace runs", {
  expect_warning(QuasR:::freeDiskSpace("nonexistent"))
  expect_is(QuasR:::freeDiskSpace("."), "numeric")
})

test_that("truncString works as expected", {
  expect_identical(QuasR:::truncString("AAAAAA", w = 5), "...AA")
})

test_that("getListOfBiocParallelParam works as expected", {
  requireNamespace("BiocParallel")
  bpp <- QuasR:::getListOfBiocParallelParam(list(BiocParallel::SerialParam()))
  expect_is(bpp, "list")
  expect_length(bpp, 2L)
  expect_is(bpp[[1]], "BiocParallelParam")
  
  bpp2 <- QuasR:::getListOfBiocParallelParam(clObj)
  expect_error(QuasR:::getListOfBiocParallelParam(bpp2[2:1]))
})
