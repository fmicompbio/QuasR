
context("createAlignments-functions")

samf  <- system.file("unitTests", "cases", "ex1.sam.gz", package = "Rsamtools")
bout1 <- tempfile(tmpdir = "extdata")
bamf  <- Rsamtools::asBam(file = samf, destination = bout1)
bout2 <- tempfile(tmpdir = "extdata")

test_that("samToSortedBamParallel correctly digests its arguments", {
    expect_error(QuasR:::samToSortedBamParallel(file = "nonexistent.sam", destination = bout2, p = 2L, cacheDir = NULL))
    expect_error(QuasR:::samToSortedBamParallel(file = samf,              destination = bout2, p = 2L, cacheDir = "error"))
})

test_that("samToSortedBamParallel works as expected", {
    bamout <- QuasR:::samToSortedBamParallel(file = samf, destination = bout2, p = 2L, cacheDir = NULL)
    expect_identical(bamout, paste0(bout2, ".bam"))
    expect_identical(unname(alignmentStats(bamout)),
                     unname(alignmentStats(bamf)))
})
