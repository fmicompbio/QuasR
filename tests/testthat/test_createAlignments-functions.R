
context("createAlignments-functions")

bamf     <- pChipSingle@alignments$FileName[1]
samf     <- Rsamtools::asSam(file = bamf)
boutnext <- tempfile(tmpdir = "extdata")

test_that("samToSortedBamParallel correctly digests its arguments", {
    expect_error(QuasR:::samToSortedBamParallel(file = "nonexistent.sam", destination = boutnext, p = 2L, cacheDir = NULL))
    expect_error(QuasR:::samToSortedBamParallel(file = samf,              destination = boutnext, p = 2L, cacheDir = "error"))
})

test_that("samToSortedBamParallel works as expected", {
    bamout <- QuasR:::samToSortedBamParallel(file = samf, destination = boutnext, p = 2L, cacheDir = NULL)
    message("bamf: ", bamf)
    message("samf: ", samf)
    message("bamout: ", bamout)
    expect_identical(bamout, paste0(boutnext, ".bam"))
    expect_identical(unname(alignmentStats(bamout)),
                     unname(alignmentStats(bamf)))
})
