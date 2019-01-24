
context("createAlignments-functions")

bamf     <- pChipSingle@alignments$FileName[1]
samf     <- Rsamtools::asSam(file = bamf)
boutnext <- tempfile(tmpdir = "extdata")

test_that("samToSortedBamParallel correctly digests its arguments", {
    expect_error(QuasR:::samToSortedBamParallel(file = "nonexistent.sam", destination = boutnext, p = 2L, cacheDir = NULL))
    expect_error(QuasR:::samToSortedBamParallel(file = samf,              destination = boutnext, p = 2L, cacheDir = "error"))
})

test_that("samToSortedBamParallel works as expected", {
    cat("bamf: ", bamf, "\n")
    cat("samf: ", samf, "\n")
    cat("  first 10 lines:\n", paste(readLines(samf, n = 10L), collapse = "\n"), "\n")
    cat("bamout: ", bamout, "\n")
    bamout <- QuasR:::samToSortedBamParallel(file = samf, destination = boutnext, p = 2L, cacheDir = NULL)
    expect_identical(bamout, paste0(boutnext, ".bam"))
    expect_identical(unname(alignmentStats(bamout)),
                     unname(alignmentStats(bamf)))
})
