context("qProject")

# create qProject instance
genomeFile     <- file.path("extdata", "hg19sub.fa")
sChipSingle    <- file.path("extdata", "samples_chip_single.txt")
auxFile        <- file.path("extdata", "auxiliaries.txt")
pChipSingle    <- qAlign(sChipSingle,  genomeFile, projectName = "genome")
pChipSingleAux <- qAlign(sChipSingle,  genomeFile, auxiliaryFile = auxFile,
                         projectName = "genome-aux")

test_that("qProject length works as expected", {
  expect_identical(length(pChipSingle), 2L)
  expect_identical(length(pChipSingleAux), 2L)
})

test_that("qProject has a working genome method", {
  expect_identical(normalizePath(genomeFile), normalizePath(genome(pChipSingle)))
  expect_identical(normalizePath(genomeFile), normalizePath(genome(pChipSingleAux)))
})

test_that("qProject has a working auxiliaries method", {
  expect_true(0 == nrow(auxiliaries(pChipSingle)))
  expect_true(1 == nrow(auxiliaries(pChipSingleAux)))
})

test_that("qProject has a working alignments method", {
  aln <- alignments(pChipSingle)
  expect_true(2 == nrow(aln$genome))
  expect_true(0 == nrow(aln$aux))
  
  aln <- alignments(pChipSingleAux)
  expect_true(2 == nrow(aln$genome))
  expect_true(2 == ncol(aln$aux))
})

test_that("qProject has a working show method", {
  expect_output(show(pChipSingle), "Project: genome")
  expect_output(show(pChipSingleAux), "Project: genome-aux")
})

test_that("qProject can be subset", {
  expect_is(pChipSingle[1:2], "qProject")
  expect_identical(length(pChipSingle[1]), 1L)
  expect_error(pChipSingle["Sample3"])
  expect_identical(pChipSingle[1], suppressWarnings(pChipSingle["Sample1"]))
  expect_identical(pChipSingle[2], suppressWarnings(pChipSingle["Sample2"]))
  expect_length(pChipSingleAux[2], 1L)
})
