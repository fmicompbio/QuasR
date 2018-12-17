context("qProject")

# create qProject instance
genomeFile <- file.path("extdata", "hg19sub.fa")
sampleFile <- file.path("extdata", "samples_rna_single.txt")
auxFile <- file.path("extdata", "auxiliaries.txt")
projectGenome <- qAlign(sampleFile, genomeFile, projectName = "genome", clObj = clObj)
projectGenomeAux <- qAlign(sampleFile, genomeFile, auxiliaryFile = auxFile,
                           projectName = "genome-aux", clObj = clObj)

test_that("qProject length works as expected", {
  expect_identical(length(projectGenome), 4L)
  expect_identical(length(projectGenomeAux), 4L)
})

test_that("qProject has a working genome method", {
  expect_identical(normalizePath(genomeFile), normalizePath(genome(projectGenome)))
  expect_identical(normalizePath(genomeFile), normalizePath(genome(projectGenomeAux)))
})

test_that("qProject has a working auxiliaries method", {
  expect_true(0 == nrow(auxiliaries(projectGenome)))
  expect_true(1 == nrow(auxiliaries(projectGenomeAux)))
})

test_that("qProject has a working alignments method", {
  aln <- alignments(projectGenome)
  expect_true(4 == nrow(aln$genome))
  expect_true(0 == nrow(aln$aux))
  
  aln <- alignments(projectGenomeAux)
  expect_true(4 == nrow(aln$genome))
  expect_true(4 == ncol(aln$aux))
})

test_that("qProject has a working show method", {
  expect_output(show(projectGenome), "Project: genome")
  expect_output(show(projectGenomeAux), "Project: genome-aux")
})

test_that("qProject can be subset", {
  expect_is(projectGenome[1:2], "qProject")
  expect_identical(length(projectGenome[3:4]), 2L)
  expect_warning(projectGenome["Sample1"])
  expect_identical(projectGenome[1], suppressWarnings(projectGenome["Sample1"]))
  expect_identical(projectGenome[3], suppressWarnings(projectGenome["Sample2"]))
  expect_is(projectGenomeAux[1:2], "qProject")
  expect_identical(length(projectGenomeAux[3:4]), 2L)
})
