
context("qExportWig")

test_that("qExportWig correctly digests its arguments", {
  expect_error(qExportWig("error"))
  expect_error(qExportWig(pChipSingle, collapseBySample = FALSE, strand = c("+", "-")))
  expect_error(qExportWig(pChipSingle, file = "test.wig"))
  expect_error(qExportWig(pChipSingle, createBigWig = TRUE, file = c("w1.wig","w2.wig")))
  expect_error(qExportWig(pChipSingle, file = c("w1.wig.bz2","w2.wig.bz2")))
  expect_error(qExportWig(pChipSingle, binsize = c(100, 200)))
  expect_error(qExportWig(pChipSingle, binsize = -100))
  expect_error(qExportWig(pChipSingle, shift = NA))
  expect_error(qExportWig(pChipSingle, shift = c(10, 20, 30)))
  expect_error(suppressWarnings(qExportWig(pChipSingle, scaling = c(1,2,3))))
  expect_error(qExportWig(pChipSingle, scaling = TRUE, log2p1 = "yes"))
  expect_error(qExportWig(pChipSingle, scaling = FALSE, colors = "gray", includeSecondary = "yes"))
  expect_error(qExportWig(pChipSingle, absIsizeMin = 100))
  expect_error(qExportWig(pChipSingle, scaling = FALSE, mapqMin = "0"))
  expect_error(qExportWig(pChipSingle, scaling = FALSE, mapqMax = "0"))
  expect_warning(qExportWig(pSingle, "extdata/sample.wig", scaling = FALSE, useRead = "first"))
  expect_warning(qExportWig(pSingle, "extdata/sample.wig", scaling = FALSE, pairedAsSingle = TRUE))
})


test_that("qExportWig works as expected", {
  requireNamespace("rtracklayer")
  
  # paired and halfInsert shift
  wigfiles <- tempfile(fileext = rep(".wig",2), tmpdir = "extdata")
  expect_warning(res <- qExportWig(pRnaPaired, wigfiles, scaling = FALSE, shift = 1L))
  wig <- lapply(wigfiles, function(wf) suppressWarnings(rtracklayer::import.wig(wf)))
  res <- qCount(pRnaPaired, wig[[1]], shift = "halfInsert", collapseBySample = FALSE)
  expect_equal(mcols(wig[[1]])$score, unname(res[,2]/2))
  
  # check that we get the same files with parallelization
  wigfilespar <- tempfile(fileext = rep(".wig",2), tmpdir = "extdata")
  expect_warning(respar <- qExportWig(pRnaPaired, wigfilespar, scaling = FALSE, shift = 1L, clObj = clObj))
  expect_equal(unname(md5sum(wigfiles)), unname(md5sum(wigfilespar)))

  # includeSecondary, strand, scaling
  auxGrPlus <- auxGr; strand(auxGrPlus) <- "+"
  res1 <- qExportWig(pPhiX, wigfiles[1], scaling = FALSE, strand = "+", includeSecondary = TRUE)
  wigcnt1 <- sum(as.numeric(readLines(res1)[-c(1, 2)]))
  res2 <- qExportWig(pPhiX, wigfiles[1], scaling = 1e6, includeSecondary = FALSE)
  wigcnt2 <- sum(as.numeric(readLines(res2)[-c(1, 2)]))
  cnt1 <- qCount(pPhiX, auxGrPlus, useRead = "first", orientation = "same", includeSecondary = TRUE)[1,2]
  cnt2 <- qCount(pPhiX, auxGr, useRead = "first", includeSecondary = FALSE)[1,2]
  expect_identical(wigcnt1, cnt1)
  expect_equal(wigcnt2, cnt2 / 696 * 1e6, tolerance = 0.001)
  
  # bigwig
  bwfiles <- tempfile(fileext = rep(".bw", 2), tmpdir = "extdata")
  expect_identical(qExportWig(pPhiX, bwfiles[1], scaling = FALSE, useRead = "first", createBigWig = TRUE),
                   bwfiles[1])
  expect_identical(qExportWig(pPhiX, bwfiles[2], scaling = FALSE, useRead = "last",  createBigWig = TRUE),
                   bwfiles[2])
  
  # check that we get the same files with parallelization
  bwfilespar <- tempfile(fileext = rep(".bw",1), tmpdir = "extdata")
  resparbw <- qExportWig(pPhiX, bwfilespar, scaling = FALSE, useRead = "first", 
                         createBigWig = TRUE, clObj = clObj)
  expect_equal(unname(md5sum(bwfiles[1])), unname(md5sum(bwfilespar)))
})
