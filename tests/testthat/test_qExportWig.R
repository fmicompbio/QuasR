
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
})


test_that("qExportWig works as expected", {
  requireNamespace("rtracklayer")
  
  # paired and halfInsert shift
  wigfiles <- tempfile(fileext = rep(".wig",2))
  res <- qExportWig(pRnaPaired, wigfiles, scaling = FALSE)
  wig <- lapply(wigfiles, function(wf) suppressWarnings(rtracklayer::import.wig(wf)))
  res <- qCount(pRnaPaired, wig[[1]], shift = "halfInsert", collapseBySample = FALSE)
  expect_equal(mcols(wig[[1]])$score, unname(res[,2]/2))

  # includeSecondary
  res1 <- qExportWig(pPhiX, wigfiles[1], scaling = FALSE, includeSecondary = TRUE)
  wigcnt1 <- sum(as.numeric(readLines(res1)[-(1:2)]))
  res2 <- qExportWig(pPhiX, wigfiles[1], scaling = FALSE, includeSecondary = FALSE)
  wigcnt2 <- sum(as.numeric(readLines(res2)[-(1:2)]))
  cnt1 <- qCount(pPhiX, auxGr, useRead = "first", includeSecondary = TRUE)[1,2]
  cnt2 <- qCount(pPhiX, auxGr, useRead = "first", includeSecondary = FALSE)[1,2]
  expect_identical(wigcnt1, cnt1)
  expect_identical(wigcnt2, cnt2)
})
