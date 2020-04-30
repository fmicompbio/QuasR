
context("qProfile")

test_that("qProfile correctly digests its arguments", {
  expect_error(qProfile("error", gtfGr))
  expect_error(qProfile(pSingle, "error"))
  expect_error(qProfile(pSingle, qGenome, mask = "error"))
  expect_error(qProfile(pSingle, qGenome, mask = qGenome))
  expect_error(qProfile(pChipSingle, gtfGr, upstream = "error"))
  expect_error(qProfile(pChipSingle, gtfGr, downstream = "error"))
  expect_error(qProfile(pChipSingle, gtfGr, shift = "halfInsert"))
  expect_error(qProfile(pChipSingle, gtfGr, shift = "error"))
  expect_error(qProfile(pChipSingle, gtfGr, collapseBySample = "error"))
  expect_error(qProfile(pChipSingle, gtfGr, includeSpliced = "error"))
  expect_error(qProfile(pChipSingle, gtfGr, includeSecondary = "error"))
  expect_error(qProfile(pSingle, gtfGr))
  expect_warning(qProfile(pSingle, qTiles, useRead = "first"))
  expect_error(qProfile(pSingle, qGenome, upstream = c(100, 100)))
  expect_error(qProfile(pSingle, qGenome, downstream = c(100, 100)))
  expect_error(qProfile(pSingle, qGenome, binSize = "error"))
  expect_warning(qProfile(pSingle, qGenome, upstream = 100001))
})

test_that("qProfile works as expected", {
  requireNamespace("GenomicRanges")
  qreg <- GenomicRanges::resize(gtfGr, fix = "start", width = 200)
  qreg <- qreg[!duplicated(qreg)]
  names(qreg) <- NULL
  
  # no shift
  pr  <- qProfile(pChipSingle, qreg, upstream = 0, downstream = 199, shift = c(0, 0), clObj = clObj)
  cnt <- qCount(pChipSingle, qreg)
  expect_equal(sum(cnt[,2]), unname(sum(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))

  # shift
  pr  <- qProfile(pChipSingle, qreg, upstream = 0, downstream = 199, shift = 50)
  cnt <- qCount(pChipSingle, qreg, shift = 50)
  expect_equal(sum(cnt[,2]), unname(rowSums(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))
  
  # smart shift, useRead
  pr  <- qProfile(pRnaPaired, qreg, upstream = 0, downstream = 199, shift = "halfInsert", useRead = "first")
  cnt <- qCount(pRnaPaired, qreg, shift = "halfInsert", useRead = "first")
  expect_equal(sum(cnt[,2]), unname(rowSums(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))

  # includeSecondary
  cnt1 <- qCount(pPhiX, auxGr, includeSecondary = TRUE)[1,2]
  cnt2 <- qCount(pPhiX, auxGr, includeSecondary = FALSE)[1,2]
  pr1  <- qProfile(pPhiX, auxGr, upstream = 0, downstream = 5386, includeSecondary = TRUE)
  pr2  <- qProfile(pPhiX, auxGr, upstream = 0, downstream = 5386, includeSecondary = FALSE)
  expect_equal(cnt1, sum(pr1[[2]]))
  expect_equal(cnt2, sum(pr2[[2]]))
  
  # allelic
  pr1 <- qProfile(pChipSingleSnps, qreg, collapseBySample = TRUE)
  pr2 <- qProfile(pChipSingleSnps, qreg, collapseBySample = FALSE)
  expect_identical(unname(pr1), unname(pr2))
  expect_length(pr1, 7L)
  
  # binSize
  pr1 <- qProfile(pChipSingle, qreg, upstream = 100, downstream = 100)
  pr2 <- qProfile(pChipSingle, qreg, upstream = 100, downstream = 100, binSize = 10)
  expect_identical(colnames(pr2[[1]]), as.character(seq(-100, 100, by = 10)))
  expect_identical(sapply(pr1[-1], sum), sapply(pr2[-1], sum))
})
