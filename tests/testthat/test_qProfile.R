
context("qProfile")

test_that("qProfile works as expected", {
  requireNamespace("GenomicRanges")
  qreg <- GenomicRanges::resize(gtfGr, fix = "start", width = 200)
  qreg <- qreg[!duplicated(qreg)]
  names(qreg) <- NULL
  
  # no shift
  pr  <- qProfile(pChipSingle, qreg, upstream = 0, downstream = 199)
  cnt <- qCount(pChipSingle, qreg)
  expect_equal(sum(cnt[,2]), unname(rowSums(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))

  # shift
  pr  <- qProfile(pChipSingle, qreg, upstream = 0, downstream = 199, shift = 50)
  cnt <- qCount(pChipSingle, qreg, shift = 50)
  expect_equal(sum(cnt[,2]), unname(rowSums(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))
  
  # smart shift
  pr  <- qProfile(pRnaPaired, qreg, upstream = 0, downstream = 199, shift = "halfInsert")
  cnt <- qCount(pRnaPaired, qreg, shift = "halfInsert")
  expect_equal(sum(cnt[,2]), unname(rowSums(pr[[2]])))
  expect_equal(sum(cnt[,3]), unname(rowSums(pr[[3]])))

  # test includeSecondary
  cnt1 <- qCount(pPhiX, auxGr, includeSecondary = TRUE)[1,2]
  cnt2 <- qCount(pPhiX, auxGr, includeSecondary = FALSE)[1,2]
  pr1  <- qProfile(pPhiX, auxGr, upstream = 0, downstream = 5386, includeSecondary = TRUE)
  pr2  <- qProfile(pPhiX, auxGr, upstream = 0, downstream = 5386, includeSecondary = FALSE)
  expect_equal(cnt1, sum(pr1[[2]]))
  expect_equal(cnt2, sum(pr2[[2]]))
})
