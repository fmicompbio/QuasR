
context("qCount")

test_that("qCount correctly digests its arguments", {
  expect_error(qCount(pSingle, query, mapqMin = -1))
  expect_error(qCount(pSingle, query, mapqMax = 256))
})

test_that("qCount correctly uses MAPQ", {
  aln <- Rsamtools::scanBam(pSingle@alignments$FileName[1], what = "mapq")
  binbreaks <- c(seq(0, 240, by = 40), 256)
  resSoll <- table(findInterval(aln[[1]]$mapq, vec = binbreaks,
                                all.inside = TRUE, rightmost.closed = TRUE))
  res <- unlist(lapply(seq_along(binbreaks[-1]),
                       function(i) qCount(pSingle[1], qGenome,
                                          mapqMin = binbreaks[i], mapqMax = binbreaks[i + 1] - 1)[1,2]))
  expect_equal(as.vector(resSoll), res)
})

test_that("qCount correctly uses shift and selectReadPosition", {
  # smart shift
  #fr R1->left R2->right
  res <- qCount(pPaired, q01.20, selectReadPosition = "start", shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 4 ,0), c(9, 1, 10)), unname(res))
  res <- qCount(pPaired, q01.20, selectReadPosition = "end",   shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(5, 1, 7, 1, 6)), unname(res))
  #fr R2->left R1->right
  res <- qCount(pPaired, q21.40, selectReadPosition = "start", shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 4 ,0), c(9, 1, 10)), unname(res))
  res <- qCount(pPaired, q21.40, selectReadPosition = "end",  shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(5, 1, 7, 1, 6)), unname(res))
  #ff R1->left R2->right
  res <- qCount(pPaired, q41.60, selectReadPosition = "start", shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(5, 1, 3, 1, 10)), unname(res))
  res <- qCount(pPaired, q41.60, selectReadPosition = "end",   shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(9, 1, 3, 1, 6)), unname(res))
  #rr R2->left R1->right
  res <- qCount(pPaired, q61.80, selectReadPosition = "start", shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(9, 1, 3, 1, 6)), unname(res))
  res <- qCount(pPaired, q61.80, selectReadPosition = "end",   shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(5, 1, 3, 1, 10)), unname(res))
  #rf R1->left R2->right
  res <- qCount(pPaired, q81.99, selectReadPosition = "start", shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 2, 0, 2, 0), c(5, 1, 7, 1, 5)), unname(res))
  res <- qCount(pPaired, q81.99, selectReadPosition = "end",   shift = "halfInsert", orientation = "any")[,-1]
  expect_equal(rep(c(0, 4 ,0), c(9, 1, 9)), unname(res))
  
  # qCount with interger as shift
  aln <- GenomicAlignments::readGAlignments(pPaired@alignments$FileName)
  resSoll <- rep(0, 99)
  pos <- Rle(ifelse(strand(aln) == "+", start(aln), end(aln)))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "start", shift = 0, orientation = "any")[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll <- rep(0, 99)
  pos <- Rle(ifelse(strand(aln) == "+", end(aln) + 1, start(aln) - 1))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "end", shift = 1, orientation = "any")[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll <- rep(0, 99)
  pos <- Rle(ifelse(strand(aln) == "+", start(aln) - 1, end(aln) + 1))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "start", shift = -1, orientation = "any")[,-1]
  expect_equal(resSoll, unname(res))
})

test_that("qCount correctly processes alignment orientation", {
  cnt1 <- qCount(pSingle, qGenome, orientation = "any",      collapseBySample = FALSE)
  cnt2 <- qCount(pSingle, qGenome, orientation = "same",     collapseBySample = FALSE)
  cnt3 <- qCount(pSingle, qGenome, orientation = "opposite", collapseBySample = FALSE)
  expect_equal(unname(cnt1[1,]), c(800, 512, 512))
  expect_equal(unname(cnt2[1,]), c(800, 512,   0))
  expect_equal(unname(cnt3[1,]), c(800,   0, 512))
})

test_that("qCount correctly uses useRead", {
  aln <- GenomicAlignments::readGAlignmentPairs(pPaired@alignments$FileName)
  aln1 <- GenomicAlignments::first(aln)
  aln2 <- GenomicAlignments::last(aln)

  resSoll <- rep(0, 99)
  pos <- Rle(ifelse(strand(aln1) == "+", start(aln1), end(aln1)))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "start", shift = 0, orientation = "any", useRead = "first")[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll <- rep(0, 99)
  pos <- Rle(ifelse(strand(aln1) == "+", end(aln1), start(aln1)))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "end", shift = 0, orientation = "any", useRead = "first")[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll <- rep(0,99)
  pos <- Rle(ifelse(strand(aln2) == "+", start(aln2), end(aln2)))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "start", shift = 0, orientation = "any", useRead = "last")[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll <- rep(0,99)
  pos <- Rle(ifelse(strand(aln2) == "+", end(aln2), start(aln2)))
  resSoll[runValue(pos)] <- runLength(pos)
  res <- qCount(pPaired, q01.99, selectReadPosition = "end", shift = 0, orientation = "any", useRead = "last")[,-1]
  expect_equal(resSoll, unname(res))
})

test_that("qCount corretly uses maxInsertSize", {
  resSoll <- rep(0,20)
  res <- qCount(pPaired, q01.20, selectReadPosition = "start", shift = "halfInsert", maxInsertSize = 0)[,-1]
  expect_equal(resSoll, unname(res))
  
  resSoll[10] <- 4
  res <- qCount(pPaired, q01.20, selectReadPosition = "start", shift = "halfInsert", maxInsertSize = 14)[,-1]
  expect_equal(resSoll, unname(res))
})

test_that("qCount correctly uses mask", {
  mask <- qTiles[names(qTiles) == "H4"]
  region <- qTiles
  strand(region) <- "+"    
  resSoll <- matrix(0, nrow = 4, ncol = 3) 
  resSoll[, 1] <- c(200, 300, 150, 0)
  resSoll[, 2] <- resSoll[, 3] <- c(150, 200, 100, 0)
  res <- qCount(pSingle, region, mask = mask, collapseBySample = FALSE, orientation = "any")
  expect_equal(resSoll, unname(res))
})

test_that("qCount correctly works with a GRangesList query", {
  regionList <- split(qTiles, names(qTiles))
  resSoll <- cbind(c(300, 150, 150, 0), c(212, 100, 100, 0), c(221, 100, 100, 0))
  res <- qCount(pSingle, regionList, collapseBySample = FALSE, orientation = "same")
  expect_equal(resSoll, unname(res))

  strand(regionList) <- "+"
  resSoll[, 3] <- 0
  res <- qCount(pSingle, regionList, collapseBySample = FALSE, orientation = "same")
  expect_equal(resSoll, unname(res))
  
  strand(regionList) <- "-"
  resSoll <- cbind(c(300, 150, 150, 0), c(0, 0, 0, 0), c(221, 100, 100, 0))
  res <- qCount(pSingle, regionList, collapseBySample = FALSE, orientation = "same")
  expect_equal(resSoll, unname(res))
})

test_that("qCount correctly works with a TxDb query", {
  requireNamespace("GenomicFeatures")
  requireNamespace("GenomicRanges")
  
  # gene
  res1 <- qCount(pRnaSingleSpliced, gtfGr, collapseBySample = FALSE, reportLevel = NULL)
  res2 <- qCount(pRnaSingleSpliced, txdb,  collapseBySample = FALSE, reportLevel = "gene")
  expect_identical(res1, res2)

  # exon
  exs <- GenomicFeatures::exons(txdb)
  res1 <- qCount(pRnaSingleSpliced, exs,  collapseBySample = FALSE)
  res2 <- qCount(pRnaSingleSpliced, txdb, collapseBySample = FALSE, reportLevel = "exon")
  expect_identical(res1, res2[rownames(res1),])

  # promoter
  proms <- GenomicFeatures::promoters(txdb)
  res1 <- qCount(pRnaSingleSpliced, proms, collapseBySample = FALSE, reportLevel = "promoter")
  res2 <- qCount(pRnaSingleSpliced, txdb,  collapseBySample = FALSE, reportLevel = "promoter")
  rownames(res2) <- sub("^[0-9]+;", "", rownames(res2))
  expect_identical(res1, res2[rownames(res1),])

  # junction, includeSpliced
  exGr <- GenomicRanges::GRanges("chr1", IRanges(start = c(11720, 12322, 14043, 14363),
                                                 end   = c(12212, 12518, 14165, 14512)))
  inGr <- GenomicRanges::GRanges("chr1", IRanges(start = c(12213, 14166),
                                                 end   = c(12321, 14362)), strand = "+")
  res1 <- qCount(pRnaSingleSpliced, exGr, collapseBySample = FALSE)
  res2 <- qCount(pRnaSingleSpliced, exGr, collapseBySample = FALSE, includeSpliced = FALSE)
  res3 <- qCount(pRnaSingleSpliced, NULL, reportLevel = "junction", collapseBySample = FALSE)
  expect_equal(unname((res1 - res2)[c(1, 3), -1]),
               unname(as.matrix(mcols(res3)[match(inGr, res3),])))
})

test_that("qCount correctly works in allelic mode", {
  region <- qTiles
  strand(region) <- "-"
  resSoll <- matrix(c(300, 300, 300, 250, rep(rep(0, 4), 3),
                      rep(c(221, 200, 200, 171), 3)), nrow = 4)
  res <- qCount(pSingleAllelic, region, collapseBySample = FALSE, orientation = "same")
  expect_identical(resSoll, unname(res))
})

test_that("qCount correctly works with auxiliary alignments", {
  res <- qCount(pChipSingleAux, auxGr, collapseBySample = FALSE, auxiliaryName = "phiX174")
  expect_equal(c(5386, 251, 493), unname(res[1,]))

  # include secondary alignments
  expect_identical(qCount(pPhiX, auxGr, useRead = "any"  )[1,"test"], 696)
  expect_identical(qCount(pPhiX, auxGr, useRead = "first")[1,"test"], 348)
  expect_identical(qCount(pPhiX, auxGr, useRead = "last" )[1,"test"], 348)
  # exclude secondary alignments
  expect_identical(qCount(pPhiX, auxGr, useRead = "any"  , includeSecondary = FALSE)[1,'test'], 384)
  expect_identical(qCount(pPhiX, auxGr, useRead = "first", includeSecondary = FALSE)[1,'test'], 192)
  expect_identical(qCount(pPhiX, auxGr, useRead = "last" , includeSecondary = FALSE)[1,'test'], 192)
})
