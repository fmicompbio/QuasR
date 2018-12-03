# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

test_qExportWig_qQCReport <- function()
{
    # generate test results
    clObj <- makeCluster(2)

    library(rtracklayer)
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    sampleFile <- file.path("extdata", "samples_chip_single.txt")
    
    project <- qAlign(sampleFile, genomeFile, alignmentsDir = td, clObj = clObj)
    wigfiles <- tempfile(fileext = rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling = FALSE)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf)))
    res <- qCount(project, wig[[1]], collapseBySample = FALSE)
    RUnit::checkTrue(all(mcols(wig[[1]])$score == res[,2]))
    RUnit::checkTrue(all(mcols(wig[[2]])$score == res[,3]))
    
    # check qQCReport
    tmppdf <- tempfile(fileext = ".pdf")
    resqc <- qQCReport(project, pdfFilename = tmppdf, clObj = clObj)
    RUnit::checkTrue(is.list(resqc))
    RUnit::checkTrue(length(resqc) == 8L)
    RUnit::checkTrue(is.list(resqc$nuclByCycle) &&
                  is.matrix(resqc$nuclByCycle[[1]]) &&
                  is.matrix(resqc$nuclByCycle[[2]]))
    RUnit::checkTrue(all(dim(resqc$nuclByCycle) == c(5L, 36L)))
    RUnit::checkTrue(all(dim(resqc$nuclByCycle[[1]]) == c(5L, 36L)))
    RUnit::checkTrue(all(dim(resqc$nuclByCycle[[2]]) == c(5L, 36L)))
    unlink(tmppdf)
    
    # check qExportWig
    # ... arguments
    RUnit::checkException(qExportWig("error"))
    RUnit::checkException(qExportWig(project, collapseBySample = FALSE, strand = c("+", "-")))
    RUnit::checkException(qExportWig(project, file = "test.wig"))
    RUnit::checkException(qExportWig(project, createBigWig = TRUE, file = c("w1.wig","w2.wig")))
    RUnit::checkException(qExportWig(project, file = c("w1.wig.bz2","w2.wig.bz2")))
    RUnit::checkException(qExportWig(project, binsize = c(100, 200)))
    RUnit::checkException(qExportWig(project, binsize = -100))
    RUnit::checkException(qExportWig(project, shift = NA))
    RUnit::checkException(qExportWig(project, shift = c(10, 20, 30)))
    RUnit::checkException(qExportWig(project, scaling = c(1,2,3)))
    RUnit::checkException(qExportWig(project, scaling = TRUE, log2p1 = "yes"))
    RUnit::checkException(qExportWig(project, scaling = FALSE, colors = "gray", includeSecondary = "yes"))
    RUnit::checkException(qExportWig(project, absIsizeMin = 100))

    # ... paired and halfInsert shift
    sampleFile <- file.path("extdata", "samples_rna_paired.txt")
    project2 <- qAlign(sampleFile, genomeFile, alignmentsDir = td, clObj = clObj)
    wigfiles <- tempfile(fileext = rep(".wig",2))
    res <- qExportWig(project2, wigfiles, scaling = FALSE)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf)))
    res <- qCount(project2, wig[[1]], shift = "halfInsert", collapseBySample = FALSE)
    RUnit::checkTrue(all(mcols(wig[[1]])$score == res[,2]/2))

    # ... includeSecondary
    project3 <- qAlign(file.path("extdata", "phiX_paired_withSecondary_sampleFile.txt"),
                      file.path("extdata", "NC_001422.1.fa"), paired = "fr")
    query <- GRanges("phiX174", IRanges(start = 1, end = 5386))
    wigfiles <- tempfile(fileext = rep(".wig.gz",2))
    res1 <- qExportWig(project3, wigfiles[1], scaling = FALSE, includeSecondary = TRUE)
    wigcnt1 <- sum(as.numeric(readLines(res1)[-(1:2)]))
    res2 <- qExportWig(project3, wigfiles[1], scaling = FALSE, includeSecondary = FALSE)
    wigcnt2 <- sum(as.numeric(readLines(res2)[-(1:2)]))
    cnt1 <- qCount(project3, query, useRead = "first", includeSecondary = TRUE)[1,2]
    cnt2 <- qCount(project3, query, useRead = "first", includeSecondary = FALSE)[1,2]
    RUnit::checkTrue(wigcnt1 == cnt1)
    RUnit::checkTrue(wigcnt2 == cnt2)

    stopCluster(clObj)
}
