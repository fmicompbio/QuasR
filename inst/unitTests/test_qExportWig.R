# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

test_exportWig <- function()
{
    clObj <- makeCluster(2)

    library(rtracklayer)
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    sampleFile <- file.path("extdata", "samples_chip_single.txt")
    
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    wigfiles <- tempfile(fileext=rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling=FALSE)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=FALSE)))
    res <- qCount(project, wig[[1]], collapseBySample=FALSE)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]))
    checkTrue(all(mcols(wig[[2]])$score == res[,3]))
    
    # paired and halfInsert shift
    sampleFile <- file.path("extdata", "samples_rna_paired.txt")
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    wigfiles <- tempfile(fileext=rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling=FALSE)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=FALSE)))
    res <- qCount(project, wig[[1]], shift="halfInsert", collapseBySample=FALSE)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]/2))

    # includeSecondary
    project <- qAlign(file.path("extdata", "phiX_paired_withSecondary_sampleFile.txt"),
                      file.path("extdata", "NC_001422.1.fa"), paired="fr")
    query <- GRanges("phiX174", IRanges(start=1, end=5386))
    wigfiles <- tempfile(fileext=rep(".wig.gz",2))
    res1 <- qExportWig(project, wigfiles[1], scaling=FALSE, includeSecondary=TRUE)
    wigcnt1 <- sum(as.numeric(readLines(res1)[-(1:2)]))
    res2 <- qExportWig(project, wigfiles[1], scaling=FALSE, includeSecondary=FALSE)
    wigcnt2 <- sum(as.numeric(readLines(res2)[-(1:2)]))
    cnt1 <- qCount(project, query, useRead="first", includeSecondary=TRUE)[1,2]
    cnt2 <- qCount(project, query, useRead="first", includeSecondary=FALSE)[1,2]
    checkTrue(wigcnt1==cnt1)
    checkTrue(wigcnt2==cnt2)

    stopCluster(clObj)
}
