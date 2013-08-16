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
    res <- qExportWig(project, wigfiles, scaling=F)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=F)))
    res <- qCount(project, wig[[1]], collapseBySample=F)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]))
    checkTrue(all(mcols(wig[[2]])$score == res[,3]))
    
    # paired and halfInsert shift
    sampleFile <- file.path("extdata", "samples_rna_paired.txt")
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    wigfiles <- tempfile(fileext=rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling=F)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=F)))
    res <- qCount(project, wig[[1]], shift="halfInsert", collapseBySample=F)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]/2))

    stopCluster(clObj)
}
