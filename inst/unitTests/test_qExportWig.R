test_exportWig <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    library(rtracklayer)
    td <- tempdir()
    genomeFile <- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    sampleFile <- system.file(package="QuasR", "extdata", "samples_chip_single.txt")
    
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    wigfiles <- tempfile(fileext=rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling=F)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=F)))
    res <- qCount(project, wig[[1]], collapseBySample=F)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]))
    checkTrue(all(mcols(wig[[2]])$score == res[,3]))
    
    # paired and halfInsert shift
    sampleFile <- system.file(package="QuasR", "extdata", "samples_rna_paired.txt")
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    wigfiles <- tempfile(fileext=rep(".wig",2))
    res <- qExportWig(project, wigfiles, scaling=F)
    wig <- lapply(wigfiles, function(wf) suppressWarnings(import.wig(wf, asRangedData=F)))
    res <- qCount(project, wig[[1]], shift="halfInsert", collapseBySample=F)
    checkTrue(all(mcols(wig[[1]])$score == res[,2]/2))
}
