test_md5subsum <- function()
{
    # check if md5subsum is the same when called several times
    file <- system.file(package="QuasR", "extdata", "hg19sub.fa")
    sollRes <- "cf9c426a33b0e99261f310c67f0df0b6"
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)

    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)
    
    # check if md5subsum is the same when set is seed
    file <- system.file(package="QuasR", "extdata", "bis_1_1.fa.bz2")
    sollRes <- "fb9bd7b28edc59b41757833c68ae94ff"
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)
    
    set.seed(95874)
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)

    set.seed(948620)
    file <- system.file(package="QuasR", "extdata", "chip_1_1.fq.bz2")
    sollRes <- "653105d10a200f5663ceb174027e4eb9"
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)

    set.seed(95)
    file <- system.file(package="QuasR", "extdata", "rna_1_1.fq.bz2")
    sollRes <- "89282741f79188e8b8bfb517f606c035"
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)

    file <- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
    sollRes <- "9d91fb2b59c4134ab1dc3249ed81fbce"
    res <- QuasR:::md5subsum(file)
    checkTrue(sollRes == res)
}

test_alignmentStats <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    sampleFile <- system.file(package="QuasR", "extdata", "samples_chip_single.txt")
    
    resSoll <- matrix(c(95000,95000,
                        2339,3609,
                        258,505), nrow=2, ncol=3)
    project <- qAlign(sampleFile, genomeFile, clObj=clObj)
    res <- alignmentStats(project)
    checkTrue(all(resSoll == res))
    
    resSoll <- matrix(c(95000,95000,5386,5386,
                        2339,3609,251,493,
                        258,505,7,12), nrow=4, ncol=3)
    project <- qAlign(sampleFile, genomeFile, auxFile, clObj=clObj)
    res <- alignmentStats(project)
    checkTrue(all(resSoll == res))
    
    resSoll <- matrix(c(95000,95000,
                        2339,3609,
                        258,505), nrow=2, ncol=3)
    project <- qAlign(sampleFile, genomeFile, clObj=clObj)
    res <- alignmentStats(alignments(project)$genome$FileName)
    checkTrue(all(resSoll == res))
    
    resSoll <- matrix(c(95000,2339,258), nrow=1, ncol=3)
    project <- qAlign(sampleFile, genomeFile, clObj=clObj)
    res <- alignmentStats(alignments(project)$genome$FileName[1])
    checkTrue(all(resSoll == res))
}
