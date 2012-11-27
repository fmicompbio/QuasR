test_qProfile <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    if(!"gtfRegion" %in% ls(envir=.GlobalEnv)){
        gtfRegion <<- createGtfRegion()
    }    
    
    genomeFile <- system.file(package="QuasR", "extdata", "hg19sub.fa")
    sampleFile <- system.file(package="QuasR", "extdata", "samples_chip_single.txt")
    project <- qAlign(sampleFile, genomeFile, clObj=clObj)
    
    query <- resize(gtfRegion, fix="start", width=200)
    query <- query[!duplicated(query)]
    name <- names(query)
    
    # test
    names(query) <- NULL
    pr <- qProfile(project, query, upstream=0, downstream=199)
    cnt <- qCount(project, query)
    checkTrue(all(sum(cnt[,2]) == rowSums(pr[[2]])))                  
    checkTrue(all(sum(cnt[,3]) == rowSums(pr[[3]])))

    queryShift <- query
    start(queryShift) <- ifelse(strand(query) == "+", start(query)-50, start(query)+50)
    end(queryShift) <- ifelse(strand(query) == "+", end(query)-50, end(query)+50)
    pr <- qProfile(project, query, upstream=50, downstream=149)
    cnt <- qCount(project, queryShift)
    checkTrue(all(sum(cnt[,2]) == rowSums(pr[[2]])))                  
    checkTrue(all(sum(cnt[,3]) == rowSums(pr[[3]])))
    
    # test shift
    pr <- qProfile(project, query, upstream=0, downstream=199, shift=50)
    cnt <- qCount(project, query, shift=50)
    checkTrue(all(sum(cnt[,2]) == rowSums(pr[[2]])))                  
    checkTrue(all(sum(cnt[,3]) == rowSums(pr[[3]])))
    
    # test multiple profile
#     names(query) <- name
#     pr <- qProfile(project, query, upstream=50, downstream=149)
#     cnt <- qCount(project, query)
#     cnt <- cnt[sort(rownames(cnt)),]
#     checkTrue(all(cnt[,2] == rowSums(pr[[2]])))                  
#     checkTrue(all(cnt[,3] == rowSums(pr[[3]])))

    # test smart shift
    sampleFile <- system.file(package="QuasR", "extdata", "samples_rna_paired.txt")
    project <- qAlign(sampleFile, genomeFile, clObj=clObj)
    pr <- qProfile(project, query, upstream=0, downstream=199, shift="halfInsert")
    cnt <- qCount(project, query, shift="halfInsert")
    checkTrue(all(sum(cnt[,2]) == rowSums(pr[[2]])))                  
    checkTrue(all(sum(cnt[,3]) == rowSums(pr[[3]])))
}
