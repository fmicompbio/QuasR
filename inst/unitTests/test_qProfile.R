# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!existsFunction("createFastaReads"))
    source(system.file(package="QuasR", "unitTests", "help_function.R"))

if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

gtfRegion <- createGtfRegion()
td <- tempdir()
genomeFile <- file.path("extdata", "hg19sub.fa")

.setUp <- function() { # runs before each test_...()
    # make sure clObj exists and is a working cluster object
    if(!exists("clObj", envir=.GlobalEnv) ||
       !inherits(clObj, "cluster") ||
       inherits(try(all(unlist(clusterEvalQ(clObj, TRUE))), silent=TRUE), "try-error")) {
        clObj <<- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
        clusterEvalQ(clObj, library("QuasR"))
    }
}


test_qProfile <- function()
{
    sampleFile <- file.path("extdata", "samples_chip_single.txt")
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    
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
    sampleFile <- file.path("extdata", "samples_rna_paired.txt")
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    pr <- qProfile(project, query, upstream=0, downstream=199, shift="halfInsert")
    cnt <- qCount(project, query, shift="halfInsert")
    checkTrue(all(sum(cnt[,2]) == rowSums(pr[[2]])))                  
    checkTrue(all(sum(cnt[,3]) == rowSums(pr[[3]])))

    # test includeSecondary
    project <- qAlign(file.path("extdata", "phiX_paired_withSecondary_sampleFile.txt"),
                      file.path("extdata", "NC_001422.1.fa"), paired="fr")
    query <- GRanges("phiX174", IRanges(start=1, end=5386))
    cnt1 <- qCount(project, query, includeSecondary=TRUE)[1,2]
    cnt2 <- qCount(project, query, includeSecondary=FALSE)[1,2]
    pr1 <- qProfile(project, query, upstream=0, downstream=5386, includeSecondary=TRUE)
    pr2 <- qProfile(project, query, upstream=0, downstream=5386, includeSecondary=FALSE)
    checkTrue(cnt1 == sum(pr1[[2]]))                  
    checkTrue(cnt2 == sum(pr2[[2]]))
}
