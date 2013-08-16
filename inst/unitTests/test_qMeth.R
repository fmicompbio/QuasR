# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!existsFunction("createFastaReads"))
    source(system.file(package="QuasR", "unitTests", "help_function.R"))

if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

library(Rsamtools)

genomeFile <- file.path("extdata", "hg19sub.fa")
genome <- scanFa(genomeFile)
td <- tempdir()
unMethRanges <- scanFaIndex(genomeFile)
partialMethRanges <- narrow(unMethRanges, start=5000, end=c(15000,10000,15000))

sampleFileGenomeSingle <- createReads(genomeFile, td, paired=FALSE)
sampleFileGenomePaired <- createReads(genomeFile, td, paired=TRUE)

sampleFileGenomePairedBisPartial <- createReads(genomeFile, td, paired=TRUE, bisulfit=partialMethRanges)

sampleFileGenomeSingleBisUn <- createReads(genomeFile, td, paired=FALSE, bisulfit=unMethRanges)
sampleFileGenomePairedBisUn <- createReads(genomeFile, td, paired=TRUE, bisulfit=unMethRanges)



.setUp <- function() { # runs before each test_...()
    # make sure clObj exists and is a working cluster object
    if(!exists("clObj", envir=.GlobalEnv) ||
       !inherits(clObj, "cluster") ||
       inherits(try(all(unlist(clusterEvalQ(clObj, TRUE))), silent=TRUE), "try-error")) {
        clObj <<- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
        clusterEvalQ(clObj, library("QuasR"))
    }
}


test_un_meth_paired <- function(){
    project <- qAlign(sampleFileGenomePairedBisUn, genomeFile, bisulfit="dir", alignmentsDir=td, clObj=clObj)

    meth <- qMeth(project, mode="allC")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated directed, mode=allC")
    
    meth <- qMeth(project, mode="CpG")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated directed, mode=CpG")

    meth <- qMeth(project, mode="CpGcomb")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated directed, mode=CpGcomb")
    
    project <- qAlign(sampleFileGenomePairedBisUn, genomeFile, bisulfit="undir", alignmentsDir=td, clObj=clObj)
    
    meth <- qMeth(project, mode="allC")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=allC")
    
    meth <- qMeth(project, mode="CpG")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=CpG")
    
    meth <- qMeth(project, mode="CpGcomb")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=CpGcomb")
}

test_un_meth_single <- function(){
    project <- qAlign(sampleFileGenomeSingleBisUn, genomeFile, bisulfit="undir", alignmentsDir=td, clObj=clObj)
    
    meth <- qMeth(project, mode="allC")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=allC")
    
    meth <- qMeth(project, mode="CpG")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=CpG")
    
    meth <- qMeth(project, mode="CpGcomb")
    checkTrue(all(0 == mcols(meth)[,2]), "Test un-methylated undirected, mode=CpGcomb")    
}

test_full_meth_paired <- function(){
    project <- qAlign(sampleFileGenomePaired, genomeFile, bisulfit="dir", alignmentsDir=td, clObj=clObj)

    meth <- qMeth(project, mode="allC")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated directed, mode=allC")
    
    meth <- qMeth(project, mode="CpG")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated directed, mode=CpG")
    
    meth <- qMeth(project, mode="CpGcomb")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated directed, mode=CpGcomb")
    
    meth <- qMeth(project, mode="var")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated, mode=var")
    
    project <- qAlign(sampleFileGenomePaired, genomeFile, bisulfit="undir", alignmentsDir=td, clObj=clObj)
    
    meth <- qMeth(project, mode="allC")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=allC")
    
    meth <- qMeth(project, mode="CpG")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=CpG")
    
    meth <- qMeth(project, mode="CpGcomb")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=CpGcomb")
    
    meth <- qMeth(project, mode="var")
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated, mode=var")
}

test_full_meth_single <- function(){
    project <- qAlign(sampleFileGenomeSingle, genomeFile, bisulfit="undir", alignmentsDir=td, clObj=clObj)
    
    meth <- qMeth(project, mode="allC", keepZero=FALSE)
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=allC")
    
    meth <- qMeth(project, mode="CpG", keepZero=FALSE)
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=CpG")
    
    meth <- qMeth(project, mode="CpGcomb", keepZero=FALSE)
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated undirected, mode=CpGcomb")
    
    meth <- qMeth(project, mode="var", keepZero=FALSE)
    checkTrue(all(mcols(meth)[,1] == mcols(meth)[,2]), "Test full-methylated, mode=var")
}

test_partial_meth <- function(){
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfit="dir", alignmentsDir=td, clObj=clObj)
    
    methAllC <- qMeth(project, mode="allC")
    # check all C and G found
    lapply(seq_along(genome), function(i){
        c <- ranges(matchPattern("C", genome[[i]]))
        g <- ranges(matchPattern("G", genome[[i]]))
        checkTrue(all(c == ranges(methAllC)[seqnames(methAllC) == names(genome[i]) & strand(methAllC) == "+"]))
        checkTrue(all(g == ranges(methAllC)[seqnames(methAllC) == names(genome[i]) & strand(methAllC) == "-"]))
    })
              
    methCpG <- qMeth(project, mode="CpG")
    checkTrue(all(mcols(methCpG)[,1] == mcols(methCpG)[,2] | 0 == mcols(methCpG)[,2]),
              "Test1 partial-methylated directed, mode=CpG")
    
    methCpGcomp <- qMeth(project, mode="CpGcomb")
    # check all CpG found
    lapply(seq_along(genome), function(i){
        cg <- ranges(matchPattern("CG", genome[[i]]))   
        checkTrue(all(cg == ranges(methCpGcomp)[seqnames(methCpGcomp) == names(genome[i])]))
    })
    checkTrue(all(mcols(methCpGcomp)[,1] == mcols(methCpGcomp)[,2] | 0 == mcols(methCpGcomp)[,2]),
              "Test1 partial-methylated directed, mode=CpGcomb")
    
    checkTrue(all(mcols(methAllC)[,1] == mcols(methAllC)[,2] | 0 == mcols(methAllC)[,2]),
              "Test1 partial-methylated directed, mode=allC")
    CpG <- overlapsAny(methAllC, methCpG)
    checkTrue(all(0 == mcols(methAllC)[!CpG,2]),
              "Test2 partial-methylated directed, mode=allC")
    
    unMeth <- 5000 < start(methCpG) & end(methCpG) <= 15000
    checkTrue(all(mcols(methCpG)[!unMeth,1] == mcols(methCpG)[!unMeth,2]), 
              "Test1 partial-methylated directed, mode=CpG")
    checkTrue(all(0 == mcols(methCpG)[unMeth,2]), 
              "Test2 partial-methylated directed, mode=CpG")
 
    unMeth <- 5000 < start(methCpGcomp) & end(methCpGcomp) <= 15000
    checkTrue(all(mcols(methCpGcomp)[!unMeth,1] == mcols(methCpGcomp)[!unMeth,2]), 
              "Test1 partial-methylated directed, mode=CpGcomb")
    checkTrue(all(0 == mcols(methCpGcomp)[unMeth,2]), 
              "Test2 partial-methylated directed, mode=CpGcomb")

    # check total read counts
    aln <- readGAlignmentsFromBam(project@alignments$FileName)
    cov_plus <- coverage(aln[strand(aln) == "+"])
    cov_minus <- coverage(aln[strand(aln) == "-"])
    mcols(methAllC)$cov_plus <- unlist(lapply(seq_along(cov_plus), function(i){
        as.vector(cov_plus[[i]][start(methAllC[seqnames(methAllC) == names(cov_plus)[i]])])
    }))
    mcols(methAllC)$cov_minus <- unlist(lapply(seq_along(cov_minus), function(i){
        as.vector(cov_minus[[i]][start(methAllC[seqnames(methAllC) == names(cov_minus)[i]])])
    }))
    checkTrue(all(mcols(methAllC[strand(methAllC) == "+"])[,1] == mcols(methAllC[strand(methAllC) == "+"])$cov_plus),
              "Test 1 total read counts")
    checkTrue(all(mcols(methAllC[strand(methAllC) == "-"])[,1] == mcols(methAllC[strand(methAllC) == "-"])$cov_minus),
              "Test 2 total read counts")
}
