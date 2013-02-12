test_subset_project <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")

    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    len <- length(project)
    projectS1 <- project[1:2]
    checkTrue(len/2 == length(projectS1))
    projectS2 <- project[3:4]    
    checkTrue(len/2 == length(projectS2))

    suppressWarnings(checkIdentical(project[1], project["Sample1"]))
    suppressWarnings(checkIdentical(project[3], project["Sample2"]))
    
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    len <- length(project)
    projectS1 <- project[1:2]
    checkTrue(len/2 == length(projectS1))
    projectS2 <- project[3:4]
    checkTrue(len/2 == length(projectS2))
}

test_length <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    checkTrue(4 == length(project))
    
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    checkTrue(4 == length(project))
}

test_genome <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    
    ## check without auxFile
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    checkTrue(normalizePath(genomeFile) == normalizePath(genome(project)))

    ## check with auxFile
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    checkTrue(normalizePath(genomeFile) == normalizePath(genome(project)))
}

test_auxiliary <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    
    ## check without auxFile
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    checkTrue(0 == nrow(auxiliaries(project)))
    
    ## check with auxFile
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    checkTrue(1 == nrow(auxiliaries(project)))
}

test_alignment <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    
    ## check without auxFile
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    aln <- alignments(project)
    checkTrue(4 == nrow(aln$genome))
    checkTrue(0 == nrow(aln$aux))
    
    ## check with auxFile
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    aln <- alignments(project)
    checkTrue(4 == nrow(aln$genome))
    checkTrue(4 == ncol(aln$aux))
}

test_show <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    td <- tempdir()
    genomeFile <- file.path("extdata", "hg19sub.fa")
    sampleFile <- file.path("extdata", "samples_rna_single.txt")
    auxFile <- file.path("extdata", "auxiliaries.txt")
    
    ## check without auxFile
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    show(project)
    
    ## check with auxFile
    project <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    show(project)
}
