# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

td <- tempdir()
genomeFile <- file.path("extdata", "hg19sub.fa")
sampleFile <- file.path("extdata", "samples_rna_single.txt")
auxFile <- file.path("extdata", "auxiliaries.txt")

clObj <- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
projectGenome <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
projectGenomeAux <- qAlign(sampleFile, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
stopCluster(clObj)
rm(clObj)


test_subset_project <- function()
{
    len <- length(projectGenome)
    projectGenomeS1 <- projectGenome[1:2]
    checkTrue(len/2 == length(projectGenomeS1))
    projectGenomeS2 <- projectGenome[3:4]    
    checkTrue(len/2 == length(projectGenomeS2))

    suppressWarnings(checkIdentical(projectGenome[1], projectGenome["Sample1"]))
    suppressWarnings(checkIdentical(projectGenome[3], projectGenome["Sample2"]))
    
    len <- length(projectGenomeAux)
    projectGenomeAuxS1 <- projectGenomeAux[1:2]
    checkTrue(len/2 == length(projectGenomeAuxS1))
    projectGenomeAuxS2 <- projectGenomeAux[3:4]
    checkTrue(len/2 == length(projectGenomeAuxS2))
}

test_length <- function()
{
    checkTrue(4 == length(projectGenome))
    
    checkTrue(4 == length(projectGenomeAux))
}

test_genome <- function()
{
    ## check without auxFile
    checkTrue(normalizePath(genomeFile) == normalizePath(genome(projectGenome)))

    ## check with auxFile
    checkTrue(normalizePath(genomeFile) == normalizePath(genome(projectGenomeAux)))
}

test_auxiliary <- function()
{
    ## check without auxFile
    checkTrue(0 == nrow(auxiliaries(projectGenome)))
    
    ## check with auxFile
    checkTrue(1 == nrow(auxiliaries(projectGenomeAux)))
}

test_alignment <- function()
{
    ## check without auxFile
    aln <- alignments(projectGenome)
    checkTrue(4 == nrow(aln$genome))
    checkTrue(0 == nrow(aln$aux))
    
    ## check with auxFile
    aln <- alignments(projectGenomeAux)
    checkTrue(4 == nrow(aln$genome))
    checkTrue(4 == ncol(aln$aux))
}

test_show <- function()
{
    ## check without auxFile
    show(projectGenome)
    
    ## check with auxFile
    show(projectGenomeAux)
}
