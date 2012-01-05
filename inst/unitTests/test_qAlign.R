test_qAlign_single <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    on.exit(unlink(system.file(package="QuasR", "extdata", "phage_genomes", "RbowtieIndex_phage_genomes"),
                   recursive=TRUE), add=TRUE)
    
    ## align sequence
    project <- qProject(sampleFile, genomeName, bamfileDir=td)
    checkTrue(all(is.na(project@env$alignments$genome)))
    ans <- qAlign(project)
    checkTrue(all(!is.na(ans@env$alignments$genome)))
    checkTrue(all(!is.na(project@env$alignments$genome)))
    checkEquals(ans, project) ## check if pass by reference works
    
    ## find existing bamfile
    project <- qProject(sampleFile, genomeName, bamfileDir=td)
    checkTrue(all(!is.na(project@env$alignments$genome)))
    checkEquals(ans@env$alignmentParameter, project@env$alignmentParameter)
    checkEquals(ans@env$alignments, project@env$alignments)
    ans <- qAlign(project)
    checkEquals(ans@env$alignmentParameter, project@env$alignmentParameter)
    checkEquals(ans@env$alignments, project@env$alignments)
    
    ## align to auxiliary sequence
    auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qProject(sampleFile, genomeName, auxiliaryFile=auxiliaryFile, bamfileDir=td)
    checkTrue(!is.null(project@env$auxiliaries))
    checkTrue(all(!is.na(project@env$alignments$genome)))
    checkTrue(!is.null(project@env$alignments$phiX174))
    checkTrue(all(is.na(project@env$alignments$phiX174)))
    ans <- qAlign(project)
    checkTrue(all(!is.na(ans@env$alignments$phiX174)))
    checkEquals(ans@env$alignmentParameter, project@env$alignmentParameter)
    checkEquals(ans, project)

    ## check destination of bamfile 
    checkEquals(td, unique(dirname(project@env$alignments$genome)))
    checkEquals(td, unique(dirname(project@env$alignments$phiX174)))
    ## TODO check destination of bamfileDir when NULL
    project <- qProject(sampleFile, genomeName, auxiliaryFile=auxiliaryFile)
    qAlign(project)
    checkEquals(dirname(project@env$samples$filepath), dirname(project@env$alignments$genome))
    checkEquals(dirname(project@env$samples$filepath), dirname(project@env$alignments$phiX174))
    ## TODO check paired sample
}

test_unmappedToFasta <- function()
{
    samFilename <- system.file(package="QuasR", "unitTests", "case", "test.sam")
    fastaFilename <- system.file(package="QuasR", "unitTests", "case", "test_unmapped.fa")
    bamFilename <- asBam(samFilename, tempfile(), indexDestination=FALSE, overwrite=TRUE)
    on.exit(unlink(bamFilename))
    ansFilename <- QuasR:::.unmappedToFasta(bamFilename)
    on.exit(unlink(ansFilename), add=TRUE)
    ans <- readFasta(ansFilename)
    fasta <- readFasta(fastaFilename)
    checkEquals(toString(sort(id(fasta))), toString(sort(id(ans))))
    checkEquals(toString(sort(sread(fasta))), toString(sort(sread(ans))))
    ## TODO fastaFilename <- QuasR:::.unmappedToFasta(bamFilename, tempfile())
}

test_unmappedToFastq <- function()
{
    samFilename <- system.file(package="QuasR", "unitTests", "case", "test.sam")
    fastqFilename <- system.file(package="QuasR", "unitTests", "case", "test_unmapped.fq")
    bamFilename <- asBam(samFilename, tempfile(), indexDestination=FALSE, overwrite=TRUE)
    on.exit(unlink(bamFilename))
    ansFilename <- QuasR:::.unmappedToFastq(bamFilename)
    on.exit(unlink(ansFilename), add=TRUE)
    ans <- readFastq(ansFilename)
    fastq <- readFastq(fastqFilename)
    checkEquals(toString(sort(id(fastq))), toString(sort(id(ans))))
    checkEquals(toString(sort(sread(fastq))), toString(sort(sread(ans))))
    ## TODO fastqFilename <- QuasR:::.unmappedToFastq(bamFilename, tempfile())
}