test_QuasR_BSgenome <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_dm_single.txt")
    auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    genomeName <- "BSgenome.Dmelanogaster.UCSC.dm3"
    project <- qProject(sampleFile, genome=genomeName, auxiliaryFile=auxiliaryFile, bamfileDir=td, aligner="Rbowtie")
    project <- qAlign(project)
    seqInfo <- seqlengths(Dmelanogaster)
    gr <- GRanges(seqnames=names(seqInfo), ranges=IRanges(1, width=seqInfo), seqlengths=seqInfo)
    cnt <- qCount(project, gr)
    cnt## align with existing index
    project <- qAlign(project)
}


test_QuasR_FastaDirGenome <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    project <- qProject(sampleFile, genome=genomeName, bamfileDir=td, aligner="Rbowtie")
    project <- qAlign(project)
    require(rtracklayer)
    gr <- import.gff(system.file(package="QuasR", "extdata", "NC_001422_1_genes.gff"), asRangedData=FALSE)
    cnt <- qCount(project, gr)
    ## align with existing index
    project <- qAlign(project)
}

test_QuasR_FastaFileGenome <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes", "NC_001422.1.fa")
    project <- qProject(sampleFile, genome=genomeName, bamfileDir=td, aligner="Rbowtie")
    project <- qAlign(project)
    require(rtracklayer)
    gr <- import.gff(system.file(package="QuasR", "extdata", "NC_001422_1_genes.gff"), asRangedData=FALSE)
    cnt <- qCount(project, gr)
    ## align with existing index
    project <- qAlign(project)
}
