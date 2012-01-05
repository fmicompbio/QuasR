test_qCount <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    on.exit(unlink(system.file(package="QuasR", "extdata", "phage_genomes", "RbowtieIndex_phage_genomes"),
                   recursive=TRUE), add=TRUE)
    project <- qProject(sampleFile, genome=genomeName, bamfileDir=td)
    project <- qAlign(project)
    checkTrue(all(!is.na(project@env$alignments)))
    
    require(rtracklayer)
    gr <- import.gff(system.file(package="QuasR", "extdata", "NC_001422_1_genes.gff"), asRangedData=FALSE)
    checkTrue(is(gr, "GRanges"))
    ans <- qCount(project, gr, collapseSamples=TRUE)
    checkTrue(is.data.frame(ans))
    checkEquals(length(levels(project@env$samples$name)), ncol(ans))
    
    ans <- qCount(project, gr, collapseSamples=FALSE)
    checkTrue(is.data.frame(ans))
    checkEquals(nrow(project@env$samples), ncol(ans))
    
    checkEquals(as.character(1:length(gr)), rownames(ans))
    names(gr) <- paste("region", 1:length(gr), sep="")
    ans <- qCount(project, gr)
    checkEquals(names(gr), rownames(ans))
    
    seqlevels(gr) <- "TestSeq"
    checkException(qCount(project, gr))
}
