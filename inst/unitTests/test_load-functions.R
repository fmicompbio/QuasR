test_readSamples <- function()
{
    sampleFile <- system.file(package="QuasR", "extdata", "samples_dm_single.txt")
    checkTrue(file.exists(sampleFile))

    checkException(.readSamples(file=sampleFile),
                   msg="Function should not be exported.")
    
    ans <- QuasR:::.readSamples(file=sampleFile)
    checkTrue(is.data.frame(ans), 
              msg="Should be a data.frame.")
    checkEquals(c("name", "filepath", "filetype"), colnames(ans), 
                msg="Colnames should be 'name', 'filepath' and 'filetype'.")
    checkEquals(as.character(1:dim(ans)[1]), rownames(ans), 
                msg="Rownames are wrong.")
    checkTrue(!any(is.na(ans)), 
              msg="Should not contain NAs.")
    checkTrue(!any(ans == ""), 
              msg="Should not contain empty character vector.")
    ## TODO checkException() if basename(ans$filepath) is duplicated
    ## TODO checkExeption() if type is mix of fasta and fastq
    ## TODO check class(ans$name) ans$filepath ans$filetype
    sampleFile <- system.file(package="QuasR", "extdata", "samples_dm_paired.txt")
    ## TODO QuasR:::.readSamples(file=sampleFile, paired=FALSE)
    ## TODO wrong sample file format
}

test_loadAlignments <- function(){
    DEACTIVATED("Not implemented yet")
    .loadAlignments(qproject)
}

test_readAuxiliaries <- function(){
    auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    checkTrue(file.exists(auxiliaryFile))
 
    checkException(.readAuxiliaries(file=auxiliaryFile),
                   msg="Function should not be exported.")
    
    ans <- QuasR:::.readAuxiliaries(file=auxiliaryFile)
    checkTrue(is.data.frame(ans), 
              msg="Should be a data.frame.")
    checkEquals(c("feature", "filepath", "filetype"), colnames(ans), 
                msg="Colnames should be 'feature', 'filepath' and 'filetype'.")
    checkTrue(!any(is.na(ans)), 
              msg="Should not contain NAs.")
    checkTrue(!any(ans == ""), 
              msg="Should not contain empty character vector.")
}

test_loadBSgenome <- function(){
    checkTrue(require("BSgenome", quietly=TRUE), msg="Package 'BSgenome' is missing")
    bsgenomes <- installed.genomes()
    if(length(bsgenomes) == 0L)
        available.genomes()
    checkTrue(length(bsgenomes) != 0L)

    checkException(.loadBSgenome(bsgenomes[1]),
                   msg="Function should not be exported.")

    ans <- QuasR:::.loadBSgenome(bsgenomes[1])
    checkTrue(is(ans,"BSgenome"))

    ans <- QuasR:::.loadBSgenome(bsgenomes[1], lib.loc=.libPaths()[1])
    checkTrue(is(ans,"BSgenome"))

    ## test takes long
#     bsgenomes <- available.genomes()
#     genome <- bsgenomes[!(bsgenomes %in% installed.genomes())][1]
#     ans <- QuasR:::.loadBSgenome(genome)
#     checkTrue(is(ans,"BSgenome"))
#     remove.packages(genome)
}

test_loadFastaGenome <- function(){
    genomeDir <- system.file(package="QuasR", "extdata", "phage_genomes")   
    checkException(.loadFastaGenome(genomeDir),
                   msg="Function should not be exported.")
    ans <- QuasR:::.loadFastaGenome(genomeDir)
    checkEquals(genomeDir, ans$name)
    checkEquals(basename(genomeDir), ans$shortname)
    checkTrue(is.list(ans), msg="Should be a list.")
    checkTrue(!ans$bsgenome, msg="Should be FALSE")

    genomeFile <- system.file(package="QuasR", "extdata", "phage_genomes", "NC_001416.fna")
    ans <- QuasR:::.loadFastaGenome(genomeFile)
    checkEquals(genomeFile, ans$name)
    checkEquals(basename(genomeFile), ans$shortname)
    checkTrue(is.list(ans), msg="Should be a list.")
    checkTrue(!ans$bsgenome, msg="Should be FALSE")   
}

test_checkGenome <- function(){
    checkTrue(require("BSgenome", quietly=TRUE), msg="Package 'BSgenome' is missing")
    bsgenomes <- installed.genomes()
    if(length(bsgenomes) == 0L)
        available.genomes()
    checkTrue(length(bsgenomes) != 0L)

    checkException(.checkGenome(bsgenomes[1]),
                   msg="Function should not be exported.")
    ans <- QuasR:::.checkGenome(bsgenomes[1])
    checkEquals(bsgenomes[1], ans$name)
    checkTrue(ans$bsgenome, msg="Should be TRUE")
    checkTrue(is.list(ans), msg="Should be a list.")
    
    ans <- QuasR:::.checkGenome(bsgenomes[1], lib.loc=.libPaths()[1])
    checkEquals(bsgenomes[1], ans$name)
    checkTrue(ans$bsgenome, msg="Should be TRUE")
    checkTrue(is.list(ans), msg="Should be a list.")

    bsgenomes <- available.genomes()
    genome <- bsgenomes[!(bsgenomes %in% installed.genomes())][1]
    ans <- QuasR:::.checkGenome(genome)
    checkEquals(bsgenomes[1], ans$name)
    checkTrue(ans$bsgenome, msg="Should be TRUE")
    checkTrue(is.list(ans), msg="Should be a list.")

    ## TODO fasta-dir and fasta file
}

test_loadAligner <- function(){
    pkgname <- "Rbowtie"

    checkException(.loadAligner(pkgname),
                   msg="Function should not be exported.")

    ans <- QuasR:::.loadAligner(pkgname)
    checkEquals(pkgname, ans$pkgname)
    checkTrue(!is.null(ans$pkgversion))
    
    ans <- QuasR:::.loadAligner(pkgname, lib.loc=.libPaths()[1])
    checkEquals(pkgname, ans$pkgname)
    checkTrue(!is.null(ans$pkgversion))
}

test_loadIndex <- function(){
    # TODO check index path / name
    DEACTIVATED("Not implemented yet")
}