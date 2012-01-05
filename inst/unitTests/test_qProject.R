#.setUp <- function(){}
#.tearDown <- function(){}

test_initialize <- function()
{
    ans <- new("qProject") # just check this does not throw an error
    checkTrue(is(ans,"qProject"))
}
test_qProject_default <- function()
{
    checkException(qProject())
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    checkException(qProject(sampleFile))
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")

    ans <- qProject(sampleFile, genomeName)
    checkTrue(is(ans,"qProject"))

    ## check default values
    checkTrue(is.environment(ans@env))    
    checkTrue(is.data.frame(ans@env$samples))
    checkTrue(!ans@env$paired)
    checkEquals(3, length(ans@env$samples))
    checkTrue(is.null(ans@env$annotations))
    checkTrue(is.data.frame(ans@env$alignments))
    checkTrue(!is.null(ans@env$alignments$genome))
    checkTrue(is.list(ans@env$genome))
    checkEquals(genomeName, ans@env$genome$name)
    checkTrue(is.list(ans@env$aligner))
    checkEquals("Rbowtie", ans@env$aligner$pkgname)
    checkEquals(100, ans@env$maxHits)
    checkTrue(is.character(ans@env$alignmentParameter))
    checkTrue("" != ans@env$alignmentParameter)
    checkTrue(!ans@env$splicedAlignment)
    checkTrue(!ans@env$bisulfite)
    checkTrue(is.null(ans@env$bamfileDir))
    checkTrue(is.null(ans@env$indexLocation))
    checkTrue(is.null(ans@env$index))
    checkTrue(!is.null(ans@env$cacheDir))
    checkTrue(file.exists(ans@env$cacheDir))    
}

test_qProject_parameters <- function()
{
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_paired.txt")
    ans <- qProject(sampleFile, genomeName)
    checkTrue(ans@env$paired)
    checkEquals(4, length(ans@env$samples))
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")

    auxiliaryFile <- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    ans <- qProject(sampleFile, genomeName, auxiliaryFile=auxiliaryFile)
    checkTrue(is.data.frame(ans@env$auxiliaries))

    checkException(qProject(sampleFile, genomeName, aligner="RbowtieTest"))
#     ans <- qProject(sampleFile, genomeName, aligner="Rbwa")
#     checkEquals("Rbwa", ans@env$aligner$pkgname)

    ans <- qProject(sampleFile, genomeName, maxHits=47)
    checkEquals(47, ans@env$maxHits)

    ans <- qProject(sampleFile, genomeName, alignmentParameter="test")
    checkEquals("test", ans@env$alignmentParameter)

    checkException(qProject(sampleFile, genomeName, splicedAlignment=TRUE))
#     ans <- qProject(sampleFile, genomeName, splicedAlignment=TRUE)
#     checkTrue(!ans@env$splicedAlignment)

    checkException(qProject(sampleFile, genomeName, bisulfite=TRUE))
#     ans <- qProject(sampleFile, genomeName, bisulfite=TRUE)
#     checkTrue(!ans@env$bisulfite)

    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))
    ans <- qProject(sampleFile, genomeName, bamfileDir=td)
    checkEquals(td, ans@env$bamfileDir)
    ans <- qProject(sampleFile, genomeName, bamfileDir=".")
    checkEquals(tools::file_path_as_absolute("."), ans@env$bamfileDir)

    ans <- qProject(sampleFile, genomeName, indexLocation=td)
    checkEquals(td, ans@env$indexLocation)
    ans <- qProject(sampleFile, genomeName, indexLocation=".")
    checkEquals(tools::file_path_as_absolute("."), ans@env$indexLocation)
    
    ans <- qProject(sampleFile, genomeName, cacheDir=td)
    checkEquals(td, ans@env$cacheDir)
    ans <- qProject(sampleFile, genomeName, cacheDir=".")
    checkEquals(tools::file_path_as_absolute("."), ans@env$cacheDir)
}
    
test_show <- function(object){
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    project <- qProject(sampleFile, genomeName)
    checkTrue(is.null(show(project)))
}

test_qSaveProject_qReadProject <- function(project, filename){
    sampleFile <- system.file(package="QuasR", "extdata", "samples_phiX_single.txt")
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    project <- qProject(sampleFile, genomeName)

    filename1 <- qSaveProject(project)
    on.exit(unlink(filename1))
    checkTrue(file.exists(filename1))
    filename2 <- qSaveProject(project, tempfile())
    on.exit(unlink(filename2), add=TRUE)
    checkTrue(file.exists(filename2))

    ans <- qReadProject(filename2)
    checkEquals(project, ans)
}
