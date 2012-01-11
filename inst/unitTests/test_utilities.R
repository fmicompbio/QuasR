# test.requirePkg <- function()
# {
#   DEACTIVATED("Not implemented yet")
# }

test_baseFileName <- function()
{
  filenames <- character()
  res <- QuasR:::.baseFileName(filenames)
  checkIdentical(NULL, res)
    
  filenames <- "test.txt"
  res <- QuasR:::.baseFileName(filenames)
  checkIdentical("test", res)
  
  filenames <- c("/dirA/dirB/file.txt",
                 "dirB/file.txt",
                 "/dirA/dirB/file.test.txt",
                 "dirB/file.test.txt",
                 "file.test.txt",
                 "/dirA/dirB/file",
                 "dirB/file",
                 "file",
                 "")
  targetNames <- c("file","file","file.test","file.test",
                   "file.test","file","file","file","")  
  res <- QuasR:::.baseFileName(filenames)
  checkIdentical(targetNames, res)
}

test_fileExtension <- function()
{
  filenames <- character()
  res <- QuasR:::.fileExtension(filenames)
  checkIdentical(character(), res)
    
  filenames <- "test.txt"
  res <- QuasR:::.fileExtension(filenames)
  checkIdentical("txt", res)
  
  filenames <- c("/dirA/dirB/file.txt",
                 "dirB/file.txt",
                 "/dirA/dirB/file.test.txt",
                 "dirB/file.test.txt",
                 "file.test.txt",
                 "/dirA/dirB/file",
                 "dirB/file",
                 "file",
                 "")
  targetExt <- c("txt","txt","txt","txt","txt","","","","")  
  res <- QuasR:::.fileExtension(filenames)
  checkIdentical(targetExt, res)  
}

# test.progressReport <- function()
# {
#   DEACTIVATED("Not implemented yet")
# }

# test.multiToSingleFasta <- function()
# {
#   DEACTIVATED("Not implemented yet")
# }

test_getBamFile <- function()
{
    td <- tempfile()
    checkTrue(dir.create(td, showWarnings=FALSE, recursive=TRUE))
    on.exit(unlink(td, recursive=TRUE))

    sampleBaseName <- "test"
    
    ## default 
    genomeName <- system.file(package="QuasR", "extdata", "phage_genomes")
    alignerParameter <- "bowtie -v 2 -k 99 -m 99 --best --strata -q"
    hf <- tempfile(tmpdir=td, fileext=".sam")
    cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
        sprintf("ap:%s", alignerParameter),
        sprintf("at:%s", genomeName),
        "CL:\"Add or modifiy IH Tag\"\nsq1\t4\t*\t0\t0\t*\t*\t0\t0\tAAA\t###\tXM:i:0\tIH:i:0", sep="\t", file=hf)
    bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
    listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
    ans <- QuasR:::.getBamFile(listBamFilenames, genomeName, alignerParameter)
    checkEquals(bamfile, ans)
    
    ## change of alignment Parameter
    ## find the right file out of several files
    alignerParameter <- "bowtie -v 2 -k 5 -m 5 --best --strata -q"
    hf <- tempfile(tmpdir=td, fileext=".sam")
    cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
        sprintf("ap:%s", alignerParameter),
        sprintf("at:%s", genomeName),
        "CL:\"Add or modifiy IH Tag\"", sep="\t", file=hf)
    bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
    listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
    ans <- QuasR:::.getBamFile(listBamFilenames, genomeName, alignerParameter)
    checkEquals(bamfile, ans)

    ## change the genome
    ## find the right file out of several files
    genomeName <- "BSgenome"
    hf <- tempfile(tmpdir=td, fileext=".sam")
    cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
        sprintf("ap:%s", alignerParameter),
        sprintf("at:%s", genomeName),
        "CL:\"Add or modifiy IH Tag\"", sep="\t", file=hf)
    bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
    listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
    ans <- QuasR:::.getBamFile(listBamFilenames, genomeName, alignerParameter)
    checkEquals(bamfile, ans)

    ## create a second identical bamfile and check if exception occur
    ## uncertain situation the process should not continue
    hf <- tempfile(tmpdir=td, fileext=".sam")
    cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr2L\tLN:23011544\n@PG\tID:QuasR\tPN:QuasR\tVN:0.1.1",
        sprintf("ap:%s", alignerParameter),
        sprintf("at:%s", genomeName),
        "CL:\"Add or modifiy IH Tag\"", sep="\t", file=hf)
    bamfile <- asBam(hf, tempfile(pattern=sampleBaseName, tmpdir=td), indexDestination=F)
    listBamFilenames <- list.files(path=td, pattern=sprintf("%s.*\\.bam$", sampleBaseName), full.names=TRUE)
    checkException(QuasR:::.getBamFile(listBamFilenames, genomeName, alignerParameter))
    file.remove(bamfile)

    ## check if exception occurs if parameters are NULL
    checkException(QuasR:::.getBamFile(listBamFilenames, genomeName, NULL))
    checkException(QuasR:::.getBamFile(listBamFilenames, NULL, alignerParameter)) 
}
