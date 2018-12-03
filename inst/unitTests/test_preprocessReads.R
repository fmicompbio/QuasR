# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!existsFunction("createFastaReads"))
    source(system.file(package="QuasR", "unitTests", "help_function.R"))
if(!exists("testFastaFiles"))
    testFastaFiles <- createFastaReads()
if(!exists("testFastqFiles"))
    testFastqFiles <- createFastqReads()

test_arguments <- function() {
    RUnit::checkException(preprocessReads(1L))
    RUnit::checkException(preprocessReads(c("a","b"), "c"))
    RUnit::checkException(preprocessReads("a", 1L))
    RUnit::checkException(preprocessReads("a", testFastqFiles[1]))
    RUnit::checkException(preprocessReads("a.txt", NULL))
    RUnit::checkException(preprocessReads("a.fasta", "b.fastq"))
    RUnit::checkException(preprocessReads("a.fasta.gz", "b.fasta.gz"))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa", "in2.fa", "out2.fa", Lpattern = "AAA"))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa", 1L, "out2.fa"))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa", "in2.fa", 1L))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa", "in2.fa", c("o2a.fa","o2b.fa")))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa", "in2.fa", "out2.fq"))
    RUnit::checkException(preprocessReads("in1.fa", "out1.fa.gz", "in2.fa", "out2.fa.bz2"))
}

test_paired_fasta <- function()
{
    library(ShortRead)
    faFiles <- testFastaFiles
    
    # results for default parameter
    resSoll <- matrix(c(9,NA,NA,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           complexity=0.5)
    RUnit::checkEqualsNumeric(2, res["lowComplexity",])

    # Test minLength
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           minLength=30)
    RUnit::checkEqualsNumeric(2, res["tooShort",])
    
    # Test nBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           nBases=1)
    RUnit::checkEqualsNumeric(2, res["tooManyN",])

    # Test truncateStartBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           truncateStartBases=2)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile1)
    RUnit::checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile2)
    RUnit::checkEqualsNumeric(width(readIn2)-2, width(readOut2))

    # Test truncateEndBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           truncateEndBases=3)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile1)
    RUnit::checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile2)
    RUnit::checkEqualsNumeric(width(readIn2)-3, width(readOut2))

#     # Not supported yet
#     # Test Lpattern
#     outFile1 <- tempfile(fileext=".fa")
#     outFile2 <- tempfile(fileext=".fa")
#     res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
#                            Lpattern="ATACTG")
#     checkEqualsNumeric(2, res[,1]$matchTo5pAdapter)
# 
#     # Test Rpattern
#     outFile1 <- tempfile(fileext=".fa")
#     outFile2 <- tempfile(fileext=".fa")
#     res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
#                            Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
#     checkEqualsNumeric(2, res[,1]$matchTo3pAdapter)
}

test_single_fasta <- function()
{
    library(ShortRead)
    faFiles <- testFastaFiles
    
    # results for default parameter
    resSoll <- matrix(0, nrow=7, ncol=2)
    resSoll[,1] <- c(9,0,0,0,0,0,9)
    resSoll[,2] <- c(9,0,0,0,0,0,9)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")

    # Test Default
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, complexity=0.5)
    RUnit::checkEqualsNumeric(c(1,1), res["lowComplexity",])
    
    # Test minLength
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, minLength=30)
    RUnit::checkEqualsNumeric(c(1,1), res["tooShort",])
    
    # Test nBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, nBases=1)
    RUnit::checkEqualsNumeric(c(2,1), res["tooManyN",])

    # Test truncateStartBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, truncateStartBases=2)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile[1])
    RUnit::checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile[2])
    RUnit::checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    # TODO check seq

    # Test truncateEndBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, truncateEndBases=3)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile[1])
    RUnit::checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile[2])
    RUnit::checkEqualsNumeric(width(readIn2)-3, width(readOut2))
    # TODO check seq

    # Test Lpattern
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, Lpattern="ATACTG")
    RUnit::checkEqualsNumeric(c(7,4), res["matchTo5pAdapter",])

    # Test Rpattern
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
    RUnit::checkEqualsNumeric(c(6,5), res["matchTo3pAdapter",])
        
#     # Test combination of parameters
#     outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
#     res <- preprocessReads(faFiles, outFile, 
#                            Lpattern="ATACTG", 
#                            Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")    
#     
#     outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
#     res <- preprocessReads(faFiles, outFile, 
#                            truncateEndBases=4, 
#                            Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
#     
#     outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
#     res <- preprocessReads(faFiles, outFile, 
#                            complexity=0.5,
#                            truncateEndBases=4, 
#                            Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
}

test_paired_fastq <- function()
{
    library(ShortRead)
    fqFiles <- testFastqFiles
    
    # results for default parameter
    resSoll <- matrix(c(9,NA,NA,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           complexity=0.5)
    RUnit::checkEqualsNumeric(2, res["lowComplexity",])
    
    # Test minLength
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           minLength=30)
    RUnit::checkEqualsNumeric(2, res["tooShort",])
    
    # Test nBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           nBases=1)
    RUnit::checkEqualsNumeric(2, res["tooManyN",])
    
    # Test truncateStartBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           truncateStartBases=2)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile1)
    RUnit::checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile2)
    RUnit::checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    
    # Test truncateEndBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           truncateEndBases=3)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile1)
    RUnit::checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile2)
    RUnit::checkEqualsNumeric(width(readIn2)-3, width(readOut2))
}

test_single_fastq <- function()
{
    library(ShortRead)
    fqFiles <- testFastqFiles
    
    # results for default parameter
    resSoll <- matrix(0, nrow=7, ncol=2)
    resSoll[,1] <- c(9,0,0,0,0,0,9)
    resSoll[,2] <- c(9,0,0,0,0,0,9)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, complexity=0.5)
    RUnit::checkEqualsNumeric(c(1,1), res["lowComplexity",])
    
    # Test minLength
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, minLength=30)
    RUnit::checkEqualsNumeric(c(1,1), res["tooShort",])
    
    # Test nBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, nBases=1)
    RUnit::checkEqualsNumeric(c(2,1), res["tooManyN",])
    
    # Test truncateStartBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, truncateStartBases=2)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile[1])
    RUnit::checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile[2])
    RUnit::checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    # TODO check seq
    
    # Test truncateEndBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, truncateEndBases=3)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile[1])
    RUnit::checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile[2])
    RUnit::checkEqualsNumeric(width(readIn2)-3, width(readOut2))
    # TODO check seq
    
    # Test Lpattern
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, Lpattern="ATACTG")
    RUnit::checkEqualsNumeric(c(7,4), res["matchTo5pAdapter",])
    
    # Test Rpattern
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
    RUnit::checkEqualsNumeric(c(6,5), res["matchTo3pAdapter",])
}

test_input_fasta <- function()
{
    input <- testFastaFiles[1]
    
    # results for default parameter
    resSoll <- matrix(c(9,0,0,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")

    outFile <- tempfile(fileext=".fa.gz")
    res <- preprocessReads(input, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    outFile <- tempfile(fileext=".fa.bz2")
    res <- preprocessReads(input, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    outFile <- tempfile(fileext=".fa.xz")
    res <- preprocessReads(input, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
}

test_input_fastq <- function()
{
    str <- scan(testFastqFiles[1], character())
    
    # results for default parameter
    resSoll <- matrix(c(9,0,0,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    gzFile <- tempfile(fileext=".fq.gz")
    outFile <- tempfile(fileext=".fq.gz")
    con <- gzfile(gzFile, "w")
    cat(str, sep="\n", file=con)
    close(con)
    res <- preprocessReads(gzFile, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    bzFile <- tempfile(fileext=".fq.bz2")
    outFile <- tempfile(fileext=".fq.bz2")
    con <- bzfile(bzFile, "w")
    cat(str, sep="\n", file=con)
    close(con)
    res <- preprocessReads(bzFile, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
    
    xzFile <- tempfile(fileext=".fq.xz")
    outFile <- tempfile(fileext=".fq.xz")
    con <- xzfile(xzFile, "w")
    cat(str, sep="\n", file=con)
    close(con)
    res <- preprocessReads(xzFile, outFile)
    RUnit::checkEqualsNumeric(resSoll, res)
}

