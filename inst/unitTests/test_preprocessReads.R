# Fragment design
# Left Adapter  ATACTG
# Insert        TGTGACAGACCTCGGGGCCACATGCACTGACTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC
# Right Adapter ATCTCGTATGCCGTCTTCTGCTTG
createFastaReads <- function(){
    faFileName1 <- tempfile(fileext=".fa")
    faFile1 <-file(faFileName1, open="w")
    writeLines(">seq1", con=faFile1)
    writeLines("CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", con=faFile1) # onlyInsert
    writeLines(">seq2", con=faFile1)
    writeLines("AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", con=faFile1) # low compexity, 2N
    writeLines(">seq3", con=faFile1)
    writeLines("AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", con=faFile1) # A.., LAdapter
    writeLines(">seq4", con=faFile1)
    writeLines("ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", con=faFile1) # LAdapter
    writeLines(">seq5", con=faFile1)
    writeLines("TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", con=faFile1) # Insert
    writeLines(">seq6", con=faFile1)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=faFile1) # LAdapter, noInsert, RAdapter, A..
    writeLines(">seq7", con=faFile1)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", con=faFile1) # LAdapter, noInsert, RAdapter, A.., 2N  
    writeLines(">seq8", con=faFile1)
    writeLines("ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", con=faFile1) # LAdapter, shortInsert, RAdapter
    writeLines(">seq9", con=faFile1)
    writeLines("TGACAGACCTCGGGGCCAC", con=faFile1)               # shortLength
    close(faFile1)
    
    faFileName2 <- tempfile(fileext=".fa")
    faFile2 <-file(faFileName2, open="w")
    writeLines(">seq1", con=faFile2)
    writeLines("AAAAAAAAAAAAAAAAAAA", con=faFile2) # low compexity
    writeLines(">seq2", con=faFile2)
    writeLines("TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", con=faFile2) # low compexity, shortLength
    writeLines(">seq3", con=faFile2)
    writeLines("CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", con=faFile2) # onlyInsert
    writeLines(">seq4", con=faFile2)
    writeLines("GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=faFile2) # RAdapter
    writeLines(">seq5", con=faFile2)
    writeLines("AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", con=faFile2) # Insert
    writeLines(">seq6", con=faFile2)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=faFile2) # LAdapter, noInsert, RAdapter, A..
    writeLines(">seq7", con=faFile2)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", con=faFile2) # LAdapter, noInsert, RAdapter, A.., 2N
    writeLines(">seq8", con=faFile2)
    writeLines("AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", con=faFile2) # A.., LAdapter, shortInsert, RAdapter
    writeLines(">seq9", con=faFile2)
    writeLines("CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT", con=faFile2)   # Insert
    close(faFile1)
    return(c(faFileName1, faFileName2))
}

createFastqReads <- function(){
    fqFileName1 <- tempfile(fileext=".fq")
    fqFile1 <- file(fqFileName1, open="w")
    writeLines("@seq1", con=fqFile1)
    writeLines("CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("5..49<494*<49493####################", con=fqFile1)
    writeLines("@seq2", con=fqFile1)
    writeLines("AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("+4544944444444######################", con=fqFile1)
    writeLines("@seq3", con=fqFile1)
    writeLines("AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", con=fqFile1)
    writeLines("+", con=fqFile1)    
    writeLines("BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@/=", con=fqFile1)
    writeLines("@seq4", con=fqFile1)
    writeLines("ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABC", con=fqFile1)
    writeLines("@seq5", con=fqFile1)
    writeLines("TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("BA?AB??60>A6?0BBBB95>057;@*.<(434;+4", con=fqFile1)
    writeLines("@seq6", con=fqFile1)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>BBB@", con=fqFile1)
    writeLines("@seq7", con=fqFile1)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("####################################", con=fqFile1)    
    writeLines("@seq8", con=fqFile1)
    writeLines("ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("(33(;?B@AB43*,/;9(6</7>5;<##########", con=fqFile1)    
    writeLines("@seq9", con=fqFile1)
    writeLines("TGACAGACCTCGGGGCCAC", con=fqFile1)
    writeLines("+", con=fqFile1)
    writeLines("###################", con=fqFile1)
    close(fqFile1)
    
    fqFileName2 <- tempfile(fileext=".fq")
    fqFile2 <- file(fqFileName2, open="w")    
    writeLines("@seq1", con=fqFile2)
    writeLines("AAAAAAAAAAAAAAAAAAA", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("3@?3/>9A8@A1-*/4@BB", con=fqFile2)
    writeLines("@seq2", con=fqFile2)
    writeLines("TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("4-093<##############################", con=fqFile2)
    writeLines("@seq3", con=fqFile2)
    writeLines("CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("CB?(8=(<A/<=-(07+7&883@#############", con=fqFile2)
    writeLines("@seq4", con=fqFile2)
    writeLines("GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("5..49<494*<4949#####################", con=fqFile2)
    writeLines("@seq5", con=fqFile2)
    writeLines("AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("+45449444444########################", con=fqFile2)
    writeLines("@seq6", con=fqFile2)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@A<", con=fqFile2)
    writeLines("@seq7", con=fqFile2)
    writeLines("ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABB", con=fqFile2)
    writeLines("@seq8", con=fqFile2)
    writeLines("AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("BA?AB??60>A6?0BBBB95>057;@*.<(434;@B", con=fqFile2)
    writeLines("@seq9", con=fqFile2)
    writeLines("CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT", con=fqFile2)
    writeLines("+", con=fqFile2)
    writeLines("+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>B@4;", con=fqFile2)
    close(fqFile2)
    return(c(fqFileName1, fqFileName2))
}

test_paired_fasta <- function()
{
    faFiles <- createFastaReads()
    
    # results for default parameter
    resSoll <- matrix(c(9,NA,NA,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2)
    checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           complexity=0.5)
    checkEqualsNumeric(2, res["lowComplexity",])

    # Test minLength
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           minLength=30)
    checkEqualsNumeric(2, res["tooShort",])
    
    # Test nBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           nBases=1)
    checkEqualsNumeric(2, res["tooManyN",])

    # Test truncateStartBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           truncateStartBases=2)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile1)
    checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile2)
    checkEqualsNumeric(width(readIn2)-2, width(readOut2))

    # Test truncateEndBases
    outFile1 <- tempfile(fileext=".fa")
    outFile2 <- tempfile(fileext=".fa")
    res <- preprocessReads(faFiles[1], outFile1, faFiles[2], outFile2,
                           truncateEndBases=3)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile1)
    checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile2)
    checkEqualsNumeric(width(readIn2)-3, width(readOut2))

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
    faFiles <- createFastaReads()
    
    # results for default parameter
    resSoll <- matrix(0, nrow=7, ncol=2)
    resSoll[,1] <- c(9,0,0,0,0,0,9)
    resSoll[,2] <- c(9,0,0,0,0,0,9)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")

    # Test Default
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile)
    checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, complexity=0.5)
    checkEqualsNumeric(c(1,1), res["lowComplexity",])
    
    # Test minLength
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, minLength=30)
    checkEqualsNumeric(c(1,1), res["tooShort",])
    
    # Test nBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, nBases=1)
    checkEqualsNumeric(c(2,1), res["tooManyN",])

    # Test truncateStartBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, truncateStartBases=2)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile[1])
    checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile[2])
    checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    # TODO check seq

    # Test truncateEndBases
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, truncateEndBases=3)
    readIn1 <- readFasta(faFiles[1])
    readOut1 <- readFasta(outFile[1])
    checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFasta(faFiles[2])
    readOut2 <- readFasta(outFile[2])
    checkEqualsNumeric(width(readIn2)-3, width(readOut2))
    # TODO check seq

    # Test Lpattern
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, Lpattern="ATACTG")
    checkEqualsNumeric(c(7,4), res["matchTo5pAdapter",])

    # Test Rpattern
    outFile <- c(tempfile(fileext=".fa"), tempfile(fileext=".fa"))
    res <- preprocessReads(faFiles, outFile, Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
    checkEqualsNumeric(c(6,5), res["matchTo3pAdapter",])
        
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
    fqFiles <- createFastqReads()
    
    # results for default parameter
    resSoll <- matrix(c(9,NA,NA,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2)
    checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           complexity=0.5)
    checkEqualsNumeric(2, res["lowComplexity",])
    
    # Test minLength
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           minLength=30)
    checkEqualsNumeric(2, res["tooShort",])
    
    # Test nBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           nBases=1)
    checkEqualsNumeric(2, res["tooManyN",])
    
    # Test truncateStartBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           truncateStartBases=2)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile1)
    checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile2)
    checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    
    # Test truncateEndBases
    outFile1 <- tempfile(fileext=".fq")
    outFile2 <- tempfile(fileext=".fq")
    res <- preprocessReads(fqFiles[1], outFile1, fqFiles[2], outFile2,
                           truncateEndBases=3)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile1)
    checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile2)
    checkEqualsNumeric(width(readIn2)-3, width(readOut2))
}

test_single_fastq <- function()
{
    fqFiles <- createFastqReads()
    
    # results for default parameter
    resSoll <- matrix(0, nrow=7, ncol=2)
    resSoll[,1] <- c(9,0,0,0,0,0,9)
    resSoll[,2] <- c(9,0,0,0,0,0,9)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")
    
    # Test Default
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile)
    checkEqualsNumeric(resSoll, res)
    
    # Test complexity
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, complexity=0.5)
    checkEqualsNumeric(c(1,1), res["lowComplexity",])
    
    # Test minLength
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, minLength=30)
    checkEqualsNumeric(c(1,1), res["tooShort",])
    
    # Test nBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, nBases=1)
    checkEqualsNumeric(c(2,1), res["tooManyN",])
    
    # Test truncateStartBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, truncateStartBases=2)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile[1])
    checkEqualsNumeric(width(readIn1)-2, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile[2])
    checkEqualsNumeric(width(readIn2)-2, width(readOut2))
    # TODO check seq
    
    # Test truncateEndBases
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, truncateEndBases=3)
    readIn1 <- readFastq(fqFiles[1])
    readOut1 <- readFastq(outFile[1])
    checkEqualsNumeric(width(readIn1)-3, width(readOut1))
    readIn2 <- readFastq(fqFiles[2])
    readOut2 <- readFastq(outFile[2])
    checkEqualsNumeric(width(readIn2)-3, width(readOut2))
    # TODO check seq
    
    # Test Lpattern
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, Lpattern="ATACTG")
    checkEqualsNumeric(c(7,4), res["matchTo5pAdapter",])
    
    # Test Rpattern
    outFile <- c(tempfile(fileext=".fq"), tempfile(fileext=".fq"))
    res <- preprocessReads(fqFiles, outFile, Rpattern="ATCTCGTATGCCGTCTTCTGCTTG")
    checkEqualsNumeric(c(6,5), res["matchTo3pAdapter",])
}

test_input_fasta <- function()
{
    input <- createFastaReads()[1]
    
    # results for default parameter
    resSoll <- matrix(c(9,0,0,0,0,0,9), nrow=7, ncol=1)
    rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                           "tooShort","tooManyN","lowComplexity","totalPassed")

    outFile <- tempfile(fileext=".fa.gz")
    res <- preprocessReads(input, outFile)
    checkEqualsNumeric(resSoll, res)
    
    outFile <- tempfile(fileext=".fa.bz2")
    res <- preprocessReads(input, outFile)
    checkEqualsNumeric(resSoll, res)
    
    outFile <- tempfile(fileext=".fa.xz")
    res <- preprocessReads(input, outFile)
    checkEqualsNumeric(resSoll, res)
}

test_input_fastq <- function()
{
    str <- scan(createFastqReads()[1], character())
    
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
    checkEqualsNumeric(resSoll, res)
    
    bzFile <- tempfile(fileext=".fq.bz2")
    outFile <- tempfile(fileext=".fq.bz2")
    con <- bzfile(bzFile, "w")
    cat(str, sep="\n", file=con)
    close(con)
    res <- preprocessReads(bzFile, outFile)
    checkEqualsNumeric(resSoll, res)
    
    xzFile <- tempfile(fileext=".fq.xz")
    outFile <- tempfile(fileext=".fq.xz")
    con <- xzfile(xzFile, "w")
    cat(str, sep="\n", file=con)
    close(con)
    res <- preprocessReads(xzFile, outFile)
    checkEqualsNumeric(resSoll, res)
}