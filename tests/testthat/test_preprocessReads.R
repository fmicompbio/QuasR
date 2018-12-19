context("preprocessReads")

# prepare sample data
#   Fragment design:
#     Left Adapter  ATACTG
#     Insert        TGTGACAGACCTCGGGGCCACATGCACTGACTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC
#     Right Adapter ATCTCGTATGCCGTCTTCTGCTTG
# ... create two temporary fasta files and return file names
faFiles <- tempfile(fileext = rep(".fa", 2), tmpdir = "extdata")
writeLines(c(">seq1", "CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", # onlyInsert
             ">seq2", "AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", # low compexity, 2N
             ">seq3", "AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", # A.., LAdapter
             ">seq4", "ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", # LAdapter
             ">seq5", "TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", # Insert
             ">seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # LAdapter, noInsert, RAdapter, A..
             ">seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", # LAdapter, noInsert, RAdapter, A.., 2N
             ">seq8", "ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", # LAdapter, shortInsert, RAdapter
             ">seq9", "TGACAGACCTCGGGGCCAC"),                 # shortLength
             faFiles[1])
writeLines(c(">seq1", "AAAAAAAAAAAAAAAAAAA",                  # low compexity
             ">seq2", "TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", # low compexity, shortLength
             ">seq3", "CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", # onlyInsert
             ">seq4", "GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # RAdapter
             ">seq5", "AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", # Insert
             ">seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # LAdapter, noInsert, RAdapter, A..
             ">seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", # LAdapter, noInsert, RAdapter, A.., 2N
             ">seq8", "AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", # A.., LAdapter, shortInsert, RAdapter
             ">seq9", "CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT"),# Insert
             faFiles[2])
# ... create two temporary fastq files and return file names
fqFiles <- tempfile(fileext = rep(".fq", 2), tmpdir = "extdata")
writeLines(c("@seq1", "CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", "+", "5..49<494*<49493####################",
             "@seq2", "AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "+", "+4544944444444######################",
             "@seq3", "AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", "+", "BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@/=",
             "@seq4", "ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", "+", "BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABC",
             "@seq5", "TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", "+", "BA?AB??60>A6?0BBBB95>057;@*.<(434;+4",
             "@seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>BBB@",
             "@seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", "+", "####################################",
             "@seq8", "ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", "+", "(33(;?B@AB43*,/;9(6</7>5;<##########",
             "@seq9", "TGACAGACCTCGGGGCCAC",                  "+", "###################"),
             fqFiles[1])
writeLines(c("@seq1", "AAAAAAAAAAAAAAAAAAA",                  "+", "3@?3/>9A8@A1-*/4@BB",
             "@seq2", "TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", "+", "4-093<##############################",
             "@seq3", "CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", "+", "CB?(8=(<A/<=-(07+7&883@#############",
             "@seq4", "GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "5..49<494*<4949#####################",
             "@seq5", "AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", "+", "+45449444444########################",
             "@seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@A<",
             "@seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", "+", "BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABB",
             "@seq8", "AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", "+", "BA?AB??60>A6?0BBBB95>057;@*.<(434;@B",
             "@seq9", "CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT", "+", "+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>B@4;"),
             fqFiles[2])

test_that("arguments in preprocessReads are digested correctly", {
  expect_error(preprocessReads(1L),
               ".filename. must be of type character")
  expect_error(preprocessReads(c("a","b"), "c"),
               ".filename. and .outputFilename. must have equal length")
  expect_error(preprocessReads("a", 1L),
               ".outputFilename. must be of type character")
  expect_error(preprocessReads("a.txt", NULL),
               "unsupported file format")
  expect_error(preprocessReads("a.fasta", "b.fastq"),
               "format of .filename. and .outputFilename. must be identical")
  expect_error(preprocessReads("a.fa.gz", "b.fa.gz"),
               "compressed .fasta. input is not yet supported")
  expect_error(preprocessReads("in1.fa", "out1.fa", "in2.fa", "out2.fa", Lpattern = "AAA"),
               "Removing adapters from paired-end samples is not yet supported")
  expect_error(preprocessReads("in1.fa", "out1.fa", 1L, "out2.fa"),
               ".filenameMate. must be of type character")
  expect_error(preprocessReads("in1.fa", "out1.fa", "in2.fa", 1L),
               ".outputFilenameMate. must be of type character")
  expect_error(preprocessReads("in1.fa", "out1.fa", "in2.fa", c("o2a.fa","o2b.fa")),
               ".filenameMate. and .outputFilenameMate. must have equal length")
  expect_error(preprocessReads("in1.fa", "out1.fa", "in2.fa", "out2.fq"),
               "format of .filenameMate. and .outputFilenameMate. must be identical")
  expect_error(preprocessReads("in1.fa", "out1.fa.gz", "in2.fa", "out2.fa.bz2"),
               "compression format of .outputFilename. and .outputFilenameMate. must be identical")
})

test_that("preprocessReads correctly processes single-read files", {
  # results for default parameters
  outFilesFasta <- tempfile(fileext = rep(".fa", 2), tmpdir = "extdata")
  resSoll <- cbind(c(9,0,0,0,0,0,9), c(9,0,0,0,0,0,9))
  rownames(resSoll) <- c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                         "tooShort","tooManyN","lowComplexity","totalPassed")
  colnames(resSoll) <- basename(faFiles)
  outFilesFastq <- tempfile(fileext = rep(".fq", 2), tmpdir = "extdata")

  # default
  res <- preprocessReads(faFiles, outFilesFasta)
  expect_equal(resSoll, res)
  unlink(outFilesFasta)
  
  # complexity
  res <- preprocessReads(fqFiles, outFilesFastq, complexity = 0.5)
  expect_equal(c(1, 1), unname(res["lowComplexity",]))
  unlink(outFilesFastq)
  
  # minLength
  res <- preprocessReads(faFiles, outFilesFasta, minLength = 30)
  expect_equal(c(1, 1), unname(res["tooShort",]))
  unlink(outFilesFasta)
  
  # nBases
  res <- preprocessReads(fqFiles, outFilesFastq, nBases = 1)
  expect_equal(c(2, 1), unname(res["tooManyN",]))
  unlink(outFilesFastq)
  
  # truncateStartBases
  res <- preprocessReads(faFiles, outFilesFasta, truncateStartBases = 2)
  expect_equal(nchar(readLines(faFiles[1]))[seq(2, 18, by = 2)] - 2L,
               nchar(readLines(outFilesFasta[1])[seq(2, 18, by = 2)]))
  unlink(outFilesFasta)
  
  # truncateEndBases
  res <- preprocessReads(fqFiles, outFilesFastq, truncateEndBases = 3)
  expect_equal(nchar(readLines(fqFiles[2]))[seq(2, 36, by = 4)] - 3L,
               nchar(readLines(outFilesFastq[2])[seq(2, 36, by = 4)]))
  unlink(outFilesFastq)
  
  # Lpattern
  res <- preprocessReads(faFiles, outFilesFasta, Lpattern = "ATACTG")
  expect_equal(c(7, 4), unname(res["matchTo5pAdapter",]))
  unlink(outFilesFasta)
  
  # Rpattern
  res <- preprocessReads(fqFiles, outFilesFastq, Rpattern = "ATCTCGTATGCCGTCTTCTGCTTG")
  expect_equal(c(6, 5), unname(res["matchTo3pAdapter",]))
  unlink(outFilesFastq)
})

test_that("preprocessReads correctly processes paired-end files", {
  # results for default parameters
  outFilesFasta <- tempfile(fileext = rep(".fa", 2), tmpdir = "extdata")
  resSoll <- matrix(c(9,NA,NA,0,0,0,9), nrow = 7, ncol = 1,
                    dimnames = list(c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                                      "tooShort","tooManyN","lowComplexity","totalPassed"),
                                    paste(basename(faFiles), collapse = ":")))
  outFilesFastq <- tempfile(fileext = rep(".fq", 2), tmpdir = "extdata")
  
  # default
  res <- preprocessReads(faFiles[1], outFilesFasta[1], faFiles[2], outFilesFasta[2])
  expect_equal(resSoll, res)
  unlink(outFilesFasta)
  
  # complexity
  res <- preprocessReads(fqFiles[1], outFilesFastq[1], fqFiles[2], outFilesFastq[2],
                         complexity = 0.5)
  expect_equal(2, res["lowComplexity",])
  unlink(outFilesFastq)
  
  # minLength
  res <- preprocessReads(faFiles[1], outFilesFasta[1], faFiles[2], outFilesFasta[2],
                         minLength = 30)
  expect_equal(2, res["tooShort",])
  unlink(outFilesFasta)
  
  # nBases
  res <- preprocessReads(fqFiles[1], outFilesFastq[1], fqFiles[2], outFilesFastq[2],
                         nBases = 1)
  expect_equal(2, res["tooManyN",])
  unlink(outFilesFastq)
  
  # truncateStartBases
  res <- preprocessReads(faFiles[1], outFilesFasta[1], faFiles[2], outFilesFasta[2],
                         truncateStartBases = 2)
  expect_equal(nchar(readLines(faFiles[1]))[seq(2, 18, by = 2)] - 2L,
               nchar(readLines(outFilesFasta[1])[seq(2, 18, by = 2)]))
  expect_equal(nchar(readLines(faFiles[2]))[seq(2, 18, by = 2)] - 2L,
               nchar(readLines(outFilesFasta[2])[seq(2, 18, by = 2)]))
  unlink(outFilesFasta)
  
  # truncateEndBases
  res <- preprocessReads(fqFiles[1], outFilesFastq[1], fqFiles[2], outFilesFastq[2],
                         truncateEndBases = 3)
  expect_equal(nchar(readLines(fqFiles[1]))[seq(2, 36, by = 4)] - 3L,
               nchar(readLines(outFilesFastq[1])[seq(2, 36, by = 4)]))
  expect_equal(nchar(readLines(fqFiles[2]))[seq(2, 36, by = 4)] - 3L,
               nchar(readLines(outFilesFastq[2])[seq(2, 36, by = 4)]))
  unlink(outFilesFastq)
})

test_that("preprocessReads correctly (de-)compresses files", {
  # results for default parameters
  resSoll <- matrix(c(9,0,0,0,0,0,9), ncol = 1,
                    dimnames = list(c("totalSequences","matchTo5pAdapter","matchTo3pAdapter",
                                      "tooShort","tooManyN","lowComplexity","totalPassed"),
                                    basename(faFiles[1])))
  
  outFile <- tempfile(fileext = ".fa.gz", tmpdir = "extdata")
  res <- preprocessReads(faFiles[1], outFile)
  expect_equal(resSoll, res)
  unlink(outFile)
  
  outFile <- tempfile(fileext = ".fa.bz2", tmpdir = "extdata")
  res <- preprocessReads(faFiles[1], outFile)
  expect_equal(resSoll, res)
  unlink(outFile)
  
  outFile <- tempfile(fileext = ".fa.xz", tmpdir = "extdata")
  res <- preprocessReads(faFiles[1], outFile)
  expect_equal(resSoll, res)
  unlink(outFile)
})
