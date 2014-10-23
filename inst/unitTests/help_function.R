### helper functions used by QuasR's unit test functions

# Fragment design
# Left Adapter  ATACTG
# Insert        TGTGACAGACCTCGGGGCCACATGCACTGACTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC
# Right Adapter ATCTCGTATGCCGTCTTCTGCTTG
createFastaReads <- function(){ # create two temporary fasta files and return file names
    faFileName1 <- tempfile(fileext=".fa")
    faFile1 <-file(faFileName1, open="w")
    writeLines(c(">seq1", "CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", # onlyInsert
                 ">seq2", "AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", # low compexity, 2N
                 ">seq3", "AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", # A.., LAdapter
                 ">seq4", "ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", # LAdapter
                 ">seq5", "TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", # Insert
                 ">seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # LAdapter, noInsert, RAdapter, A..
                 ">seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", # LAdapter, noInsert, RAdapter, A.., 2N
                 ">seq8", "ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", # LAdapter, shortInsert, RAdapter
                 ">seq9", "TGACAGACCTCGGGGCCAC"),                 # shortLength
               con=faFile1)
    close(faFile1)
    
    faFileName2 <- tempfile(fileext=".fa")
    faFile2 <-file(faFileName2, open="w")
    writeLines(c(">seq1", "AAAAAAAAAAAAAAAAAAA",                  # low compexity
                 ">seq2", "TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", # low compexity, shortLength
                 ">seq3", "CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", # onlyInsert
                 ">seq4", "GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # RAdapter
                 ">seq5", "AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", # Insert
                 ">seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", # LAdapter, noInsert, RAdapter, A..
                 ">seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", # LAdapter, noInsert, RAdapter, A.., 2N
                 ">seq8", "AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", # A.., LAdapter, shortInsert, RAdapter
                 ">seq9", "CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT"),# Insert
               con=faFile2)
    close(faFile1)
    return(c(faFileName1, faFileName2))
}

createFastqReads <- function(){ # create two temporary fastq files and return file names
    fqFileName1 <- tempfile(fileext=".fq")
    fqFile1 <- file(fqFileName1, open="w")
    writeLines(c("@seq1", "CCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGA", "+", "5..49<494*<49493####################",
                 "@seq2", "AAAAANNAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "+", "+4544944444444######################",
                 "@seq3", "AAAATACTGTGTGACAGACCTCGGGGCCACATGCAC", "+", "BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@/=",
                 "@seq4", "ATACTGTGTGACAGACCTCGGGGCCACATGCACTGA", "+", "BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABC",
                 "@seq5", "TGTGACAGACCTCGGGGCCACATGCACTGACTCCTC", "+", "BA?AB??60>A6?0BBBB95>057;@*.<(434;+4",
                 "@seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>BBB@",
                 "@seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", "+", "####################################",
                 "@seq8", "ATACTGGAGGTCATCTCGTATGCCGTCTTCTGCTTG", "+", "(33(;?B@AB43*,/;9(6</7>5;<##########",
                 "@seq9", "TGACAGACCTCGGGGCCAC",                  "+", "###################"),
               con=fqFile1)
    close(fqFile1)
    
    fqFileName2 <- tempfile(fileext=".fq")
    fqFile2 <- file(fqFileName2, open="w")    
    writeLines(c("@seq1", "AAAAAAAAAAAAAAAAAAA",                  "+", "3@?3/>9A8@A1-*/4@BB",
                 "@seq2", "TCAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGG", "+", "4-093<##############################",
                 "@seq3", "CTCCTCAGCTGCCAGATGTGCAGTCCAAGCTGGGCC", "+", "CB?(8=(<A/<=-(07+7&883@#############",
                 "@seq4", "GAGGTCATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "5..49<494*<4949#####################",
                 "@seq5", "AGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGTC", "+", "+45449444444########################",
                 "@seq6", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAAAA", "+", "BCCBCCBBB60>CA;5;@BB@A6+;8@BC?0:B@A<",
                 "@seq7", "ATACTGATCTCGTATGCCGTCTTCTGCTTGAAAANN", "+", "BCBBB>ACBCCCBCC@BCC@*7@82=BBBB1>CABB",
                 "@seq8", "AAAAATACTGGAGGTCATCTCGTATGCCGTCTTCTG", "+", "BA?AB??60>A6?0BBBB95>057;@*.<(434;@B",
                 "@seq9", "CAGCTGCCAGATGTGCAGTCCAAGCTGGGCCGAGGT", "+", "+BCCCCBCCCBACB:?BBCCCCCCBCBC>;(>B@4;"),
               con=fqFile2)
    close(fqFile2)
    return(c(fqFileName1, fqFileName2))
}

createProjectPairedMode <- function(){
    require(Rsamtools)

    # Create Samfile without XV tag
    samfile <- tempfile()
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:99\n", file=samfile, append=T)
    # fr R1->left R2->right
    cat("seq1\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\n", file=samfile, append=T)
    cat("seq1\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\n", file=samfile, append=T)
    cat("seq2\t99\tchrV\t4\t255\t5M\tchrV\t13\t14\t*\t*\n", file=samfile, append=T)
    cat("seq2\t147\tchrV\t13\t255\t5M\tchrV\t4\t-14\t*\t*\n", file=samfile, append=T)
    # fr R2->left R1->right
    cat("seq3\t163\tchrV\t24\t255\t5M\tchrV\t32\t13\t*\t*\n", file=samfile, append=T)
    cat("seq3\t83\tchrV\t32\t255\t5M\tchrV\t24\t-13\t*\t*\n", file=samfile, append=T)
    cat("seq4\t163\tchrV\t24\t255\t5M\tchrV\t33\t14\t*\t*\n", file=samfile, append=T)
    cat("seq4\t83\tchrV\t33\t255\t5M\tchrV\t24\t-14\t*\t*\n", file=samfile, append=T)
    # ff R1->left R2->right
    cat("seq5\t67\tchrV\t44\t255\t5M\tchrV\t52\t13\t*\t*\n", file=samfile, append=T)
    cat("seq5\t131\tchrV\t52\t255\t5M\tchrV\t44\t-13\t*\t*\n", file=samfile, append=T)
    cat("seq6\t67\tchrV\t44\t255\t5M\tchrV\t53\t14\t*\t*\n", file=samfile, append=T)
    cat("seq6\t131\tchrV\t53\t255\t5M\tchrV\t44\t-14\t*\t*\n", file=samfile, append=T)
    # rr R2->left R1->right
    cat("seq7\t115\tchrV\t64\t255\t5M\tchrV\t72\t13\t*\t*\n", file=samfile, append=T)
    cat("seq7\t179\tchrV\t72\t255\t5M\tchrV\t64\t-13\t*\t*\n", file=samfile, append=T)
    cat("seq8\t115\tchrV\t64\t255\t5M\tchrV\t73\t14\t*\t*\n", file=samfile, append=T)
    cat("seq8\t179\tchrV\t73\t255\t5M\tchrV\t64\t-14\t*\t*\n", file=samfile, append=T)
    # rf R1->left R2->right
    cat("seq9\t83\tchrV\t84\t255\t5M\tchrV\t92\t13\t*\t*\n", file=samfile, append=T)
    cat("seq9\t163\tchrV\t92\t255\t5M\tchrV\t84\t-13\t*\t*\n", file=samfile, append=T)
    cat("seq0\t83\tchrV\t84\t255\t5M\tchrV\t93\t14\t*\t*\n", file=samfile, append=T)
    cat("seq0\t163\tchrV\t93\t255\t5M\tchrV\t84\t-14\t*\t*", file=samfile, append=T)
    
    # Create Bamfile
    bamfile <- asBam(samfile, samfile)

    # Create Sample File
    samplefile <- tempfile()
    write.table(data.frame(FileName=bamfile,
                           SampleName="Normal", stringsAsFactors=F),
                sep="\t", quote=F, row.names=F, file=samplefile)
    
    # Create Genome
    genome <- tempfile(fileext=".fa")
    cat(">chrV\n", file=genome, append=FALSE)
    cat(paste(rep("G",99), collapse=""), file=genome, append=T)
    
    # Create SNP file
    snp <- tempfile(fileext=".txt")
    cat(">chrV\t10\tC\tG", file=snp, append=FALSE)

    # qAlign
    td <- tempdir()
    project <- qAlign(samplefile, genome, paired="fr", alignmentsDir=td)
    
    return(project)
}

createProjectAllelic <- function(){
    require(Rsamtools)

    # Create samfile with XV tag
    samfile <- tempfile()
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:99\n", file=samfile, append=T)
    # fr XV=R
    cat("seq1\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:R\n", file=samfile, append=T)
    cat("seq1\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:R\n", file=samfile, append=T)
    # fr XV=U
    cat("seq2\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:U\n", file=samfile, append=T)
    cat("seq2\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:U\n", file=samfile, append=T)
    # fr XV=A
    cat("seq3\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:A\n", file=samfile, append=T)
    cat("seq3\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:A", file=samfile, append=T)

    # Create Bamfile
    bamfile <- asBam(samfile, samfile)
    
    # Create Sample File
    samplefile <- tempfile()
    write.table(data.frame(FileName=bamfile,
                           SampleName="Allelic", stringsAsFactors=F),
                sep="\t", quote=F, row.names=F, file=samplefile)
    
    # Create Genome
    genome <- tempfile(fileext=".fa")
    cat(">chrV\n", file=genome, append=FALSE)
    cat(paste(rep("G",99), collapse=""), file=genome, append=T)
    
    # Create SNP file
    snp <- tempfile(fileext=".txt")
    cat(">chrV\t10\tC\tG", file=snp, append=FALSE)     
    
    # 
    td <- tempdir()
    project <- qAlign(samplefile, genome, snpFile=snp, paired="fr", alignmentsDir=td)
    
    return(project)
}

createProjectSingleMode <- function(allelic=FALSE){
    require(Rsamtools)
    
    # Create Bam File
    samfile_plus <- tempfile()
    samfile_minus <- tempfile()
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile_plus, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:800", file=samfile_plus, append=T)
    
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile_minus, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:800", file=samfile_minus, append=T)
    
    for(pos in 1:791){
        cat("\n", file=samfile_plus, append=T)
        cat("\n", file=samfile_minus, append=T)
        cat(paste("seq1\t0\tchrV", pos, pos %% 256, "10M", "*\t0\t0\t*\t*", c("XV:A:R","XV:A:A","XV:A:U"), 
                  sep="\t", collapse="\n"), 
            file=samfile_plus, append=T)
        cat(paste("seq1\t16\tchrV", pos, pos %% 256, "10M", "*\t0\t0\t*\t*", c("XV:A:R","XV:A:A","XV:A:U"),
                  sep="\t", collapse="\n"), 
            file=samfile_minus, append=T)            
    }
    bamfile_plus <- asBam(samfile_plus, samfile_plus)
    bamfile_minus <- asBam(samfile_minus, samfile_minus)
    
    # Create Genome
    genome <- tempfile(fileext=".fa")
    cat(">chrV\n", file=genome, append=FALSE)
    cat(paste(rep("G",800), collapse=""), file=genome, append=T)
    
    # Create SNP file 
    snp <- tempfile(fileext=".txt")
    cat(">chrV\t10\tC\tG", file=snp, append=FALSE)
    
    # Create Sample File
    samplefile <- tempfile()
    write.table(data.frame(FileName=c(bamfile_plus, bamfile_minus),
                           SampleName=c("Sample", "Sample"), stringsAsFactors=F),
                sep="\t", quote=F, row.names=F, file=samplefile)            
    
    # qAlign
    td <- tempdir()
    if(!allelic)
        project <- qAlign(samplefile, genome, paired="no", alignmentsDir=td)
    else
        project <- qAlign(samplefile, genome, snpFile=snp, paired="no", alignmentsDir=td)        
    return(project)
}

createTilingRegion <- function()
{
    #        101  151  201  251  301  351  401  451  501  551  601  651  701     X = 10x 3reads(=R,U,A) 
    #   H1 :   XXXXXXXXXX          XXXXXXXXXX          XXXXXXXXXX            --> 30X => 900
    #   H2 :        XXXXXXXXXX          XXXXXXXXX           XXXXXXXXXX       --> 30X => 900
    #   H3 :             XXXXXXXXXX          XXXXXXXXXX          XXXXXXXXXX  --> 30X => 900
    #   H4 :                  XXXXXXXXXX          XXXXXXXXXX          XXXXX  --> 24X => 750
    
    width=100L
    gap=100L
    times=600L/(width+gap)
    h1 <- GRanges("chrV", successiveIRanges(rep(width,times),gap, from=101))
    names(h1) <- rep("H1", length(h1))
    h2 <- GRanges("chrV", successiveIRanges(rep(width,times),gap, from=151))
    names(h2) <- rep("H2", length(h2))
    h3 <- GRanges("chrV", successiveIRanges(rep(width,times),gap, from=201))
    names(h3) <- rep("H3", length(h3))
    h4 <- GRanges("chrV", successiveIRanges(rep(width,times),gap, from=251))
    names(h4) <- rep("H4", length(h4))
    end(h4[3]) <- 700
    tilingRegion <- sort(c(h1, h2, h3, h4))    
    return(tilingRegion)
}

createGtfRegion <- function()
{
    require("rtracklayer")
    annotationFile <- file.path("extdata", "hg19sub_annotation.gtf")
    gtfRegion <- import.gff(annotationFile, format="gtf", asRangedData=F,
                            feature.type="exon")
    names(gtfRegion) <- mcols(gtfRegion)$gene_name
    return(gtfRegion)
}

createGenomeRegion <- function()
{
    genomeFile <- file.path("extdata", "hg19sub.fa")
    genomeRegion <- scanFaIndex(genomeFile)
    return(genomeRegion)
}

createAuxRegion <- function()
{
    auxGenomeFile <- file.path("extdata", "NC_001422.1.fa")
    auxRegion <- scanFaIndex(auxGenomeFile)
    return(auxRegion)
}

createTranscriptDb <- function()
{
    require("GenomicFeatures")
    annotationFile <- file.path("extdata", "hg19sub_annotation.gtf")
    genomeRegion <- createGenomeRegion()
    chrominfo <- data.frame(chrom=as.character(seqnames(genomeRegion)),
                            length=end(genomeRegion),
                            is_circular=rep(FALSE, length(genomeRegion)))
    txdb <- makeTxDbFromGFF(annotationFile, 
                            format="gtf", 
                            exonRankAttributeName="exon_number",  
                            gffGeneIdAttributeName="gene_name", 
                            chrominfo=chrominfo,
                            dataSource="Ensembl",
                            species="Homo sapiens")
    return(txdb)
}

createReads <- function(genomeFile, destDir, paired=TRUE, snpFile=NULL, bisulfite=NULL, transcript=NULL, format="fastq"){
    library(Rsamtools)
    genome <- scanFa(genomeFile)
    destination <- tempfile(tmpdir=destDir)
    filename1 <- ifelse(format=="fastq", 
                        paste(destination,"_1.fq", sep=""),
                        paste(destination,"_1.fa", sep=""))
    filename2 <- ifelse(format=="fastq", 
                        paste(destination,"_2.fq", sep=""),
                        paste(destination,"_2.fa", sep=""))
    if(paired){
        suppressWarnings(file.remove(c(filename1, filename2)))
    }else{
        suppressWarnings(file.remove(filename1))
    }
    
    if(!is.null(snpFile)){
        snps <- read.table(snpFile, colClasses=c("factor","numeric","character","character"))
        snps[,3] <- toupper(snps[,3]) # convert ref nucleotides to upper case if not already the case
        snps[,4] <- toupper(snps[,4]) # convert alt nucleotides to upper case if not already the case
        genomeSNP <- lapply(split(snps, snps[,1]), function(s){
            chr <- as.character(s[1,1])
            as.character(replaceLetterAt(genome[[chr]], s[,2], s[,4]))
        })
        genomeSNP <- DNAStringSet(unlist(genomeSNP, use.names=T))
        names(genomeSNP) <- paste(names(genomeSNP), "A", sep="-")
        names(genome) <- paste(names(genome), "R", sep="-")
        genome <- append(genome, genomeSNP)
    } else {
        names(genome) <- paste(names(genome), "R", sep="-")
    }
    
    if(!is.null(bisulfite)){
        CpG <- lapply(seq_along(genome), function(i){
            seq <- genome[[i]]
            seqnames <- unlist(strsplit(names(genome[i]), "-"))[1]
            ranges <- ranges(bisulfite[seqnames(bisulfite) == seqnames])
            mask <- Mask(mask.width=length(seq), start=start(ranges), end=end(ranges))
            masks(seq) <- mask
            start(matchPattern("CG", seq))
        })
        names(CpG) <- names(genome)
        
        genomeCtoT <- chartr("C","T", genome)
        genomeGtoA <- chartr("G","A", genome)
        genome <- lapply(seq_along(CpG), function(i){
            seqCtoT <- genomeCtoT[[names(CpG)[i]]]
            seqCtoT[CpG[[i]]] <- "C"
            seqCtoT <- as(seqCtoT, "XStringSet")
            names(seqCtoT) <- paste(names(CpG)[i], "CtoT", sep="-")
            
            seqGtoA <- genomeGtoA[[names(CpG)[i]]]
            seqGtoA[CpG[[i]]+1] <- "G"
            seqGtoA <- as(seqGtoA, "XStringSet")
            names(seqGtoA) <- paste(names(CpG)[i], "GtoA", sep="-")
            
            append(seqCtoT, seqGtoA)                        
        })
        genome <- do.call(c, genome)
    }
    
    set.seed(111)   
    sampleRange <- 5L
    width=50L + sampleRange
    gap=10L #-width+coverage
    insertSize <- 90L
    append <- FALSE
    
    lapply(seq_along(genome), function(i){
        nreads=floor(length(genome[[i]])/(width+gap))
        
        # create read ranges
        r1 <- successiveIRanges(rep(width, nreads), gapwidth=gap)
        start(r1) <- start(r1) + sample(length(r1))%%sampleRange
        r2 <- successiveIRanges(rep(width, nreads), gapwidth=gap, from=width+insertSize)
        inGenome <- end(r2) < length(genome[[i]])
        # extract and write sequence
        read1 <- Views(genome[[i]], r1[inGenome])
        read2 <- Views(genome[[i]], r2[inGenome])
        ids <- paste("seq", names(genome)[i], start(read1), end(read1), start(read2), end(read2), sep="_")
        read1 <- as(read1,"XStringSet")
        read2 <- reverseComplement(as(read2,"XStringSet"))
        # swap half of the pairs
        if(!is.null(bisulfite)){
            names(read1) <- paste(ids, "b", "l", sep="_")
            names(read2) <- paste(ids, "b", "r", sep="_")
        }        
        swap <- sample(length(read1))%%2
        sread1 <- append(read1[swap!=1], read2[swap==1])
        sread2 <- append(read2[swap!=1], read1[swap==1])
        sids <- c(ids[swap!=1], ids[swap==1])
        # add qname
        if(is.null(bisulfite)){
            names(sread1) <- paste(sids, swap, "/1", sep="_")
            names(sread2) <- paste(sids, swap, "/2", sep="_")
        }
        if(paired){
            writeXStringSet(sread1, filename1, append=TRUE, format=format)
            writeXStringSet(sread2, filename2, append=TRUE, format=format)
        }else{
            writeXStringSet(sread1, filename1, append=TRUE, format=format)
            writeXStringSet(sread2, filename1, append=TRUE, format=format)            
        }
    })
    
    samplefile <- paste(destination, "sample.txt", sep="")
    
    if(!paired){
        write.table(data.frame(FileName=filename1,
                               SampleName="Normal", stringsAsFactors=F),
                    sep="\t", quote=F, row.names=F, file=samplefile)        
    } else {
        write.table(data.frame(FileName1=filename1,
                               FileName2=filename2,
                               SampleName="Normal", stringsAsFactors=F),
                    sep="\t", quote=F, row.names=F, file=samplefile)        
    }
    
    return(samplefile)
}
