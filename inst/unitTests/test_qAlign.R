source(system.file(package="QuasR", "unitTests", "help_function.R"))

test_downloadBSgenome <- function()
{
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }

    library(Rsamtools)
    
    td <- tempfile()
    checkTrue(dir.create(td, recursive=T), "Temporary directory could not be created")

    ###############################
    # Create Bam File             #
    ###############################
    
    samfile <- tempfile(tmpdir=td)
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:99\n", file=samfile, append=T)
    cat("seq1\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\n", file=samfile, append=T)
    cat("seq1\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*", file=samfile, append=T)
  
    bamfile <- asBam(samfile, samfile)
    
    ###############################
    # Create Sample File          #
    ###############################
    
    samplefile <- tempfile(tmpdir=td)
    write.table(data.frame(FileName=bamfile,
                           SampleName="Normal", stringsAsFactors=F),
                sep="\t", quote=F, row.names=F, file=samplefile)         
    
    ###############################
    # qAlign                      #
    ###############################

#     checkTrue(suppressWarnings(require("BSgenome.Dmelanogaster.UCSC.dm3", quietly=TRUE)), 
#               "BSgenome already installed in the temporary directory")
#     project <- NULL
#     project <- qAlign(samplefile, "BSgenome.Dmelanogaster.UCSC.dm3", paired="fr", clObj=clObj)
#     checkTrue(is(project, "qProject"), "qAlign (1) did not run successfully")
    
# 	checkTrue(!suppressWarnings(require("BSgenome.Ecoli.NCBI.20080805", lib.loc=td, quietly=TRUE)), 
#               "BSgenome already installed in the temporary directory")
#     project <- NULL
#     project <- qAlign(samplefile, "BSgenome.Ecoli.NCBI.20080805", paired="fr", lib.loc=td, clObj=clObj)
#     checkTrue(is(project, "qProject"), "qAlign (1) did not run successfully")
#     
#     checkTrue(require("BSgenome.Ecoli.NCBI.20080805", lib.loc=td, quietly=TRUE), 
#               "BSgenome should be installed now in the temporary directory")
#     project <- NULL
#     project <- qAlign(samplefile, "BSgenome.Ecoli.NCBI.20080805", paired="fr", lib.loc=td, clObj=clObj)
#     checkTrue(is(project, "qProject"), "qAlign (2) did not run successfully")
    
#     detach("package:BSgenome.Ecoli.NCBI.20080805", unload=T)
# 	remove.packages("BSgenome.Ecoli.NCBI.20080805", lib=td)
    unlink(td, recursive=T)
    
}

test_normal_paired <- function(){
    if(!"sampleFileGenomePaired" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePaired <<- createReads(genomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePaired, genomeFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_normal_single <- function(){
    if(!"sampleFileGenomeSingle" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomeSingle <<- createReads(genomeFile, td, paired=FALSE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingle, genomeFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_bisulfit_paired <- function(){
    if(!"sampleFileGenomePairedBisPartial" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        unMethRanges <<- scanFaIndex(genomeFile)
        partialMethRanges <<- narrow(unMethRanges, start=5000, end=c(15000,10000,15000))
        sampleFileGenomePairedBisPartial <<- createReads(genomeFile, td, paired=TRUE, bisulfit=partialMethRanges)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    ## Dir
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfit="dir", clObj=clObj)

    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")

    
    ## Undir
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfit="undir", clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")    
}

test_bisulfit_single <- function(){
    if(!"sampleFileGenomeSingleBisPartial" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        unMethRanges <<- scanFaIndex(genomeFile)
        partialMethRanges <<- narrow(unMethRanges, start=5000, end=c(15000,10000,15000))
        sampleFileGenomeSingleBisPartial <<- createReads(genomeFile, td, paired=FALSE, bisulfit=partialMethRanges)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    ## Dir
    project <- qAlign(sampleFileGenomeSingleBisPartial, genomeFile, bisulfit="dir", clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")
    
    ## Undir
    project <- qAlign(sampleFileGenomeSingleBisPartial, genomeFile, bisulfit="undir", clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(readInfo[8]=="l", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")
}

test_allelic_paired <- function(){
    if(!"sampleFileGenomePairedAllele" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        snpFile <- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
        td <<- tempdir()
        sampleFileGenomePairedAllele <<- createReads(genomeFile, td, paired=TRUE, snpFile=snpFile)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePairedAllele, genomeFile, snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T, 
                                   param=ScanBamParam(what="qname", tag="XV"))
    # extract read information from qname
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(7,8)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,7] == seqnames(aln)),
              "Test seqname")
    # load snp list
    snp <- read.delim(snpFile, header=F, colClasses=c("factor", "numeric","character","character"))
    snp <- GRanges(seqnames=snp$V1, ranges=IRanges(start=snp$V2, width=1))
    # check XV tag
    r_aln <- aln[mcols(aln)$XV == "R"]
    r_idx <- overlapsAny(r_aln, snp)
    if(!all(r_idx)){
        rSnp <- r_aln[r_idx]
        rNoSnp <- r_aln[!r_idx]
        checkTrue(all(mcols(rNoSnp)$qname %in% mcols(rSnp)$qname), 
                  "Test XV tag: All read with tag XV=R should overlap a snp or mate read should overlap a snp")
    }
    a_aln <- aln[mcols(aln)$XV == "A"]
    a_idx <- overlapsAny(a_aln, snp)
    if(!all(a_idx)){
        aSnp <- a_aln[a_idx]
        aNoSnp <- a_aln[!a_idx]
        checkTrue(all(mcols(aNoSnp)$qname %in% mcols(aSnp)$qname), 
                  "Test XV tag: All read with tag XV=A should overlap a snp or mate read should overlap a snp")
    }
    u_aln <- aln[mcols(aln)$XV == "U"]
    u_idx <- overlapsAny(u_aln, snp)
    checkTrue(all(!u_idx),
              "Test XV tag: No read with tag XV=U should overlap a snp")

    # alignments with wrong XV tag
    uSnp <- u_aln[u_idx]
}

test_allelic_single <- function(){
    if(!"sampleFileGenomeSingleAllele" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        snpFile <- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
        td <<- tempdir()
        sampleFileGenomeSingleAllele <<- createReads(genomeFile, td, paired=FALSE, snpFile=snpFile)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingleAllele, genomeFile, snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(what="qname", tag="XV"))
    # extract read information from qname
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
    
    # load snp list
    snp <- read.delim(snpFile, header=F, colClasses=c("factor", "numeric","character","character"))
    snp <- GRanges(seqnames=snp$V1, ranges=IRanges(start=snp$V2, width=1))
    # check XV tag
    r_aln <- aln[mcols(aln)$XV == "R"]
    r_idx <- overlapsAny(r_aln, snp)
    checkTrue(all(r_idx),
              "Test XV tag: All read with tag XV=R should overlap a snp")
    u_aln <- aln[mcols(aln)$XV == "U"]
    u_idx <- overlapsAny(u_aln, snp)
    checkTrue(all(!u_idx),
              "Test XV tag: No read with tag XV=U should overlap a snp")
    a_aln <- aln[mcols(aln)$XV == "A"]
    a_idx <- overlapsAny(a_aln, snp)
    checkTrue(all(a_idx),
              "Test XV tag: All read with tag XV=A should overlap a snp")
    
    # alignments with wrong XV tag
    rNoSnp <- r_aln[!r_idx]
    aNoSnp <- a_aln[!a_idx]
    uSnp <- u_aln[u_idx]
}

test_spliced_paired <- function(){
    if(!"sampleFileGenomePaired" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePaired <<- createReads(genomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePaired, genomeFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_spliced_single <- function(){
    if(!"sampleFileGenomeSingle" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomeSingle <<- createReads(genomeFile, td, paired=FALSE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingle, genomeFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_aux_normal_paired <- function(){
    if(!"sampleFileAuxPaired" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxPaired <<- createReads(auxGenomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(alnAux)=="+", readInfo[,3], readInfo[,5]) == start(alnAux)),
              "Test left read position")
    checkTrue(all(ifelse(strand(alnAux)=="+", readInfo[,4], readInfo[,6]) == end(alnAux)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")
}

test_aux_normal_single <- function(){
    if(!"sampleFileAuxSingle" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxSingle <<- createReads(auxGenomeFile, td, paired=FALSE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(alnAux)=="+", readInfo[,3], readInfo[,5]) == start(alnAux)),
              "Test left read position")
    checkTrue(all(ifelse(strand(alnAux)=="+", readInfo[,4], readInfo[,6]) == end(alnAux)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")    
}

test_aux_bisulfit_undir_paired <- function(){
    if(!"sampleFileAuxPaired" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxPaired <<- createReads(auxGenomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, bisulfite="undir", clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")             
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")
}

test_aux_bisulfit_undir_single <- function(){
    if(!"sampleFileAuxSingle" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxSingle <<- createReads(auxGenomeFile, td, paired=FALSE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, bisulfite="undir", clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")             
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")
}

test_aux_spliced_paired <- function(){
    if(!"sampleFileAuxPaired" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxPaired <<- createReads(auxGenomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")             
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")    
}

test_aux_spliced_single <- function(){
    if(!"sampleFileAuxSingle" %in% ls(envir=.GlobalEnv)){
        auxGenomeFile <<- system.file(package="QuasR", "extdata", "NC_001422.1.fa")
        td <<- tempdir()
        sampleFileAuxSingle <<- createReads(auxGenomeFile, td, paired=FALSE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
    auxFile <<- system.file(package="QuasR", "extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readBamGappedAlignments(alignments(project)$aux[1,1], use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(alnAux),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")             
    checkTrue(all(readInfo[,8] == seqnames(alnAux)),
              "Test seqname")    
}

test_normal_paired_fasta <- function(){
    if(!"sampleFileGenomePairedFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=TRUE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, clObj=clObj)

    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_normal_single_fasta <- function(){    
    if(!"sampleFileGenomeSingleFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=FALSE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_bisulfit_dir_paired_fasta <- function(){
    if(!"sampleFileGenomePairedFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=TRUE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, bisulfit="dir", clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")            
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_bisulfit_dir_single_fasta <- function(){
    if(!"sampleFileGenomeSingleFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=FALSE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, bisulfit="dir", clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")          
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_allelic_paired_fasta <- function(){
    if(!"sampleFileGenomePairedFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=TRUE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    snpFile <<- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T, 
                                   param=ScanBamParam(what="qname", tag="XV"))
    # extract read information from qname
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_allelic_single_fasta <- function(){
    if(!"sampleFileGenomeSingleFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomeSingleFasta <<- createReads(genomeFile, td, paired=FALSE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    snpFile <<- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(what="qname", tag="XV"))
    # extract read information from qname
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,3], readInfo[,5]) == start(aln)),
              "Test left read position")
    checkTrue(all(ifelse(strand(aln)=="+", readInfo[,4], readInfo[,6]) == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_spliced_paired_fasta <- function(){
    if(!"sampleFileGenomePairedFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=TRUE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_spliced_single_fasta <- function(){
    if(!"sampleFileGenomeSingleFasta" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePairedFasta <<- createReads(genomeFile, td, paired=FALSE, format="fasta")
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, splicedAlignment=TRUE, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T)
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(8,9)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")
    checkTrue(all(readInfo[,8] == seqnames(aln)),
              "Test seqname")
}

test_bisulfit_dir_allelic_paired <- function(){
    if(!"sampleFileGenomePairedAllele" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        snpFile <- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
        td <<- tempdir()
        sampleFileGenomePairedAllele <<- createReads(genomeFile, td, paired=TRUE, snpFile=snpFile)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }

    project <- qAlign(sampleFileGenomePairedAllele, genomeFile, bisulfit="dir", snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag=c("XV","NM")))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(7,8)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,7] == seqnames(aln)),
              "Test seqname")    
    # load snp list
    snp <- read.delim(snpFile, header=F, colClasses=c("factor", "numeric","character","character"))
    names(snp) <- c("seqnames","position","R","A")
    snp <- GRanges(seqnames=snp$seqnames, ranges=IRanges(start=snp$position, width=1), Ref=snp$R, Alt=snp$A)
    # check XV tag
    r_aln <- aln[mcols(aln)$XV == "R"]
    r_idx <- overlapsAny(r_aln, snp)
    if(!all(r_idx)){
        rSnp <- r_aln[r_idx]
        rNoSnp <- r_aln[!r_idx]
        checkTrue(all(mcols(rNoSnp)$qname %in% mcols(rSnp)$qname), 
                  "Test XV tag: All read with tag XV=R should overlap a snp or mate read should overlap a snp")
    }
    a_aln <- aln[mcols(aln)$XV == "A"]
    a_idx <- overlapsAny(a_aln, snp)
    if(!all(a_idx)){
        aSnp <- a_aln[a_idx]
        aNoSnp <- a_aln[!a_idx]
        checkTrue(all(mcols(aNoSnp)$qname %in% mcols(aSnp)$qname), 
                  "Test XV tag: All read with tag XV=A should overlap a snp or mate read should overlap a snp")
    }
}

test_bisulfit_dir_allelic_single <- function(){
    if(!"sampleFileGenomeSingleAllele" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        snpFile <- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
        td <<- tempdir()
        sampleFileGenomeSingleAllele <<- createReads(genomeFile, td, paired=FALSE, snpFile=snpFile)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    
    project <- qAlign(sampleFileGenomeSingleAllele, genomeFile, bisulfit="dir", snpFile=snpFile, clObj=clObj)
    
    aln <- readBamGappedAlignments(project@alignments$FileName, use.names=T,
                                   param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(7,8)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,7] == seqnames(aln)),
              "Test seqname")
    # load snp list
    snp <- read.delim(snpFile, header=F, colClasses=c("factor", "numeric","character","character"))
    names(snp) <- c("seqnames","position","R","A")
    snp <- GRanges(seqnames=snp$seqnames, ranges=IRanges(start=snp$position, width=1), Ref=snp$R, Alt=snp$A)
    # check XV tag
    r_aln <- aln[mcols(aln)$XV == "R"]
    r_idx <- overlapsAny(r_aln, snp)
    if(!all(r_idx)){
        rSnp <- r_aln[r_idx]
        rNoSnp <- r_aln[!r_idx]
        checkTrue(all(mcols(rNoSnp)$qname %in% mcols(rSnp)$qname), 
                  "Test XV tag: All read with tag XV=R should overlap a snp or mate read should overlap a snp")
    }
    a_aln <- aln[mcols(aln)$XV == "A"]
    a_idx <- overlapsAny(a_aln, snp)
    if(!all(a_idx)){
        aSnp <- a_aln[a_idx]
        aNoSnp <- a_aln[!a_idx]
        checkTrue(all(mcols(aNoSnp)$qname %in% mcols(aSnp)$qname), 
                  "Test XV tag: All read with tag XV=A should overlap a snp or mate read should overlap a snp")
    }
    u_aln <- aln[mcols(aln)$XV == "U"]
    u_idx <- overlapsAny(u_aln, snp)
    uSnp <- u_aln[u_idx]
    checkTrue(all(uSnp[strand(uSnp) == "+"] %in% snp[mcols(snp)$Ref == "C"]),
              "Test XV tag: Read with tag XV=U should overlap only snp of type C to T")
    checkTrue(all(uSnp[strand(uSnp) == "-"] %in% snp[mcols(snp)$Ref == "G"]),
              "Test XV tag: Read with tag XV=U should overlap only snp of type G to A")
}

test_maxHits <- function(){
    if(!"sampleFileGenomePaired" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        sampleFileGenomePaired <<- createReads(genomeFile, td, paired=TRUE)
    }
    if(!"clObj" %in% ls(envir=.GlobalEnv)){
        clObj <<- makeCluster(2)
    }
    project <- qAlign(sampleFileGenomePaired, genomeFile, maxHits=100, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])

    snpFile <- system.file(package="QuasR", "extdata", "hg19sub_snp.txt")
    project <- qAlign(sampleFileGenomePaired, genomeFile, snpFile=snpFile, maxHits=100, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])
    
    if(!"sampleFileGenomePairedBisPartial" %in% ls(envir=.GlobalEnv)){
        genomeFile <<- system.file(package="QuasR", "extdata", "hg19sub.fa")
        td <<- tempdir()
        unMethRanges <<- scanFaIndex(genomeFile)
        partialMethRanges <<- narrow(unMethRanges, start=5000, end=c(15000,10000,15000))
        sampleFileGenomePairedBisPartial <<- createReads(genomeFile, td, paired=TRUE, bisulfit=partialMethRanges)
    }
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfit="undir", maxHits=100, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])
}

# Not supported
# test_aux_allelic_paired <- function(){}
# test_aux_allelic_single <- function(){}
# test_aux_bisulfit_dir_allelic_paired <- function(){}
# test_aux_bisulfit_dir_allelic_single <- function(){}
# test_aux_bisulfit_undir_allelic_paired <- function(){}
# test_aux_bisulfit_undir_allelic_single <- function(){}
# test_aux_spliced_allelic_paired <- function(){}
# test_aux_spliced_allelic_single <- function(){}

# Not tested yet
# test_bisulfit_undir_allelic_paired <- function(){}
# test_bisulfit_undir_allelic_single <- function(){}
# test_spliced_allelic_paired <- function(){}
# test_spliced_allelic_single <- function(){}
# test_aux_bisulfit_dir_paired <- function(){}
# test_aux_bisulfit_dir_single <- function(){}
# test_normal_paired_ff <- function(){}
# test_normal_paired_rf <- function(){}
