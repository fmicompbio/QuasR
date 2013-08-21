# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!existsFunction("createFastaReads"))
    source(system.file(package="QuasR", "unitTests", "help_function.R"))

if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

genomeFile <- file.path("extdata", "hg19sub.fa")
snpFile <- file.path("extdata", "hg19sub_snp.txt")
auxGenomeFile <- file.path("extdata", "NC_001422.1.fa")
td <- tempdir()

sampleFileGenomeSingle <- createReads(genomeFile, td, paired=FALSE)
sampleFileGenomePaired <- createReads(genomeFile, td, paired=TRUE)

unMethRanges <- scanFaIndex(genomeFile)
partialMethRanges <- narrow(unMethRanges, start=5000, end=c(15000,10000,15000))

sampleFileGenomeSingleBisPartial <- createReads(genomeFile, td, paired=FALSE, bisulfite=partialMethRanges)
sampleFileGenomePairedBisPartial <- createReads(genomeFile, td, paired=TRUE, bisulfite=partialMethRanges)

sampleFileGenomeSingleAllele <- createReads(genomeFile, td, paired=FALSE, snpFile=snpFile)
sampleFileGenomePairedAllele <- createReads(genomeFile, td, paired=TRUE, snpFile=snpFile)

sampleFileAuxSingle <- createReads(auxGenomeFile, td, paired=FALSE)
sampleFileAuxPaired <- createReads(auxGenomeFile, td, paired=TRUE)

sampleFileGenomeSingleFasta <- createReads(genomeFile, td, paired=FALSE, format="fasta")
sampleFileGenomePairedFasta <- createReads(genomeFile, td, paired=TRUE, format="fasta")

#.setUp <- function() { # runs before each test_...()
#    # make sure clObj exists and is a working cluster object
#    if(!exists("clObj", envir=.GlobalEnv) ||
#       !inherits(clObj, "cluster") ||
#       inherits(try(all(unlist(clusterEvalQ(clObj, TRUE))), silent=TRUE), "try-error")) {
#        clObj <<- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
#        clusterEvalQ(clObj, library("QuasR"))
#    }
#}

# tests sporadically fail if clObj is not freshly generated for each test_...()
.setUp <- function() {
    if(exists("clObj", envir=.GlobalEnv) && inherits(clObj, "cluster") &&
       !inherits(try(all(unlist(clusterEvalQ(clObj, TRUE))), silent=TRUE), "try-error")) {
        stopCluster(get("clObj", envir=.GlobalEnv))
    }
    clObj <<- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
}

.tearDown <- function() {
    stopCluster(get("clObj", envir=.GlobalEnv))
    rm(clObj, envir=.GlobalEnv)
}


test_normal_paired <- function(){
    project <- qAlign(sampleFileGenomePaired, genomeFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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
    project <- qAlign(sampleFileGenomeSingle, genomeFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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

test_bisulfite_paired_dir <- function(){
    ## Dir
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfite="dir", alignmentsDir=td, clObj=clObj)

    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
                                  param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")            
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")
}

test_bisulfite_paired_undir <- function(){
    ## Undir
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfite="undir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
                                  param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
              "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")              
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")    
}

test_bisulfite_single_dir <- function(){
    ## Dir
    project <- qAlign(sampleFileGenomeSingleBisPartial, genomeFile, bisulfite="dir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
                                  param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkVar <- readInfo[,3] == start(aln) | readInfo[,5] == start(aln)
    checkTrue(all(checkVar), 
              paste("Test left read position:", readInfo[!checkVar,8][1],
                    readInfo[!checkVar,3][1],readInfo[!checkVar,5][1], start(aln)[!checkVar][1]))
    #checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
    #          "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")            
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")
}

test_bisulfite_single_undir <- function(){
    ## Undir
    project <- qAlign(sampleFileGenomeSingleBisPartial, genomeFile, bisulfite="undir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
                                  param=ScanBamParam(tag="NM"))
    readInfo <- as.data.frame(do.call(rbind, strsplit(names(aln),"_")), stringsAsFactors=F)
    readInfo[,c(9,10,11)] <- do.call(rbind,strsplit(readInfo[,2], "-"))
    # check start, end and seqname
    checkVar <- readInfo[,3] == start(aln) | readInfo[,5] == start(aln)
    checkTrue(all(checkVar), 
              paste("Test left read position:", readInfo[!checkVar,8][1],
                    readInfo[!checkVar,3][1],readInfo[!checkVar,5][1], start(aln)[!checkVar][1]))
    #checkTrue(all(readInfo[,3] == start(aln) | readInfo[,5] == start(aln)),
    #          "Test left read position")
    checkTrue(all(readInfo[,4] == end(aln) | readInfo[,6] == end(aln)),
              "Test right read position")           
    checkTrue(all(readInfo[,9] == seqnames(aln)),
              "Test seqname")
    
    # check if strand corrspond to converted genome
    nm0 <- mcols(aln)$NM == 0
    checkTrue(all(ifelse(as.vector(strand(aln[nm0])=="+"), "CtoT", "GtoA") == readInfo[nm0,11]),
              "Test if strand correspond to the coverted genome")
}

test_allelic_paired <- function(){
    project <- qAlign(sampleFileGenomePairedAllele, genomeFile, snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T, 
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
    project <- qAlign(sampleFileGenomeSingleAllele, genomeFile, snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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
    project <- qAlign(sampleFileGenomePaired, genomeFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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
    project <- qAlign(sampleFileGenomeSingle, genomeFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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

test_aux_bisulfite_undir_paired <- function(){
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, bisulfite="undir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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

test_aux_bisulfite_undir_single <- function(){
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, bisulfite="undir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxPaired, genomeFile, auxiliaryFile=auxFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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
    #genomeFile <<- file.path("extdata", "hg19sub.fa")
    auxFile <<- file.path("extdata", "auxiliaries.txt")
    project <- qAlign(sampleFileAuxSingle, genomeFile, auxiliaryFile=auxFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(alignments(project)$genome$FileName, use.names=T)
    alnAux <- readGAlignmentsFromBam(alignments(project)$aux[1,1], use.names=T)
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
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, alignmentsDir=td, clObj=clObj)

    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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

test_bisulfite_dir_paired_fasta <- function(){
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, bisulfite="dir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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

test_bisulfite_dir_single_fasta <- function(){
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, bisulfite="dir", alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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
    snpFile <<- file.path("extdata", "hg19sub_snp.txt")
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T, 
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
    snpFile <<- file.path("extdata", "hg19sub_snp.txt")
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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
    project <- qAlign(sampleFileGenomePairedFasta, genomeFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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
    project <- qAlign(sampleFileGenomeSingleFasta, genomeFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T)
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

test_bisulfite_dir_allelic_paired <- function(){
    project <- qAlign(sampleFileGenomePairedAllele, genomeFile, bisulfite="dir", snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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

test_bisulfite_dir_allelic_single <- function(){
    project <- qAlign(sampleFileGenomeSingleAllele, genomeFile, bisulfite="dir", snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    
    aln <- readGAlignmentsFromBam(project@alignments$FileName, use.names=T,
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
    checkTrue(all(overlapsAny(uSnp[strand(uSnp) == "+"], snp[mcols(snp)$Ref == "C"])),
              "Test XV tag: Read with tag XV=U should overlap only snp of type C to T")
    checkTrue(all(overlapsAny(uSnp[strand(uSnp) == "-"], snp[mcols(snp)$Ref == "G"])),
              "Test XV tag: Read with tag XV=U should overlap only snp of type G to A")
}

test_maxHits_simple <- function(){
    project <- qAlign(sampleFileGenomePaired, genomeFile, maxHits=100, alignmentsDir=td, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])
}

test_maxHits_snps <- function(){
    project <- qAlign(sampleFileGenomePaired, genomeFile, snpFile=snpFile, maxHits=100, alignmentsDir=td, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])
}

test_maxHits_bisulfite <- function(){
    project <- qAlign(sampleFileGenomePairedBisPartial, genomeFile, bisulfite="undir", maxHits=100, alignmentsDir=td, clObj=clObj)
    checkTrue(0 == alignmentStats(project)[,"unmapped"])
}

# Not supported
# test_aux_allelic_paired <- function(){}
# test_aux_allelic_single <- function(){}
# test_aux_bisulfite_dir_allelic_paired <- function(){}
# test_aux_bisulfite_dir_allelic_single <- function(){}
# test_aux_bisulfite_undir_allelic_paired <- function(){}
# test_aux_bisulfite_undir_allelic_single <- function(){}
# test_aux_spliced_allelic_paired <- function(){}
# test_aux_spliced_allelic_single <- function(){}

# Not tested yet
# test_bisulfite_undir_allelic_paired <- function(){}
# test_bisulfite_undir_allelic_single <- function(){}
# test_spliced_allelic_paired <- function(){}
# test_spliced_allelic_single <- function(){}
# test_aux_bisulfite_dir_paired <- function(){}
# test_aux_bisulfite_dir_single <- function(){}
# test_normal_paired_ff <- function(){}
# test_normal_paired_rf <- function(){}
