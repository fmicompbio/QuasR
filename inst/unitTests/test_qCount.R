# initialization of QuasR test environment
# allows: runTestFile("test_file.R", rngKind="default", rngNormalKind="default", verbose=1L)
if(!existsFunction("createFastaReads"))
    source(system.file(package="QuasR", "unitTests", "help_function.R"))

if(!file.exists("./extdata"))
    file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

projectSingle <- createProjectSingleMode()
projectPaired <- createProjectPairedMode()
projectAllelic <- createProjectAllelic()
projectSingleAllelic <- createProjectSingleMode(allelic=TRUE)
tilingRegion <- createTilingRegion()
gtfRegion <- createGtfRegion()
td <- tempdir()
genomeFile <- file.path("extdata", "hg19sub.fa")
sampleFile <- file.path("extdata", "samples_rna_single.txt")
snpFile <- file.path("extdata", "hg19sub_snp.txt")

.setUp <- function() { # runs before each test_...()
    # make sure clObj exists and is a working cluster object
    if(!exists("clObj", envir=.GlobalEnv) ||
       !inherits(clObj, "cluster") ||
       inherits(try(all(unlist(clusterEvalQ(clObj, TRUE))), silent=TRUE), "try-error")) {
        clObj <<- makeCluster(getOption("QuasR_nb_cluster_nodes",2))
        clusterEvalQ(clObj, library("QuasR"))
    }
}


## Test alignment selection based on mapping quality
test_mapq <- function() {
  project <- projectSingle
  
  query <- GRanges(c("chrV"), IRanges(start=1, width=800))

  checkException(qCount(project, query, mapqMin=-1))
  checkException(qCount(project, query, mapqMax=256))
  
  mapq <- scanBam(alignments(project)$genome$FileName[1])[[1]]$mapq
  mapq[is.na(mapq)] <- 255
  
  for(mymin in seq(10,250,by=40))
    checkTrue(qCount(project[1], query, mapqMin=mymin)[1,2] == sum(mapq >= mymin))
  
  for(mymax in seq(10,250,by=40))
    checkTrue(qCount(project[1], query, mapqMax=mymax)[1,2] == sum(mapq <= mymax))
  
  for(mycut in seq(10,250,by=40))
    checkIdentical(length(mapq) - qCount(project[1], query, mapqMax=mycut)[1,2],
                   qCount(project[1], query, mapqMin=mycut+1)[1,2])
  
}

## Test parameter shift, selectReadPosition
test_shift <- function() {
    project <- projectPaired
    
    ## qCount with SmartShift
    #fr R1->left R2->right
    query <- GRanges(c("chrV"), IRanges(start=1:20, width=1), "+")
    resSoll <- rep(0,20)
    resSoll[10] <- 4
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 1: qCount with smartShift")
    
    resSoll <- rep(0,20)
    resSoll[c(6,14)] <- 2       
    res <- qCount(project, query, selectReadPosition="end", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 2: qCount with smartShift")
    
    #fr R2->left R1->right
    query <- GRanges("chrV", IRanges(start=21:40, width=1), "+")
    resSoll <- rep(0,20)
    resSoll[10] <- 4
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 3: qCount with smartShift")
    
    resSoll <- rep(0,20)
    resSoll[c(6,14)] <- 2       
    res <- qCount(project, query, selectReadPosition="end", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 4: qCount with smartShift")
    
    #ff R1->left R2->right
    query <- GRanges("chrV", IRanges(start=41:60, width=1), "+")
    resSoll <- rep(0,20)
    resSoll[c(6,10)] <- 2
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 5: qCount with smartShift")
    
    resSoll <- rep(0,20)
    resSoll[c(10,14)] <- 2       
    res <- qCount(project, query, selectReadPosition="end", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 6: qCount with smartShift")
    
    #rr R2->left R1->right
    query <- GRanges("chrV", IRanges(start=61:80, width=1), "+")
    resSoll <- rep(0,20)
    resSoll[c(10,14)] <- 2
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 7: qCount with smartShift")
    
    resSoll <- rep(0,20)
    resSoll[c(6,10)] <- 2       
    res <- qCount(project, query, selectReadPosition="end", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 8: qCount with smartShift")
    
    #rf R1->left R2->right
    query <- GRanges("chrV", IRanges(start=81:99, width=1), "+")
    resSoll <- rep(0,19)
    resSoll[c(6,14)] <- 2
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 9: qCount with smartShift")
    
    resSoll <- rep(0,19)
    resSoll[10] <- 4      
    res <- qCount(project, query, selectReadPosition="end", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 10: qCount with smartShift")
    
    
    ## qCount with interger as shift
    aln <- GenomicAlignments::readGAlignmentsFromBam(project@alignments$FileName)
    
    query <- GRanges(c("chrV"), IRanges(start=1:99, width=1), "+")
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln), end(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 1: qCount with shift and selectReadPosition")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln)+1, end(aln)-1))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=1, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 2: qCount with shift and selectReadPosition")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln)-1, end(aln)+1))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=-1, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 3: qCount with shift and selectReadPosition")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln), start(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 4: qCount with shift and selectReadPosition")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln)+1, start(aln)-1))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=1, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 5: qCount with shift and selectReadPosition")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln)-1, start(aln)+1))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=-1, orientation="any")[,-1]
    checkTrue(all(resSoll == res), "Test 6: qCount with shift and selectReadPosition")
}

test_shift_allelic <- function(){
    project <- projectAllelic
    query <- GRanges(c("chrV"), IRanges(start=1:20, width=1), "+")
    
    # no shift
    resSoll <- rep(0,20)
    resSoll[c(4,16)] <- 1
    res <- qCount(project, query, selectReadPosition="start", orientation="any")[,-1]
    checkTrue(all(resSoll == res[,1]), "Test 1: qCount allele specific")
    checkTrue(all(resSoll == res[,2]), "Test 2: qCount allele specific")
    checkTrue(all(resSoll == res[,3]), "Test 3: qCount allele specific")
    colname <- paste(rep(project@alignments$SampleName, each=3), c("R","U","A"), sep="_")
    checkTrue(all(colname == colnames(res)), "Test 4: qCount allele specific")
    
    # smart shift
    resSoll <- rep(0,20)
    resSoll[10] <- 2
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", orientation="any")[,-1]
    checkTrue(all(resSoll == res[,1]), "Test 5: qCount allele specific")
    checkTrue(all(resSoll == res[,2]), "Test 6: qCount allele specific")
    checkTrue(all(resSoll == res[,3]), "Test 7: qCount allele specific")
    
    # shift
    resSoll <- rep(0,20)
    resSoll[c(6,14)] <- 1
    res <- qCount(project, query, selectReadPosition="start", shift=2, orientation="any")[,-1]
    checkTrue(all(resSoll == res[,1]), "Test 8: qCount allele specific")
    checkTrue(all(resSoll == res[,2]), "Test 9: qCount allele specific")
    checkTrue(all(resSoll == res[,3]), "Test 10: qCount allele specific")
}
   
test_orientation <- function() {
    project <- projectPaired
    aln <- GenomicAlignments::readGAlignmentsFromBam(project@alignments$FileName)
    
    query <- GRanges(c("chrV"), IRanges(start=1:99, width=1), "+")
    resSoll <- rep(0,99)
    pos <- Rle(start(aln[strand(aln)=="+"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 1: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(end(aln[strand(aln)=="-"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 2: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(end(aln[strand(aln)=="+"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 3: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(start(aln[strand(aln)=="-"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 4: qCount with orientation and query strand")
    
    query <- GRanges(c("chrV"), IRanges(start=1:99, width=1), "-")
    resSoll <- rep(0,99)
    pos <- Rle(start(aln[strand(aln)=="+"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 5: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(end(aln[strand(aln)=="-"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 6: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(end(aln[strand(aln)=="+"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 7: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(start(aln[strand(aln)=="-"]))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 8: qCount with orientation and query strand")
    
    query <- GRanges(c("chrV"), IRanges(start=1:99, width=1), "*")
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln), end(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 9: qCount with orientation and query strand")
    
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 10: qCount with orientation and query strand")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln), start(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="same")[,-1]
    checkTrue(all(resSoll == res), "Test 11: qCount with orientation and query strand")
    
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="opposite")[,-1]
    checkTrue(all(resSoll == res), "Test 12: qCount with orientation and query strand")
}

test_useRead <- function() {
    project <- projectPaired
    aln <- GenomicAlignments::readGAlignmentsFromBam(project@alignments$FileName, 
                                  param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=T, isSecondMateRead=F)))
    
    query <- GRanges(c("chrV"), IRanges(start=1:99, width=1), "+")
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln), end(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="any", useRead="first")[,-1]
    checkTrue(all(resSoll == res), "Test 1: qCount with useRead")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln), start(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="any", useRead="first")[,-1]
    checkTrue(all(resSoll == res), "Test 2: qCount with useRead")
    
    aln <- GenomicAlignments::readGAlignmentsFromBam(project@alignments$FileName, 
                                   param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=F, isSecondMateRead=T)))
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", start(aln), end(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="start", shift=0, orientation="any", useRead="last")[,-1]
    checkTrue(all(resSoll == res), "Test 3: qCount with useRead")
    
    resSoll <- rep(0,99)
    pos <- Rle(ifelse(strand(aln)=="+", end(aln), start(aln)))
    resSoll[runValue(pos)] <- runLength(pos)
    res <- qCount(project, query, selectReadPosition="end", shift=0, orientation="any", useRead="last")[,-1]
    checkTrue(all(resSoll == res), "Test 4: qCount with useRead")
}

test_maxInsertSize <- function() {
    project <- projectPaired

    query <- GRanges(c("chrV"), IRanges(start=1:20, width=1), "*")
    resSoll <- rep(0,20)
    #resSoll[c(10, 30, 46, 50, 70, 74, 86, 94)] <- c(4,4,2,2,2,2,2,2)
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", maxInsertSize=0)[,-1]
    checkTrue(all(resSoll == res), "Test 1: qCount with maxInsertSize")
    
    resSoll[10] <- 4
    res <- qCount(project, query, selectReadPosition="start", shift="halfInsert", maxInsertSize=14)[,-1]
    checkTrue(all(resSoll == res), "Test 2: qCount with maxInsertSize")
}

test_query_GRanges <- function() {
    project <- projectSingle
    
    ## NO masking
    ## reduce region by query rownames
    region <- tilingRegion
    strand(region) <- "*"    
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(300,300,300,250)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 1: qCount orientation=same")
    
    strand(region) <- "+"
    resSoll[,3] = 0
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 2: qCount orientation=same")
    
    strand(region) <- "-"
    resSoll[,2] = 0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 3: qCount orientation=same")
    
    ## NO reduce region by query rownames
    names(region) <- NULL
    strand(region) <- "*"
    resSoll <- matrix(0, nrow=12, ncol=3, byrow=T) 
    resSoll[,1] =  c(rep(100,11),50)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 4: qCount orientation=same")  
    
    strand(region) <- "+"
    resSoll[,3] =  0
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 5: qCount orientation=same")
    
    strand(region) <- "-"
    resSoll[,2] =  0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 6: qCount orientation=same")
    

    ## Masking Test 1
    mask <- tilingRegion[names(tilingRegion) == "H4"]
    
    ## reduce region by query rownames
    region <- tilingRegion
    strand(region) <- "+"    
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(200,300,150,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 7: qCount with masking and orientation=any")
    
    strand(region) <- "+"
    resSoll[,3] = 0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 8: qCount with masking and orientation=same")
    
    strand(region) <- "-"
    resSoll[,2] = 0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 9: qCount with masking and orientation=same")
    
    ## NO reduce region by query rownames
    names(region) <- NULL
    strand(region) <- "+"
    resSoll <- matrix(0, nrow=12, ncol=3, byrow=T) 
    resSoll[,1] =  c(100,100,50,0,50,100,50,0,50,100,50,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 10: qCount with masking and orientation=any")  
    
    strand(region) <- "+"
    resSoll[,3] =  0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 11: qCount with masking and orientation=same")
    
    strand(region) <- "-"
    resSoll[,2] =  0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 12: qCount with masking and orientation=same")
    
        
    ## Masking Test 2
    mask <- GRanges(seqnames="chrV", IRanges(c(361,401), c(390,700)))

    ## reduce region by query rownames
    region <- tilingRegion
    strand(region) <- "+"    
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(170,120,100,100)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 13: qCount with masking and orientation=any")
    
    resSoll[,3] = 0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 14: qCount with masking and orientation=same")
    
    ## NO reduce region by query rownames
    names(region) <- NULL
    strand(region) <- "+"
    resSoll <- matrix(0, nrow=12, ncol=3, byrow=T) 
    resSoll[,1] =  c(100,100,100,100,70,20,0,0,0,0,0,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 15: qCount with masking and orientation=any")  
    
    resSoll[,3] =  0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 16: qCount with masking and orientation=same")
}    
    
test_query_GRangesList <- function() {
    project <- projectSingle
    
    region <- tilingRegion
    strand(region) <- "*"
    regionList <- split(region, names(region))
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(300,150,150,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 1: qCount orientation=same")
    
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll[,3] = 0
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 2: qCount orientation=same")

    strand(region) <- "-"
    regionList <- split(region, names(region))
    resSoll[,2] = 0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 3: qCount orientation=same")

    ## Masking Test 1
    mask <- tilingRegion[names(tilingRegion) == "H4"]
    
    region <- tilingRegion
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(200,150,0,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]   
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRangesList Test 4: qCount with masking and orientation=any")
    
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll[,3] = 0
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 5: qCount with masking and orientation=same")
    
    strand(region) <- "-"
    regionList <- split(region, names(region))
    resSoll[,2] = 0
    resSoll[,3] = 3*resSoll[,1]
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 6: qCount with masking and orientation=same")
    
    ## Masking Test 2
    mask <- GRanges(seqnames="chrV", IRanges(c(361,401), c(390,700)))
    
    region <- tilingRegion
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll <- matrix(0, nrow=4, ncol=3, byrow=T) 
    resSoll[,1] = c(170,50,50,0)
    resSoll[,c(2,3)] = 3*resSoll[,1]   
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRangesList Test 7: qCount with masking and orientation=any")
    
    regionList <- split(region, names(region))
    resSoll[,3] = 0
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 8: qCount with masking and orientation=same")

}

test_query_GRanges_allelic <- function() {
    project <- projectSingleAllelic
    
    ## NO masking
    ## reduce region by query rownames
    region <- tilingRegion
    strand(region) <- "*"    
    resSoll <- matrix(0, nrow=4, ncol=7, byrow=T) 
    resSoll[,1] = c(300,300,300,250)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 1: qCount allele specific orientation=same")
    
    strand(region) <- "+"
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 2: qCount allele specific orientation=same")
    
    strand(region) <- "-"
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 3: qCount allele specific orientation=same")
    
    ## NO reduce region by query rownames
    names(region) <- NULL
    strand(region) <- "*"
    resSoll <- matrix(0, nrow=12, ncol=7, byrow=T) 
    resSoll[,1] =  c(rep(100,11),50)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 4: qCount allele specific orientation=same")  
    
    strand(region) <- "+"
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 5: qCount allele specific orientation=same")
    
    strand(region) <- "-"
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, region, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 6: qCount allele specific orientation=same")
    
    ## Masking Test 1
    mask <- tilingRegion[names(tilingRegion) == "H4"]
  
    ## reduce region by query rownames
    region <- tilingRegion
    strand(region) <- "+"    
    resSoll <- matrix(0, nrow=4, ncol=7, byrow=T) 
    resSoll[,1] = c(200,300,150,0)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 7: qCount allele specific with masking and orientation=same")
    
    strand(region) <- "+"
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 8: qCount allele specific with masking and orientation=same")
    
    strand(region) <- "-"
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 9: qCount allele specific with masking and orientation=same")
    
    ## NO reduce region by query rownames
    names(region) <- NULL
    strand(region) <- "+"
    resSoll <- matrix(0, nrow=12, ncol=7, byrow=T) 
    resSoll[,1] =  c(100,100,50,0,50,100,50,0,50,100,50,0)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRanges Test 10: qCount allele specific with masking and orientation=any")  
    
    strand(region) <- "+"
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 11: qCount allele specific with masking and orientation=same")
    
    strand(region) <- "-"
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, region, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRanges Test 12: qCount allele specific with masking and orientation=same")
}

test_query_GRangesList_allelic <- function() {
    project <- projectSingleAllelic

    region <- tilingRegion
    strand(region) <- "*"
    regionList <- split(region, names(region))
    resSoll <- matrix(0, nrow=4, ncol=7, byrow=T) 
    resSoll[,1] = c(300,150,150,0)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 1: qCount allele specific orientation=same")
    
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 2: qCount allele specific orientation=same")
    
    strand(region) <- "-"
    regionList <- split(region, names(region))
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, regionList, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 3: qCount allele specific orientation=same")
  
    ## Masking Test 1
    mask <- tilingRegion[names(tilingRegion) == "H4"]
    
    region <- tilingRegion
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll <- matrix(0, nrow=4, ncol=7, byrow=T) 
    resSoll[,1] = c(200,150,0,0)
    resSoll[,c(2,3,4,5,6,7)] = resSoll[,1]
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="any")
    checkTrue(all(resSoll == res), "GRangesList Test 4: qCount allele specific with masking and orientation=any")
    
    strand(region) <- "+"
    regionList <- split(region, names(region))
    resSoll[,c(5,6,7)] = 0
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 5: qCount allele specific with masking and orientation=same")
    
    strand(region) <- "-"
    regionList <- split(region, names(region))
    resSoll[,c(2,3,4)] = 0
    resSoll[,c(5,6,7)] = resSoll[,1]
    res <- qCount(project, regionList, mask=mask, collapseBySample=F, orientation="same")
    checkTrue(all(resSoll == res), "GRangesList Test 6: qCount allele specific with masking and orientation=same")
}

test_query_transcriptDB <- function() {
    project <- qAlign(sampleFile, genomeFile, splicedAlignment=TRUE, alignmentsDir=td, clObj=clObj)
    txdb <- createTranscriptDb()
    mask <- gtfRegion[mcols(gtfRegion)$gene_name == "TNFRSF18"]
    
    
    ## TranscriptDB vs GRanges
    # Gene
    res <- qCount(project, txdb, collapseBySample=F, reportLevel=NULL)
    resTxdb <- qCount(project, txdb, collapseBySample=F, reportLevel="gene")
    checkTrue(all(resTxdb == res), "TranscriptDB vs GRanges Test 1")
    
    region <- gtfRegion
    names(region) <- mcols(gtfRegion)$gene_id
    resGr <- qCount(project, region, collapseBySample=F)
    resGr <- resGr[sort(rownames(resGr)),]
    checkTrue(all(resTxdb == resGr), "TranscriptDB vs GRanges Test 2")
    
    # Exon
    resTxdb <- qCount(project, txdb, collapseBySample=F, reportLevel="exon")
    
    region <- exons(txdb)
    resGr <- qCount(project, region, collapseBySample=F)
    resGr <- resGr[sort(rownames(resGr)),]
    checkTrue(all(resTxdb == resGr), "TranscriptDB vs GRanges Test 3")    
    
    # Promoter
    resTxdb <- qCount(project, txdb, collapseBySample=F, reportLevel="promoter")
    
    region <- promoters(txdb, columns=c("tx_id","tx_name"))
    names(region) <- paste(mcols(region)$tx_id,mcols(region)$tx_name, sep=";")
    resGr <- qCount(project, region, collapseBySample=F)
    resGr <- resGr[sort(rownames(resGr)),]
    checkTrue(all(resTxdb == resGr), "TranscriptDB vs GRanges Test 4")
    
    # junction, includeSpliced
    exGr <- GRanges(c("chr1","chr1","chr1","chr1"),
                    IRanges(start=c(11720,12322,14043,14363), end=c(12212,12518,14165,14512)))
    resE <- qCount(project, exGr, collapseBySample=FALSE)
    resEU <- qCount(project, exGr, collapseBySample=FALSE, includeSpliced=FALSE)
    inGr <- GRanges(c("chr1","chr1"), IRanges(start=c(12213,14166), end=c(12321,14362)), strand=c("+","+"))
    resJ <- qCount(project, NULL, reportLevel="junction", collapseBySample=FALSE)
    checkTrue(all((resE - resEU)[c(1,3),-1] == as.matrix(mcols(resJ[match(inGr, resJ)]))), "junction/includeSpliced Test 1")
    
    ## TranscriptDB vs GRanges with masked region
    resTxdb <- qCount(project, txdb, collapseBySample=F, mask=mask, reportLevel="gene")
    
    region <- gtfRegion
    names(region) <- mcols(gtfRegion)$gene_id
    resGr <- qCount(project, region, collapseBySample=F, mask=mask)
    resGr <- resGr[sort(rownames(resGr)),]
    checkTrue(all(resTxdb == resGr), "TranscriptDB vs GRanges Test 5")

    ## Collapse by sample
    resTxdbCS <- qCount(project, txdb, collapseBySample=T, mask=mask, reportLevel="gene")
    res <- cbind(rowSums(resTxdb[,c(2,3)]), rowSums(resTxdb[,c(4,5)]))
    checkTrue(all(res == resTxdbCS[,c(2,3)]), "TranscriptDB collapse by Sample Test")    
    
}

test_collapseBySample_GRanges <- function() {
    ## Non Allelic
    project <- qAlign(sampleFile, genomeFile, alignmentsDir=td, clObj=clObj)
    res <- qCount(project, gtfRegion, collapseBySample=T)
    
    projectS1 <- project[1:2]
    resS1 <- qCount(projectS1, gtfRegion, collapseBySample=T)
    
    projectS2 <- project[3:4]
    resS2 <- qCount(projectS2, gtfRegion, collapseBySample=T)

    checkTrue(all(res[,2:3] == cbind(resS1[,2], resS2[,2])),
              "Test collapseBySample; Collapse counts are not equal.")

    ## Allelic
    project <- qAlign(sampleFile, genomeFile, snpFile=snpFile, alignmentsDir=td, clObj=clObj)
    res <- qCount(project, gtfRegion, collapseBySample=T)
    
    projectS1 <- project[1:2]
    resS1 <- qCount(projectS1, gtfRegion, collapseBySample=T)
    
    projectS2 <- project[3:4]
    resS2 <- qCount(projectS2, gtfRegion, collapseBySample=T)
    
    checkTrue(all(res[,2:7] == cbind(resS1[,2:4], resS2[,2:4])),
              "Test collapseBySample; Collapse counts are not equal for allelic.")
}

test_auxiliaryName <- function() {
    project <- qAlign(file.path("extdata", "samples_chip_single.txt"), genomeFile,
                      auxiliaryFile=file.path("extdata", "auxiliaries.txt"), alignmentsDir=td, clObj=clObj)
    
    ## Test aux counts
    auxRegion <- createAuxRegion()
    res <- qCount(project, auxRegion, collapseBySample=F, auxiliaryName="phiX174")
    resSoll <- c(251,493)
    checkTrue(all(resSoll == res[,2:3]))
}
