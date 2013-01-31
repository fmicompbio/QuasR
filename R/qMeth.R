# proj      : qProject object
# query     : NULL (whole genome)
#             GRanges object (only quantify C's in these regions)
# mode      : "allC"            : all C's (+/- strands separate)
#             "CpG"             : only C's in CpG context (+/- strands separate)
#             "CpGcomb"(default): only C's in CpG context (+/- strands collapsed)
#             "var"             : variant detection (all C's, +/- strands separate)
# collapseBySample : combine (sum) counts from bamfiles with the same sample name
# collapseByQueryRegion : combine (sum) counts for C's per query region
# asGRanges : return value as GRanges object or data.frame
# mask      : mask genomic regions (e.g. unmappable regions)
# reference : source of bam files ("genome" or AuxName)
# keepZero  : return C's with total==0 in results
# clObj     : cluster object for parallelization
#
# value     : if asGRanges==TRUE, GRanges object with one region per quantified C and two metadata columns per sample (_T, _M)
#             else, data.frame with one row per quantified C and in the columns the coordinates of the C, as well as two count (_T, _M) per sample

# TODO:
# - support for gapped aligments (parsing of cigar strings at C level)
# - ignore presumable PCR duplicates (multiple alignment pairs with same external coordinates)
# - ignore overlapping read pairs (insert size smaller than twice the read length)

qMeth <-
    function(proj,
             query=NULL,
             mode=c("CpGcomb","CpG","allC","var"),
             collapseBySample=TRUE,
             collapseByQueryRegion=FALSE,
             asGRanges=TRUE,
             mask=NULL,
             reference="genome",
             keepZero=TRUE,
             clObj=NULL) {
        ## setup variables from 'proj' -----------------------------------------------------------------------
        # 'proj' is correct type?
        if(!inherits(proj, "qProject", which=FALSE))
            stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")
        if(proj@bisulfite=="no")
            stop("'proj' is not a bisufite-seq project")
        if(proj@splicedAlignment)
            stop("'spliceAlignment==TRUE' is not supported by qMeth")

        samples <- proj@alignments$SampleName
        nsamples <- length(samples)
        if(reference=="genome") {
            bamfiles <- proj@alignments$FileName
            referenceFormat <- proj@genomeFormat
            referenceSource <- proj@genome
        } else if(!is.na(i <- match(reference, rownames(proj@auxAlignments)))) {
            bamfiles <- unlist(proj@auxAlignments[i, ], use.names=FALSE)
            referenceFormat <- "file"
            referenceSource <- proj@aux[i,'FileName']
        } else {
            stop("unknown 'reference', should be one of: ",
                 paste(sprintf("'%s'", c("genome",rownames(proj@auxAlignments))), collapse=", "))
        }
        mode <- match.arg(mode)

        
        ## validate parameters -------------------------------------------------------------------------------
        # 'query' is correct type?
        if(is.null(query)) {
            tr <- scanBamHeader(bamfiles[1])[[1]]$targets
            query <- GRanges(names(tr), IRanges(start=1, end=tr))
        } else if(!inherits(query,"GRanges")) {
            stop("'query' must be either NULL or an object of type 'GRanges'")
        }

        if(!is.logical(keepZero) || length(keepZero)!=1)
            stop("'keepZero' must be either TRUE or FALSE")
        if(!is.logical(asGRanges) || length(asGRanges)!=1)
            stop("'asGRanges' must be either TRUE or FALSE")

        if(!keepZero && ((collapseBySample && length(unique(samples))>1) || (!collapseBySample && length(samples)>1)))
            stop("'keepZero' must be TRUE if there are multiple non-collapsable samples")

        if(mode=="var" && collapseByQueryRegion)
            stop("'collapseByQueryRegion' must be FALSE for variant detection mode")
        if(mode=="var" && !is.na(proj@snpFile))
            stop("allele-specific mode cannot be combined with variant detection mode")


        # all query chromosomes present in all bamfiles?
        trTab <- table(unlist(lapply(scanBamHeader(bamfiles), function(bh) names(bh$targets))))
        trCommon <- names(trTab)[trTab==length(bamfiles)]
        if(any(f <- !(seqlevels(query) %in% trCommon)))
            stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                         paste(seqlevels(query)[f],collapse=", ")))


        ## apply 'mask' to query -----------------------------------------------------------------------------
        if(!is.null(mask)) {
            stop("'mask' (masking of query regions, e.g. unmappable genomic regions) is not implemented yet")
        }
        
        
        ## setup tasks for parallelization -------------------------------------------------------------------
        ## TODO: create several chunks per chromosome?
        taskIByQuery <- split(seq.int(length(query)), as.factor(seqnames(query)))
        taskIByQuery <- taskIByQuery[ sapply(taskIByQuery,length)>0 ]
        nChunkQuery <- length(taskIByQuery)

        if(collapseBySample) {
            taskBamfiles <- split(bamfiles, samples)
            sampleNames <- names(taskBamfiles)
        } else {
            taskBamfiles <- as.list(bamfiles)
            sampleNames <- displayNames(proj)
        }
        nChunkBamfile <- length(taskBamfiles)

        nChunk <- nChunkQuery * nChunkBamfile
        taskIByQuery <- rep(taskIByQuery, nChunkBamfile)
        taskBamfiles <- rep(taskBamfiles, each=nChunkQuery)

        if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
            message("preparing to run on ", min(nChunk,length(clObj)), " nodes...", appendLF=FALSE)
            ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
            if(!all(sapply(ret, function(x) "QuasR" %in% x)))
                stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
            myapply <- function(...)
                parallel::clusterApplyLB(clObj, ...)
            message("done")
        } else {
            myapply <- function(...)
                lapply(...)
        }


        ## quantify methylation  -----------------------------------------------------------------------------
        if(!is.na(proj@snpFile)) {
            # allele-bis-mode: 6 columns per sample in output (TR, MR, TU, MU, TA, MA) -------------
            resL <- myapply(seq_len(nChunk),
                            function(i) quantifyMethylationBamfilesRegionsSingleChromosomeAllele(taskBamfiles[[i]],
                                                                                                 query[taskIByQuery[[i]]],
                                                                                                 collapseByQueryRegion,
                                                                                                 mode,
                                                                                                 referenceFormat,
                                                                                                 referenceSource,
                                                                                                 proj@snpFile,
                                                                                                 keepZero))
            # combine and reorder chunks
            # ...rbind chunks for the same sample
            resL <- lapply(split(seq_len(nChunk),rep(seq_len(nChunkBamfile),each=nChunkQuery)),
                           function(i) do.call(rbind, resL[i]))
            names(resL) <- sampleNames

            if(any(unlist(lapply(resL,nrow),use.names=FALSE)!=nrow(resL[[1]])))
                stop("error while combining partial results (chunks are incompatable)")
                
            # ...cbind TR/MR/TU/MU/TA/MA columns for different samples
            res <- cbind(resL[[1]][,c("chr","start","end","strand")],
                         do.call(cbind, lapply(resL, "[", c("TR","MR","TU","MU","TA","MA"))),
                         stringsAsFactors=FALSE)
            colnames(res)[5:ncol(res)] <- sprintf("%s_%s",rep(sampleNames,each=6),c("TR","MR","TU","MU","TA","MA"))


        } else if(mode == "var") {
            # variant detection mode: 2 columns per sample in output (total and match) -------------
            resL <- myapply(seq_len(nChunk),
                            function(i) detectVariantsBamfilesRegionsSingleChromosome(taskBamfiles[[i]],
                                                                                      query[taskIByQuery[[i]]],
                                                                                      referenceFormat,
                                                                                      referenceSource,
                                                                                      keepZero))
            # combine and reorder chunks
            # ...rbind chunks for the same sample
            resL <- lapply(split(seq_len(nChunk),rep(seq_len(nChunkBamfile),each=nChunkQuery)),
                           function(i) do.call(rbind, resL[i]))
            names(resL) <- sampleNames

            if(any(unlist(lapply(resL,nrow),use.names=FALSE)!=nrow(resL[[1]])))
                stop("error while combining partial results (chunks are incompatable)")
                
            # ...cbind T/M columns for different samples
            res <- cbind(resL[[1]][,c("chr","start","end","strand")],
                         do.call(cbind, lapply(resL, "[", c("T","M"))),
                         stringsAsFactors=FALSE)
            colnames(res)[5:ncol(res)] <- sprintf("%s_%s",rep(sampleNames,each=2),c("T","M"))

        } else {
            # normal mode: 2 columns per sample in output (T and M) -------------
            resL <- myapply(seq_len(nChunk),
                            function(i) quantifyMethylationBamfilesRegionsSingleChromosome(taskBamfiles[[i]],
                                                                                           query[taskIByQuery[[i]]],
                                                                                           collapseByQueryRegion,
                                                                                           mode,
                                                                                           referenceFormat,
                                                                                           referenceSource,
                                                                                           keepZero))
            # combine and reorder chunks
            # ...rbind chunks for the same sample
            resL <- lapply(split(seq_len(nChunk),rep(seq_len(nChunkBamfile),each=nChunkQuery)),
                           function(i) do.call(rbind, resL[i]))
            names(resL) <- sampleNames

            if(any(unlist(lapply(resL,nrow),use.names=FALSE)!=nrow(resL[[1]])))
                stop("error while combining partial results (chunks are incompatable)")
                
            # ...cbind T/M columns for different samples
            res <- cbind(resL[[1]][,c("chr","start","end","strand")],
                         do.call(cbind, lapply(resL, "[", c("T","M"))),
                         stringsAsFactors=FALSE)
            colnames(res)[5:ncol(res)] <- sprintf("%s_%s",rep(sampleNames,each=2),c("T","M"))
        }

        if(asGRanges) {
            if(referenceFormat=="file") {
                si <- seqinfo(scanFaIndex(referenceSource))
            } else {
                library(referenceSource, character.only=TRUE)
                gnmObj <- get(ls(sprintf("package:%s", referenceSource)))
                si <- seqinfo(gnmObj)
            }
            res <- GRanges(seqnames=res$chr, IRanges(start=res$start, end=res$end),
                           strand=res$strand, seqinfo=si, res[,5:ncol(res)])
        }

        ## return results
        return(res)
    }


# detect variants for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, never collapsed)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+2*nSamples vectors: chr, start, end, strand of C, counts of T (total) and M (match) reads
detectVariantsBamfilesRegionsSingleChromosome <-
    function(bamfiles, regions, referenceFormat, reference, keepZero) {
        ## verify parameters
        if(length(chr <- as.character(unique(seqnames(regions)))) != 1)
            stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")

        
        ## collapse regions
        regionsStart <- as.integer(min(start(regions)))
        regionsEnd   <- as.integer(max(end(regions)))
        regionsGr    <- GRanges(chr, IRanges(start=regionsStart, end=regionsEnd))

        
        ## get sequence string from...
        #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
        if(referenceFormat=="file") { # genome file
            chrLen <- as.integer(seqlengths(scanFaIndex(reference))[chr])
            seqstr <- as.character(scanFa(reference, regionsGr)[[1]])

        } else {                      # BSgenome object
            library(reference, character.only=TRUE)
            referenceObj <- get(ls(sprintf("package:%s", reference))) # access the BSgenome
            chrLen <- as.integer(length(referenceObj[[chr]]))
            seqstr <- getSeq(referenceObj, regionsGr, as.character=TRUE)
        }
        #message("done")


        ## call CPP function (multiple bam files, single region)
        #message("detecting single nucleotide variations...", appendLF=FALSE)
        resL <- .Call("detectSNVs", bamfiles, chr, chrLen, regionsStart, seqstr, keepZero, PACKAGE="QuasR")
        #message("done")


        ## filter out C's that do not fall into 'regions'
        #message("processing results...", appendLF=FALSE)
        ov <- findOverlaps(query=GRanges(chr, IRanges(start=resL$position,width=1)), subject=regions)
        resL <- lapply(resL, "[", unique(queryHits(ov)))

        res <- data.frame(chr=resL$chr,
                          start=resL$position,
                          end=resL$position,
                          strand=rep("*",length(resL$chr)),
                          T=resL$nTotal,
                          M=resL$nMatch,
                          stringsAsFactors=FALSE)
        #message("done")

        res
    }


# quantify methylation for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, may be collapsed if collapseByRegion==TRUE)
#  - mode (defines which and how C's are quantified)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+2*nSamples vectors: chr, start, end, strand of C, counts of T (total) and M (methylated) reads
quantifyMethylationBamfilesRegionsSingleChromosome <-
    function(bamfiles, regions, collapseByRegion, mode=c("CpGcomb","CpG","allC"), referenceFormat, reference, keepZero) {
        ## verify parameters
        if(length(chr <- as.character(unique(seqnames(regions)))) != 1)
            stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")
        mode <- c("CpGcomb"=0L,"CpG"=1L,"allC"=2L)[match.arg(mode)]
        Cwidth <- ifelse(mode==0, 2L, 1L)

        
        ## collapse regions
        regionsStart <- as.integer(min(start(regions)))
        regionsEnd   <- as.integer(max(end(regions)))
        regionsGr    <- GRanges(chr, IRanges(start=regionsStart, end=regionsEnd))

        
        ## get sequence string from...
        #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
        if(referenceFormat=="file") { # genome file
            chrLen <- as.integer(seqlengths(scanFaIndex(reference))[chr])
            seqstr <- as.character(scanFa(reference, regionsGr)[[1]])

        } else {                      # BSgenome object
            library(reference, character.only=TRUE)
            referenceObj <- get(ls(sprintf("package:%s", reference))) # access the BSgenome
            chrLen <- as.integer(length(referenceObj[[chr]]))
            seqstr <- getSeq(referenceObj, regionsGr, as.character=TRUE)
        }
        #message("done")


        ## call CPP function (multiple bam files, single region)
        #message("quantifying methylation...", appendLF=FALSE)
        resL <- .Call("quantifyMethylation", bamfiles, chr, chrLen, regionsStart, seqstr, mode, keepZero, PACKAGE="QuasR")
        #message("done")


        ## collapse by region
        #message("processing results...", appendLF=FALSE)
        ov <- findOverlaps(query=GRanges(chr, IRanges(start=resL$position,width=Cwidth)), subject=regions)

        if(collapseByRegion) {
            tmpT <- tmpM <- numeric(length(regions))
            tmp <- tapply(resL$T[queryHits(ov)], subjectHits(ov), sum)
            tmpT[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$M[queryHits(ov)], subjectHits(ov), sum)
            tmpM[as.integer(names(tmp))] <- tmp
            res <- data.frame(chr=rep(chr,length(regions)),
                              start=start(regions),
                              end=end(regions),
                              strand=as.character(strand(regions)),
                              T=tmpT,
                              M=tmpM,
                              stringsAsFactors=FALSE)

        } else {
            ## filter out C's that do not fall into 'regions'
            resL <- lapply(resL, "[", unique(queryHits(ov)))

            res <- data.frame(chr=resL$chr,
                              start=resL$position,
                              end=resL$position+Cwidth-1L,
                              strand=resL$strand,
                              T=resL$T,
                              M=resL$M,
                              stringsAsFactors=FALSE)
        }
        #message("done")

        res
    }



# quantify ALLELE-SPECIFIC methylation for:
#  - multiple bamfiles (will allways be collapsed)
#  - multiple regions (all on single chromosome, may be collapsed if collapseByRegion==TRUE)
#  - mode (defines which and how C's are quantified)
#  - referenceFormat and reference (access to sequence at 'regions')
# return a data.frame or GRanges object with 4+6*nSamples vectors: chr, start, end, strand of C, counts of TR, TU, TA (total) and MR, MU, MA (methylated) reads
quantifyMethylationBamfilesRegionsSingleChromosomeAllele <-
    function(bamfiles, regions, collapseByRegion, mode=c("CpGcomb","CpG","allC"), referenceFormat, reference, snpFile, keepZero) {
        ## verify parameters
        if(length(chr <- as.character(unique(seqnames(regions)))) != 1)
            stop("all regions need to be on the same chromosome for 'quantifyMethylationBamfilesRegionsSingleChromosome'")
        mode <- c("CpGcomb"=0L,"CpG"=1L,"allC"=2L)[match.arg(mode)]
        Cwidth <- ifelse(mode==0, 2L, 1L)

        
        ## collapse regions
        regionsStart <- as.integer(min(start(regions)))
        regionsEnd   <- as.integer(max(end(regions)))
        regionsGr    <- GRanges(chr, IRanges(start=regionsStart, end=regionsEnd))

        
        ## get sequence string from...
        #message("loading reference sequence (", chr, ")...", appendLF=FALSE)
        if(referenceFormat=="file") { # genome file
            chrLen <- as.integer(seqlengths(scanFaIndex(reference))[chr])
            seqstr <- as.character(scanFa(reference, regionsGr)[[1]])

        } else {                      # BSgenome object
            library(reference, character.only=TRUE)
            referenceObj <- get(ls(sprintf("package:%s", reference))) # access the BSgenome
            chrLen <- as.integer(length(referenceObj[[chr]]))
            seqstr <- getSeq(referenceObj, regionsGr, as.character=TRUE)
        }
        #message("done")


        ## call CPP function (multiple bam files, single region)
        #message("quantifying methylation...", appendLF=FALSE)
        resL <- .Call("quantifyMethylationAllele", bamfiles, chr, chrLen, regionsStart, seqstr, mode, keepZero, PACKAGE="QuasR")
        #message("done")


        ## filter out CpGs that overlap SNPs (may not be possible to discriminate allele from methylation status)
        #message("removing C's overlapping SNPs...", appendLF=FALSE)
        snpL <- scan(snpFile, what=list(chr="", pos=1L, R="", A=""), quiet=TRUE)
        snp <- GRanges(snpL$chr, IRanges(start=snpL$pos, width=nchar(snpL$R)))
        ikeep <- !overlapsAny(GRanges(chr, IRanges(start=resL$position,width=Cwidth)), snp)
        resL <- lapply(resL, "[", ikeep)
        #message(sprintf("removed %d C's, done",sum(!ikeep)))
        

        ## collapse by region
        #message("processing results...", appendLF=FALSE)
        ov <- findOverlaps(query=GRanges(chr, IRanges(start=resL$position,width=Cwidth)), subject=regions)

        if(collapseByRegion) {
            tmpTR <- tmpTU <- tmpTA <- tmpMR <- tmpMU <- tmpMA <- numeric(length(regions))
            tmp <- tapply(resL$TR[queryHits(ov)], subjectHits(ov), sum)
            tmpTR[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$TU[queryHits(ov)], subjectHits(ov), sum)
            tmpTU[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$TA[queryHits(ov)], subjectHits(ov), sum)
            tmpTA[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MR[queryHits(ov)], subjectHits(ov), sum)
            tmpMR[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MU[queryHits(ov)], subjectHits(ov), sum)
            tmpMU[as.integer(names(tmp))] <- tmp
            tmp <- tapply(resL$MA[queryHits(ov)], subjectHits(ov), sum)
            tmpMA[as.integer(names(tmp))] <- tmp
            res <- data.frame(chr=rep(chr,length(regions)),
                              start=start(regions),
                              end=end(regions),
                              strand=as.character(strand(regions)),
                              TR=tmpTR,
                              MR=tmpMR,
                              TU=tmpTU,
                              MU=tmpMU,
                              TA=tmpTA,
                              MA=tmpMA,
                              stringsAsFactors=FALSE)

        } else {
            ## filter out C's that do not fall into 'regions'
            resL <- lapply(resL, "[", unique(queryHits(ov)))

            res <- data.frame(chr=resL$chr,
                              start=resL$position,
                              end=resL$position+Cwidth-1L,
                              strand=resL$strand,
                              TR=resL$TR,
                              MR=resL$MR,
                              TU=resL$TU,
                              MU=resL$MU,
                              TA=resL$TA,
                              MA=resL$MA,
                              stringsAsFactors=FALSE)
        }
        #message("done")

        res
    }

