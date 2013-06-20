## query : GRanges : value has one row per unique name in GRanges
## query : GRangesList : value has one row per name of GRangesList
## query : TranscriptDB : value has has one row per 'reportLevel'
## query : NULL & reportLevel=="junction": value has has one row per 'junction'

qCount <-
    function(proj,
             query,
             reportLevel=c(NULL,"gene","exon","promoter","junction"),
             selectReadPosition=c("start", "end"),
             shift=0L,
             orientation=c("any","same","opposite"),
             useRead=c("any","first","last"),
             auxiliaryName=NULL,
             mask=NULL,
             collapseBySample=TRUE,
             includeSpliced=TRUE,
             maxInsertSize=500L,
             clObj=NULL) {
        ## setup variables from 'proj' -----------------------------------------------------------------------
        ## 'proj' is correct type?
        if(!inherits(proj, "qProject", which=FALSE))
            stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")

        samples <- proj@alignments$SampleName
        nsamples <- length(samples)
        bamfiles <-
            if(is.null(auxiliaryName))
                proj@alignments$FileName
            else if(!is.na(i <- match(auxiliaryName, proj@aux$AuxName)))
                unlist(proj@auxAlignments[i,], use.names=FALSE)
            else
                stop("unknown 'auxiliaryName', should be one of: NULL, ",
                     paste(sprintf("'%s'", proj@aux$AuxName), collapse=", "))

        
        ## validate parameters -------------------------------------------------------------------------------
        reportLevel <- match.arg(reportLevel)
        selectReadPosition <- match.arg(selectReadPosition)
        orientation <- match.arg(orientation)
        useRead <- match.arg(useRead)

        ## check shift
        if(length(shift) == 1 && shift == "halfInsert") {
            if(proj@paired == "no") {
                stop("'shift=\"halfInsert\"' can only be used for paired-end experiments")
            } else {
                shifts <- rep(-1000000L, nsamples)
                broaden <- as.integer(ceiling(maxInsertSize/2))
            }
        } else {
            if(!is.numeric(shift) || (length(shift)>1 && length(shift)!=nsamples))
                stop(sprintf("'shift' must be 'halfInsert', a single integer or an integer vector with %d values",nsamples))
            else if(length(shift) == 1)
                shifts <- rep(as.integer(shift), nsamples)
            else
                shifts <- as.integer(shift)
            broaden <- 0L
        }

        ## all query chromosomes present in all bamfiles?
        trTab <- table(unlist(lapply(scanBamHeader(bamfiles), function(bh) names(bh$targets))))
        trCommon <- names(trTab)[trTab==length(bamfiles)]
        if(!is.null(query) && any(f <- !(seqlevels(query) %in% trCommon)))
            stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                         paste(seqlevels(query)[f],collapse=", ")))

        ## 'query' is correct type?
        if(reportLevel == "junction") {
            if(proj@splicedAlignment != TRUE && proj@samplesFormat != "bam")
                stop("reportLevel=\"junction\" cannot be used for non-spliced alignments")
            if(!is.null(query))
                warning("ignoring 'query' for reportLevel=\"junction\"")

            ### reportLevel == "junction" -------------------------------------------------------------------------------
            ## setup tasks for parallelization -------------------------------------------------------------------
            if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
                message("preparing to run on ", length(clObj), " nodes...", appendLF=FALSE)
                ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
                if(!all(sapply(ret, function(x) "QuasR" %in% x)))
                    stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
                taskTargets <- rep(trCommon, nsamples)
                bamfiles <- rep(bamfiles, each=length(trCommon))
                myapply <- function(...) {
                    ret <- clusterMap(clObj, ..., SIMPLIFY=FALSE, .scheduling="dynamic")
                    # fuse
                    if(!is.na(proj@snpFile)){ # ret is a list of list(id,R,U,A)
                        ret <- lapply(split(seq_along(bamfiles),bamfiles), function(i) list(id=do.call(c, lapply(ret[i], "[[", "id")),
                                                                                            R=do.call(c, lapply(ret[i], "[[", "R")),
                                                                                            U=do.call(c, lapply(ret[i], "[[", "U")),
                                                                                            A=do.call(c, lapply(ret[i], "[[", "A"))))
                    } else {                  # ret is a list of named vectors
                        ret <- lapply(split(seq_along(bamfiles),bamfiles), function(i) do.call(c, unname(ret[i])))
                    }
                    ret
                }
                message("done")
            } else {
                taskTargets <- rep(list(NULL),length(bamfiles))
                myapply <- function(...) {
                    ret <- mapply(..., SIMPLIFY=FALSE)
                    ret
                }
            }

            ## count junctions ----------------------------------------------------------------------------------
            message("counting junctions...", appendLF=FALSE)
            resL <- myapply(countJunctionsOneBamfile,
                            bamfile=bamfiles,
                            targets=taskTargets,
                            MoreArgs=list(allelic=!is.na(proj@snpFile)))
            message("done")

            ## make result rectangular and collapse (sum) counts by sample if necessary
            if(!is.na(proj@snpFile)){ # ret is a list of list(id,R,U,A)
                allJunctions <- unique(Reduce(c, lapply(resL, "[[", "id")))
                res <- matrix(0, nrow=length(allJunctions), ncol=3*length(resL), dimnames=list(allJunctions, NULL))
                for(i in 1:length(resL)) # make res a matrix with 3 columns per sample
                    res[resL[[i]]$id,((i-1)*3+1):((i-1)*3+3)] <- do.call(cbind, resL[[i]][c("R","U","A")])
                if(nsamples > length(unique(samples))) {
                    if(collapseBySample) {
                        message("collapsing counts by sample...", appendLF=FALSE)
                        iBySample <- split(seq_len(nsamples),samples)[unique(samples)]
                        res <- do.call(cbind, lapply(iBySample, function(i)
                                                     cbind(R=rowSums(res[,(i-1)*3+1,drop=FALSE]),
                                                           U=rowSums(res[,(i-1)*3+2,drop=FALSE]),
                                                           A=rowSums(res[,(i-1)*3+3,drop=FALSE]))))
                        colnames(res) <- paste(rep(names(iBySample),each=3), c("R","U","A"), sep="_")
                        message("done")
                        
                    } else {
                        # unify non-collapsed identical sample names
                        colnames(res) <- paste(rep(displayNames(proj),each=3), c("R","U","A"), sep="_")
                    }
                } else {
                    colnames(res) <- paste(rep(samples, each=3), c("R","U","A"), sep="_")
                }
                
            } else {                  # ret is a list of named vectors
                allJunctions <- unique(Reduce(c, lapply(resL, names)))
                res <- matrix(0, nrow=length(allJunctions), ncol=length(resL), dimnames=list(allJunctions, NULL))
                for(i in 1:length(resL))
                    res[names(resL[[i]]),i] <- resL[[i]]
                if(nsamples > length(unique(samples))) {
                    if(collapseBySample) {
                        message("collapsing counts by sample...", appendLF=FALSE)
                        iBySample <- split(seq_len(nsamples),samples)[unique(samples)]
                        res <- do.call(cbind, lapply(iBySample, function(i) rowSums(res[,i,drop=FALSE])))
                        message("done")

                    } else {
                        # unify non-collapsed identical sample names
                        colnames(res) <- displayNames(proj)
                    }
                } else {
                    colnames(res) <- samples
                }
            }

            ## make GRanges object
            res2 <- GRanges(seqnames=sub("^(.+):[0-9]+:[0-9]+:.$","\\1",allJunctions),
                            ranges=IRanges(start=as.integer(sub("^.+:([0-9]+):[0-9]+:.$","\\1",allJunctions)),
                            end=as.integer(sub("^.+:[0-9]+:([0-9]+):.$","\\1",allJunctions))),
                            strand=sub("^.+:[0-9]+:[0-9]+:(.)$","\\1",allJunctions))
            mcols(res2) <- res
            
            ## return results
            return(res2)

           
        } else {
            ### reportLevel != "junction" -------------------------------------------------------------------------------
            if(!inherits(query,c("GRanges","GRangesList","TranscriptDb")))
                stop("'query' must be either an object of type 'GRanges', 'GRangesList' or 'TranscriptDb', or NULL for reportLevel=\"junction\"")
                
            ## 'useRead' set but not a paired-end experiment?
            if(useRead != "any" && proj@paired == "no")
                warning("ignoring 'useRead' for single read experiments")


            ## preprocess query -------------------------------------------------------------------------------
            ##    --> create 'flatquery', 'querynames', 'querylengths' and 'zeroquerynames'
            ##    GRanges query -------------------------------------------------------------------------------
            if(inherits(query,"GRanges")) {
                if(!is.null(names(query)) && length(query) > length(unique(names(query)))) {
                    # remove redundancy from 'query' by names
                    tmpquery <- reduce(split(query, names(query))[unique(names(query))])
                    flatquery <- unlist(tmpquery, use.names=FALSE)
                    querynames <- rep(names(tmpquery), elementLengths(tmpquery))
                    rm(tmpquery)
                } else {
                    flatquery <- query
                    querynames <- if(is.null(names(query))) as.character(seq_len(length(query))) else names(query)
                }
                querylengths <- width(flatquery)
                zeroquerynames <- character(0)
            
                ##    GRangesList query ---------------------------------------------------------------------------
            } else if(inherits(query,"GRangesList")) {
                if(any(i <- elementLengths(query)==0)) {
                    warning(sprintf("removing %d elements from 'query' with zero regions: %s",sum(i),paste(names(query)[i],collapse=", ")))
                    query <- query[-i]
                }
                # hierarchically remove redundancy from 'query'
                message("hierarchically removing redundancy from 'query'...", appendLF=FALSE)
                if(orientation=="any")
                    strand(query) <- endoapply(strand(query), function(x) Rle(factor("*", levels=c("+","-","*")), lengths=length(x)))
                query <- reduce(query)
                if(length(query)>1) {
                    cumquery <- query[[1]]
                    for(i in 2:length(query)) {
                        query[[i]] <- setdiff(query[[i]], cumquery)
                        cumquery <- c(query[[i]], cumquery)
                    }
                }
                message("done")
                flatquery <- unlist(query, use.names=FALSE)
                querynames <- rep(if(is.null(names(query))) as.character(seq_len(length(query))) else names(query),
                                  elementLengths(query))
                querylengths <- unlist(width(query), use.names=FALSE)
                zeroquerynames <- (if(is.null(names(query))) as.character(seq_len(length(query))) else names(query))[elementLengths(query)==0]
            
                ##    TranscriptDb query --------------------------------------------------------------------------
            } else if(inherits(query,"TranscriptDb")) {
                if(is.null(reportLevel))
                    stop("'reportLevel' must be set to a non-NULL value for 'query' of type 'TranscriptDb'")
                message(sprintf("extracting %s regions from TranscriptDb...",reportLevel), appendLF=FALSE)
                if(reportLevel == "gene") {
                    tmpquery <- reduce(exonsBy(query, by="gene"))
                    flatquery <- unlist(tmpquery, use.names=FALSE)
                    querynames <- rep(names(tmpquery), elementLengths(tmpquery))
                    querylengths <- unlist(width(tmpquery), use.names=FALSE)
                    rm(tmpquery)
                
                } else if (reportLevel == "exon") {
                    flatquery <- exons(query, columns="exon_id")
                    querynames <- as.character(mcols(flatquery)$exon_id)
                    querylengths <- width(flatquery)
                
                } else if (reportLevel == "promoter") {
                    flatquery <- promoters(query, columns=c("tx_id","tx_name"))
                    querynames <- paste(as.character(mcols(flatquery)$tx_id),as.character(mcols(flatquery)$tx_name),sep=";")
                    querylengths <- width(flatquery)
                }
                zeroquerynames <- character(0)
                message("done")
            }
            if(length(flatquery)==0)
                stop("'query' is empty - nothing to do")
            ## from now on, only use 'flatquery' (GRanges object) with names in 'querynames' and lengthes in 'querylengths'
        

            ## apply 'mask' to flatquery -------------------------------------------------------------------
            if(!is.null(mask)) {
                if(!inherits(mask,"GRanges"))
                    stop("'mask' must be an object of type 'GRanges'")
                message("removing 'mask' ranges from 'query'...", appendLF=FALSE)
                strand(mask) <- "*"
                mask <- reduce(mask)
                ov <- findOverlaps(flatquery, mask)
                qOM <- unique(queryHits(ov))
                gr1 <- GRanges(seqnames=Rle(qOM), ranges=IRanges(start=start(flatquery)[qOM], end=end(flatquery)[qOM]))
                gr2 <- GRanges(seqnames=Rle(queryHits(ov)), ranges=IRanges(start=start(mask)[subjectHits(ov)], end=end(mask)[subjectHits(ov)]))
                SD <- setdiff(gr1,gr2)
                notOverlappingMaskInd <- which(!((1:length(flatquery)) %in% qOM))
                completelyMaskedInd <- qOM[!(qOM %in% as.numeric(as.character(unique(seqnames(SD)))))]
                
                # 'zeroquery' contains regions that are completely masked (will get zero count and length)
                #zeroquery <- flatquery[completelyMaskedInd]
                tmpnames <- querynames[completelyMaskedInd]
                zeroquerynames <- c(zeroquerynames, tmpnames[!(tmpnames %in% querynames[-completelyMaskedInd])])
                #zeroquerylengths <- rep(0L, length(completelyMaskedInd))

                # masked 'flatquery' maybe split into several non-masked pieces
                flatquery <- c(flatquery[notOverlappingMaskInd],
                               GRanges(seqnames=seqnames(flatquery)[as.numeric(as.character(seqnames(SD)))],
                                       ranges=ranges(SD), strand=strand(flatquery)[as.numeric(as.character(seqnames(SD)))],
                                       seqlengths=seqlengths(flatquery)),
                               ignore.mcols=TRUE)
                querynames <- querynames[c(notOverlappingMaskInd, as.numeric(as.character(seqnames(SD))))]
                querylengths <- width(flatquery)
                message("done")
            }
        
            ## setup tasks for parallelization -------------------------------------------------------------------
            ## TODO: if sum(width(flatquery)) close to sum(seqlengths(genome)) -> select variant counting algorithm (sequential walk through bamfiles)
            if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
                message("preparing to run on ", length(clObj), " nodes...", appendLF=FALSE)
                ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
                if(!all(sapply(ret, function(x) "QuasR" %in% x)))
                    stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
                taskIByFlatQuery <- splitIndices(nx=length(flatquery), ncl=ceiling(length(clObj) /nsamples *2))
                if(inherits(taskIByFlatQuery, "integer", which=FALSE))
                    taskIByFlatQuery <- list(taskIByFlatQuery) # make sure taskIByFlatQuery is a list, even if ceiling(length(clObj) /nsamples *2)==1
                taskSamples <- rep(samples, each=length(taskIByFlatQuery))
                taskBamfiles <- rep(bamfiles, each=length(taskIByFlatQuery))
                flatquery <- lapply(taskIByFlatQuery, function(i) flatquery[i])
                shifts <- rep(shifts, each=length(taskIByFlatQuery))
                myapply <- function(...) {
                    ret <- clusterMap(clObj, ..., SIMPLIFY=FALSE, .scheduling="dynamic")
                    ## fuse
                    iBySample <- split(seq_along(ret),names(ret))[unique(names(ret))]
                    names(ret) <- NULL
                    if(!is.na(proj@snpFile)){
                        ret <- do.call(cbind, lapply(iBySample, function(i) do.call(rbind, ret[i])))
                        postfix <- substring(colnames(ret), nchar(colnames(ret)))
                        ## rename
                        dimnames(ret) <- list(querynames, 
                                              paste(rep(samples, each=3), postfix, sep="_"))
                    } else {
                        ret <- do.call(cbind, lapply(iBySample, function(i) do.call(c, ret[i])))
                        ## rename
                        dimnames(ret) <- list(querynames, samples)
                    }
                    ret
                }
                message("done")
            } else {
                taskSamples <- samples
                taskBamfiles <- bamfiles
                flatquery <- list(flatquery)
                myapply <- function(...) {
                    ret <- do.call(cbind, mapply(..., SIMPLIFY=FALSE))
                    ## rename
                    if(!is.na(proj@snpFile))
                        dimnames(ret) <- list(querynames, 
                                              paste(rep(samples, each=3),
                                                    substring(colnames(ret), nchar(colnames(ret))), sep="_"))
                    else
                        dimnames(ret) <- list(querynames, samples)
                    ret
                }
            }
        
            ## count alignments ----------------------------------------------------------------------------------
            message("counting alignments...", appendLF=FALSE)
            res <- myapply(countAlignments,
                           bamfile=taskBamfiles,
                           regions=flatquery,
                           shift=shifts,
                           MoreArgs=list(
                           selectReadPosition=selectReadPosition,
                           orientation=orientation,
                           useRead=useRead,
                           broaden=broaden,
                           allelic=!is.na(proj@snpFile),
                           includeSpliced=includeSpliced))
            message("done")

        
            ## collapse (sum) counts by sample if necessary
            if(nsamples > length(unique(samples))) {
                if(collapseBySample) {
                    message("collapsing counts by sample...", appendLF=FALSE)
                    if(is.na(proj@snpFile))
                        iBySample <- split(seq_len(nsamples),samples)[unique(samples)]
                    else
                        iBySample <- split(seq_len(ncol(res)),colnames(res))[unique(colnames(res))]
                    res <- do.call(cbind, lapply(iBySample, function(i) rowSums(res[,i,drop=FALSE])))
                    message("done")

                } else {
                    # unify non-collapsed identical sample names
                    if(is.na(proj@snpFile))
                        colnames(res) <- displayNames(proj)
                    else
                        colnames(res) <- paste(rep(displayNames(proj), each=3), substring(colnames(res), nchar(colnames(res))), sep="_")
                }
            }

            ## add the region width as first column
            res <- cbind(width=querylengths, res)
            rm(querylengths)
        
            ## collapse (sum) counts by 'querynames'
            if(length(querynames)>length(unique(querynames))) {
                message("collapsing counts by query name...", appendLF=FALSE)
                iByQuery <- split(seq_len(nrow(res)),querynames)[unique(querynames)]
                res <- do.call(rbind, lapply(iByQuery, function(i) colSums(res[i,1:ncol(res),drop=FALSE])))
                rownames(res) <- querynames <- names(iByQuery)
                message("done")
            }
            if(length(zeroquerynames)>length(unique(zeroquerynames)))
                zeroquerynames <- unique(zeroquerynames)

            ## combine with zeroquery and reorder according to 'query'
            res2 <- matrix(0, nrow=length(querynames)+length(zeroquerynames), ncol=ncol(res),
                           dimnames=list(if(inherits(query,"TranscriptDb"))
                           sort(c(querynames,zeroquerynames))
                           else if(is.null(names(query)))
                           as.character(seq_len(length(query)))
                           else unique(names(query)),
                           colnames(res)))
            res2[rownames(res),] <- res

            ## return results
            return(res2)
        }
    }


## count junctions (with the C-function) for single bamfile and optionally selected target sequences
## return a named vector with junction elements (names of the form "chromosome:first_intronic_base:last_intronic_base:strand")
countJunctionsOneBamfile <- function(bamfile, targets, allelic) {
    tryCatch({ # try catch block goes through the whole function
        # prepare region vectors
        bh <- scanBamHeader(bamfile)[[1]]$targets
        tid <- seq_along(bh) - 1L
        if(!is.null(targets)) {
            if(any(is.na(i <- match(targets,names(bh)))))
                stop(sprintf("some targets not found in bamfile '%s': %s",bamfile,targets[which(is.na(i))]))
            bh <- bh[i]
            tid <- tid[i]
        }
        start <- rep(0L, length(tid)) # samtool library has 0-based inclusiv start
        end <- unname(bh) ## samtool library has 0-based exclusiv end
        # count junctions
        count <- .Call(countJunctions, bamfile, tid, start, end, allelic, PACKAGE="QuasR")
        return(count)

    }, error = function(ex) {
        emsg <- paste("Internal error on", Sys.info()['nodename'], "query bamfile", bamfile, 
                      "\n Error message is:", ex$message)
        stop(emsg)
    })
    
}


## count alignments (with the C-function) for single bamfile, single shift, and single set of regions
## return a numeric vector with length(regions) elements (same order as regions)
countAlignments <- function(bamfile, regions, shift, selectReadPosition, orientation,
                            useRead, broaden, allelic, includeSpliced)
{
    tryCatch({ # try catch block goes through the whole function
        
        ## translate seqnames to tid and create region data.frame
        seqnamesBamHeader <- names(scanBamHeader(bamfile)[[1]]$targets)
        
        ## prepare region vectors
        #tid <- IRanges::as.vector(IRanges::match(seqnames(regions), seqnamesBamHeader)) - 1L
        tid <- as.vector(match(seqnames(regions), seqnamesBamHeader) - 1L) 
        start <- start(regions) - 1L ## samtool library has 0-based inclusiv start
        end <- end(regions) ## samtool library has 0-based exclusiv end
        
        ## swap strand for 'orientation="opposite"' 
        if(orientation == "any")
            strand <- rep("*", length(regions))
        else if(orientation == "opposite")
            strand <- c("+"="-", "-"="+", "*"="*")[as.character(strand(regions))]
        else # orientation == "same"
            strand <- as.character(strand(regions))
        
        ## translate useRead parameter
        BAM_FREAD1 <- 64L
        BAM_FREAD2 <- 128L
        if(useRead == "any")
            readBitMask <- BAM_FREAD1 + BAM_FREAD2
        else if (useRead == "first")
            readBitMask <- BAM_FREAD1
        else if (useRead == "last")
            readBitMask <- BAM_FREAD2
        
        ## get counts
        if(!allelic) {
            count <- .Call(countAlignmentsNonAllelic, bamfile, tid, start, end, strand,
                           selectReadPosition, readBitMask, shift, broaden, includeSpliced,
                           PACKAGE="QuasR")
        } else {
            count <- as.matrix(as.data.frame(.Call(countAlignmentsAllelic, bamfile, tid, start, end, strand,
                                                   selectReadPosition, readBitMask, shift, broaden, includeSpliced,
                                                   PACKAGE="QuasR")))
        }
     
        return(count)
    }, error = function(ex) {
        reg <- regions[c(1, length(regions))]
        emsg <- paste("Internal error on", Sys.info()['nodename'], "query bamfile", bamfile,"with regions\n", 
                      paste(seqnames(reg), start(reg), "-" , end(reg), strand(reg), collapse="\n\t...\n"), 
                      "\n Error message is:", ex$message)
        stop(emsg)
    })
}



## Counts alignments in a given set of regions which are located in a subspace of the genome
## using a per-base-coverage vector approach
##     shift the read, broaden fetch region,
## return a numeric vector with length(regions) elements (same order as regions)
countAlignmentsSubregionsC <- function(bamfile, regions, selectReadPosition, shift=0L, broaden=0L, includeSpliced=TRUE)
{
    ## check if region are located on one chromosme
    seqName <- unique(seqnames(regions))
    if(length(seqName) > 1L)
        stop("regions should only be located on one chromosome")
    
    ## check broaden and shift parameter
    if(broaden < 0)
        stop("'broaden' should not be negative") 
    #    if(shift > 0 && selectReadPosition="midwithin")
    #        stop("'shift' parameter must be zero if 'selectReadPosition' is set to midwithin")
    #    if(broaden > 0 && (selectReadPosition="startwithin" || selectReadPosition="endwithin"))
    #        stop("'broaden' parameter must be zero if 'selectReadPosition' is set to startwithin or endwithin")
    
    ## translate seqName to tid
    seqnamesList <- names(scanBamHeader(bamfile)[[1]]$targets)
    tidList <- as.integer(seq_along(seqnamesList)-1)
    tid <- tidList[ match(seqName, seqnamesList) ]
    
    ## convert grange to data.frame 
    ## with 0-based start inclusive
    ## with 0-based end exclusive
    regionsTable <- data.frame(start=as.integer(start(regions)-1), ## samtool library has 0-based start
                               end=as.integer(end(regions)),
                               strand=as.character(strand(regions)),
                               stringsAsFactors=FALSE
    )
    
    ## call c-function
    cnt <- .Call(countAlignmentsSubregions, bamfile, bamfile, tid, min(regionsTable$start), max(regionsTable$end),
                 regionsTable, as.integer(shift), as.integer(broaden), selectReadPosition, includeSpliced,
                 PACKAGE="QuasR")
    
    return(cnt)
}
