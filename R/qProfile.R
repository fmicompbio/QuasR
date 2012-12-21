## query : GRanges : value is list of matrices:
##                   one list element per sample with matrix:
##                            one(normal) or three(allelic) rows per unique query name and max(width(query)) columns
##                   e.g.   res[['Sample1']]['HCP_TSS','-100']

qProfile <-
    function(proj,
             query,
             upstream=1000,
             downstream=upstream,
             selectReadPosition=c("start", "end"),
             shift=0L,
             orientation=c("any","same","opposite"),
             useRead=c("any","first","last"),
             auxiliaryName=NULL,
             mask=NULL,
             collapseBySample=TRUE,
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
        if(!is.numeric(upstream))
            stop("'upstream' must be of type 'numeric'")
        if(!is.numeric(downstream))
            stop("'downstream' must be of type 'numeric'")
        upstream <- as.integer(upstream)
        downstream <- as.integer(downstream)
        selectReadPosition <- match.arg(selectReadPosition)
        orientation <- match.arg(orientation)
        useRead <- match.arg(useRead)

        ## check shift
        if(shift == "halfInsert") {
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
        if(any(f <- !(seqlevels(query) %in% trCommon)))
            stop(sprintf("sequence levels in 'query' not found in alignment files: %s",
                         paste(seqlevels(query)[f],collapse=", ")))

        ## 'useRead' set but not a paired-end experiment?
        if(useRead != "any" && proj@paired == "no")
            warning("ignoring 'useRead' for single read experiments")


        ## preprocess query -------------------------------------------------------------------------------
        ##    --> extract 'querynames', 'refpos', 'queryWin', 'maxUp', 'maxDown'
        ##        split into 'queryWinL', 'refposL', 'querynamesIntL' by name
        ##    GRanges query -------------------------------------------------------------------------------
        if(inherits(query,"GRanges")) {
            # define 'querynames'
            if(!is.null(names(query)))
                querynames <- names(query)
            else
                querynames <- rep("query",length(query)) # combine all regions
            
            # define 'refpos' and 'queryWin'
            if(length(upstream)==1)
                upstream <- rep(upstream,length(query))
            else if(length(upstream)!=length(query))
                stop(sprintf("the length of 'upstream' (%d) must be one or the same as 'query' (%d)",
                             length(upstream), length(query)))
            if(length(downstream)==1)
                downstream <- rep(downstream,length(query))
            else if(length(downstream)!=length(query))
                stop(sprintf("the length of 'downstream' (%d) must be one or the same as 'query' (%d)",
                             length(downstream), length(query)))

            plusStrand <- as.character(strand(query)) != "-"
            refpos <- ifelse(plusStrand, start(query), end(query))
            queryWin <- GRanges(seqnames(query),
                                IRanges(start=pmax(1, ifelse(plusStrand, refpos -upstream, refpos -downstream)),
                                        end=ifelse(plusStrand, refpos +downstream, refpos +upstream)),
                                strand=strand(query))

            # define 'maxUp' and 'maxDown'
            maxUp <- max(upstream)
            maxDown <- max(downstream)
            if((maxUp+maxDown+1) > 100000L)
                warning(sprintf("profiling over large region (%d basepairs)",maxUp+maxDown+1))

            # split 'queryWin' --> 'queryWinL' and 'refpos' --> 'refposL' by 'querynames'
            if(!is.null(clObj) & inherits(clObj,"cluster",which=FALSE)) {
                profileChunkL <- split(seq_len(length(queryWin)), factor(querynames,levels=unique(querynames)))
                approxNumRegionsPerChunk <- length(queryWin) /(length(clObj)/length(bamfiles))
                profileChunkToTask <- round(cumsum(sapply(profileChunkL,length)) /approxNumRegionsPerChunk)
                taskL <- lapply(split(seq_along(profileChunkL), profileChunkToTask), function(i) do.call(c,profileChunkL[i]))
            } else {
                taskL <- list(seq_len(length(queryWin))) # single task per bam file
            }
            queryWinL <- lapply(taskL, function(i) queryWin[i])
            refposL <- lapply(taskL, function(i) refpos[i])
            querynamesInt <- as.integer(factor(querynames,levels=unique(querynames)))
            querynamesIntL <- lapply(taskL, function(i) querynamesInt[i])

            # calculate coverage
            cvgL <- lapply(seq_along(queryWin), function(i) (maxUp-upstream[i]+1):(maxUp+downstream[i]+1))
            cvg <- do.call(rbind, lapply(split(seq_along(queryWin),factor(querynames,levels=unique(querynames))),
                                         function(i) tabulate(unlist(cvgL[i]))))
            colnames(cvg) <- as.character(seq(-maxUp, maxDown, by=1))
                           
        } else {
            stop("'query' must be an object of type 'GRanges'")
        }
        ## from now on, only use 'queryWinL' (named list of GRanges objects) with reference positions in 'refposL'
        
            
        ## apply 'mask' to queryWinL ----------------------------------------------------------------------------
        if(!is.null(mask)) {
            if(!inherits(mask,"GRanges"))
                stop("'mask' must be an object of type 'GRanges'")
            stop("'mask' is not yet supported for qProfile")
        }        

        
        ## setup tasks for parallelization -------------------------------------------------------------------
        if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
            message("preparing to run on ", length(clObj), " nodes...", appendLF=FALSE)
            ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
            if(!all(sapply(ret, function(x) "QuasR" %in% x)))
                stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
            myapply <- function(...) clusterMap(clObj, ..., SIMPLIFY=FALSE, .scheduling="dynamic")
            message("done")
        } else {
            myapply <- function(...) mapply(..., SIMPLIFY=FALSE)
        }

        
        ## count alignments ----------------------------------------------------------------------------------
        message("profiling alignments...", appendLF=FALSE)
        res <- myapply(profileAlignments,
                       bamfile=rep(bamfiles, each=length(queryWinL)),
                       queryids=querynamesIntL,
                       regions=queryWinL,
                       refpos=refposL,
                       shift=rep(shifts, each=length(queryWinL)),
                       MoreArgs=list(
                           selectReadPosition=selectReadPosition,
                           orientation=orientation,
                           useRead=useRead,
                           broaden=broaden,
                           allelic=!is.na(proj@snpFile),
                           maxUp=maxUp,
                           maxDown=maxDown))
        message("done")
        
        ## fuse by input file, rename and collapse by sample
        if(!is.na(proj@snpFile)) {
            res <- lapply(split(seq_along(res), rep(displayNames(proj), each=length(queryWinL))),
                          function(i) {
                              tmpR <- do.call(rbind, lapply(res[i],"[[","R"))
                              tmpU <- do.call(rbind, lapply(res[i],"[[","U"))
                              tmpA <- do.call(rbind, lapply(res[i],"[[","A"))
                              rownames(tmpR) <- rownames(tmpU) <- rownames(tmpA) <- unique(querynames)
                              list(R=tmpR, U=tmpU, A=tmpA)
                          })
            if(collapseBySample) {
                res <- lapply(split(seq_along(res), factor(samples,levels=unique(samples))),
                              function(i) list(R=Reduce("+", lapply(res[i],"[[","R")),
                                               U=Reduce("+", lapply(res[i],"[[","U")),
                                               A=Reduce("+", lapply(res[i],"[[","A"))))
            }
            nms <- paste(rep(names(res),each=3),c("R","U","A"),sep="_")
            res <- do.call(c, res)
            names(res) <- nms
        } else {
            res <- lapply(split(seq_along(res), rep(displayNames(proj), each=length(queryWinL))),
                          function(i) {
                              tmp <- do.call(rbind, res[i])
                              rownames(tmp) <- unique(querynames)
                              tmp
                          })
            if(collapseBySample) {
                res <- lapply(split(seq_along(res), factor(samples,levels=unique(samples))),
                              function(i) Reduce("+", res[i]))
            }
        }

        ## add the region coverage as first elemnt
        res <- c(list(coverage=cvg), res)
        
        ## return results
        return(res)
    }

## profile alignments (with the C-function) for single bamfile, multiple regions, single shift
## return a numeric vector with maxWidth elements corresponding to the positions within the query regions
profileAlignments <- function(bamfile, queryids, regions, refpos, shift, selectReadPosition,
                              orientation, useRead, broaden, allelic, maxUp, maxDown)
{
    tryCatch({ # try catch block contains whole function
        
        # translate seqnames to tid and create region data.frame
        seqnamesBamHeader <- names(scanBamHeader(bamfile)[[1]]$targets)
        
        # prepare region vectors
        tid <- as.vector(match(seqnames(regions), seqnamesBamHeader) - 1L) 
        s <- start(regions) - 1L # Samtools library has 0-based inclusive start
        e <- end(regions) # Samtools library has 0-based exclusive end
        rp <- refpos - 1L # Samtools library has 0-based inclusive start
        
        ## swap selstrand for 'orientation="opposite"'
        regstrand <- as.character(strand(regions))
        if(orientation == "any")
            selstrand <- rep("*", length(regions))
        else if(orientation == "opposite")
            selstrand <- c("+"="-", "-"="+", "*"="*")[regstrand]
        else # orientation == "same"
            selstrand <- regstrand
        
        ## translate useRead parameter
        BAM_FREAD1 <- 64L
        BAM_FREAD2 <- 128L
        if(useRead == "any")
            readBitMask <- BAM_FREAD1 + BAM_FREAD2
        else if (useRead == "first")
            readBitMask <- BAM_FREAD1
        else if (useRead == "last")
            readBitMask <- BAM_FREAD2

        ## count alignments by position
        if(!allelic) {
            count <- t(.Call(profileAlignmentsNonAllelic, bamfile, queryids, tid, s, e, rp, selstrand, regstrand,
                             selectReadPosition, readBitMask, shift, broaden, maxUp, maxDown, PACKAGE="QuasR"))
            colnames(count) <- as.character(seq(-maxUp, maxDown, by=1))
        } else {
            count <- lapply(.Call(profileAlignmentsAllelic, bamfile, queryids, tid, s, e, rp, selstrand, regstrand,
                                  selectReadPosition, readBitMask, shift, broaden, maxUp, maxDown, PACKAGE="QuasR"),
                            function(x) {
                                x <- t(x)
                                colnames(x) <- as.character(seq(-maxUp, maxDown, by=1))
                                x
                            })
        }
        
        return(count)

    }, error = function(ex) {
        reg <- regions[c(1, length(regions))]
        emsg <- paste("Internal error on ", Sys.info()['nodename'], ", bamfile ", bamfile," with regions\n\t", 
                      paste(seqnames(reg), ":", start(reg), "-" , end(reg), ":", strand(reg), sep="", collapse="\n\t...\n\t"), 
                      "\n Error message is: ", ex$message, sep="")
        stop(emsg)
    })
}
