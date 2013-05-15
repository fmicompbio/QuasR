# proj        : qProject object, fasta file, fastq file or bam file
# pdfFilename : Name of the output file. If NULL then the graph are displayed 
# chunkSize   : chunk/yield size of the sample
# clObj       : cluster object for parallelization

calcQaInformation <- function(filename, label, filetype, chunkSize){
    if(any(filetype == "fasta") && any(compressedFileFormat(filename) != "none")){
        qa <- NULL
        warning("compressed 'fasta' input is not yet supported")
    }else{
        old.o <- options("srapply_fapply"=NULL) # do not run in parallel within qa
        qa <- switch(as.character(filetype),
                     fastq = {
                         f <- FastqSampler(filename, n=chunkSize)
                         reads <- yield(f)
                         close(f)
                         qa(reads, label)
                     },
                     fasta = {
                         reads <- readFasta(as.character(filename), nrec=chunkSize)
                         qa(reads, label)
                     }
                     )
        options(old.o)
    }
    return(qa)
}

calcMmInformation <- function(filename, genome, chunkSize){
    
    # get bamfile index statistics
    stats <- as.data.frame(.Call("idxstatsBam", filename, PACKAGE="QuasR"))
    # get chromosome with mapped reads
    trg <- stats$mapped
    names(trg) <- stats$seqname
    selChr <- names(trg)[ trg != 0 ]

    # get sequence length
    if(is(genome, "BSgenome")) {
        # BSgenome
        ref <- genome
        seqlen <- seqlengths(ref)[selChr]
    } else {
        # Fasta File
        ref <- FaFile(genome)
        seqInfo <- scanFaIndex(ref)
        seqlen <- seqlengths(seqInfo)[selChr]
    }

    # if no mapped alignments then return array of NAs
    if(length(seqlen) == 0){
        mmDist <- list(array(NA,
                             dim=c(5,5,30),
                             dimnames=list(ref=c("A","C","G","T","N"),
                                           read=c("A","C","G","T","N"))),
                       array(NA,
                             dim=c(5,5,30),
                             dimnames=list(ref=c("A","C","G","T","N"),
                                           read=c("A","C","G","T","N"))),
                       rep(NA,100),
                       rep(NA,2))
        return(mmDist)
    }
        
    # create query regions
    if(sum(trg[selChr]) <= chunkSize){
        # number of mapped reads is smaller or equal than chunkSize
        # query all chromosome with mapped reads
        gr <- GRanges(seqnames=names(seqlen), ranges=IRanges(1, seqlen))
    } else {
        # number of mapped reads is bigger than chunkSize
        #  total number of alignments: sum(trg[selChr])
        #  fraction of alignments to be sampled: chunkSize / sum(trg[selChr])
        #  assuming uniform coverage, fraction of chromosome to retrieve: fraction of alignments to be samples * chromosome length
        reflen <- ceiling(chunkSize / sum(trg[selChr]) * seqlen[selChr])
        # use only chromosoms with the most mapped reads
        #seq <- names(sort(trg[selChr], decreasing=T)[1:ceiling(length(trg[selChr]) * chunkSize / sum(trg[selChr]))])
        seq <- selChr
        gr <- unlist(GRangesList(lapply(seq, function(s){
            GRanges(seqnames=s, ranges=breakInChunks(seqlen[s], reflen[s]))
        })))
        rm(seq,reflen)
    }

    # get nucleotide alignment frequencies
    cntFor <- 0
    maxLen <- 0
    # allocate result vector for "nucleotideAlignmentFrequencies" function
    mmDist <- list(integer(25*1000), # mismatch distribution read1 
                   integer(25*1000), # mismatch distribution read2
                   integer(10000), # fragment length distribution
                   integer(2)) # uniqueness (unique/total)
    set.seed(0)
    for(s in sample(length(gr))){
        refseq <- as.character(getSeq(ref, gr[s], as.character=F))
        reftid <- as.integer(match(seqnames(gr[s]), names(trg)) - 1)
        refstart <- start(gr[s])
        len <- .Call("nucleotideAlignmentFrequencies", filename, refseq, reftid, refstart, mmDist, as.integer(chunkSize), PACKAGE="QuasR")
        if(len > maxLen)
            maxLen <- len
        if(sum(mmDist[[1]][1:25]) >= chunkSize || sum(mmDist[[2]][1:25]) >= chunkSize)
            break
        cntFor <- cntFor + 1
    }

    mmDist[[1]] <- array(mmDist[[1]],
                    dim=c(5,5,maxLen),
                    dimnames=list(ref=c("A","C","G","T","N"),
                                  read=c("A","C","G","T","N"),
                                  pos=1:maxLen))
    mmDist[[2]] <- array(mmDist[[2]],
                    dim=c(5,5,maxLen),
                    dimnames=list(ref=c("A","C","G","T","N"),
                                  read=c("A","C","G","T","N"),
                                  pos=1:maxLen))
    names(mmDist[[3]]) <- c(1:(length(mmDist[[3]])-1), paste(">", length(mmDist[[3]])-1, sep=""))
    names(mmDist) <- c("NucleotidMismatchRead1","NucleotidMismatchRead2","FragmentLength", "Uniqueness")
    return(mmDist)
}


qQCReport <- function(input, pdfFilename=NULL, chunkSize=1e6L, clObj=NULL, ...)
{
    # 'proj' is correct type?
    if(inherits(input, "qProject", which=FALSE)){
        filetype <- input@samplesFormat
        if(input@paired == "no"){
            readFilename <- as.character(input@reads$FileName)
            label <- sprintf("%i. %s", 1:nrow(input@reads), basename(as.character(input@reads$FileName)))
            mapLabel <- label
        } else {
            readFilename <- as.character(rbind(input@reads$FileName1, input@reads$FileName2))
            label <- sprintf("%i(R%i). %s",rep(1:nrow(input@reads),each=2), 1:2, basename(as.character(rbind(input@reads$FileName1, input@reads$FileName2))))
            mapLabel <- sprintf("%i. %s",1:nrow(input@reads), basename(as.character(input@reads$FileName1)))

        }
        alnFilename <- input@alignments$FileName
        if(input@genomeFormat == "BSgenome"){
           require(input@genome, character.only=TRUE, quietly=TRUE)
           genome <- get(ls(sprintf("package:%s", input@genome)))
        } else {
            genome <- input@genome
        }
    }else if(is.character(input)){
        filetype <- unique(consolidateFileExtensions(input, compressed=TRUE))
        if(length(filetype) > 1L)
            stop("parameter 'input' must consist of unique filetype")
        label <- sprintf("%i. %s", 1:length(input), basename(input))
        mapLabel <- label
        if(all(file.exists(input))){
            if(filetype != "bam"){
                readFilename <- input
                alnFilename <- NULL
                genome <- NULL
            } else {
                readFilename <- NULL
                alnFilename <- input
                genome <- NULL
            }
        } else {
            stop("could not find the files '", paste(input[!file.exists(input)], collapse=" "), "'")           
        }        
    } else {
        stop("'input' must be an object of type 'qProject' (returned by 'qAlign') or filenames")
    }

    # make sure the cluster nodes are ready
    if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
        message("preparing to run on ", min(length(readFilename),length(clObj)), " nodes...", appendLF=FALSE)
        ret <- clusterEvalQ(clObj, library("QuasR")) # load libraries on nodes
        if(!all(sapply(ret, function(x) "QuasR" %in% x)))
            stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
        message("done")
    }

    message("collecting quality control data")
    # FASTQ/A quality control
    if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE) & length(readFilename)>1) {
        qc1L <- parLapplyLB(clObj,
                            seq(readFilename),
                            function(i){ 
                                calcQaInformation(readFilename[i], label[i], filetype, chunkSize) 
                                })
    } else {
        qc1L <- lapply(seq(readFilename),
                       function(i){
                           calcQaInformation(readFilename[i], label[i], filetype, chunkSize)
                           })
    }
    qa <- do.call(rbind, qc1L)

    # BAM quality control, mismatch distribution
    if(!is.null(alnFilename) && !is.null(genome)){ 
        # eventual TODO: if input is a file list of bam, then this is not processed, because the genome is missing
        if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE) & length(alnFilename)>1) {
            distL <- parLapply(clObj, 
                              alnFilename,
                              calcMmInformation,
                              genome, chunkSize)
        } else {
            distL <- lapply(alnFilename,
                           calcMmInformation,
                           genome, chunkSize)
        }

        if(input@paired == "no") {
            unique <- lapply(distL,"[[", 4)
            frag <- NULL
            mm <- lapply(distL,"[[", 1)
        } else {
            unique <- lapply(distL,"[[", 4)
            frag <- lapply(distL,"[[", 3)
            names(frag) <- mapLabel
            mm <- do.call(c, lapply(distL,"[", c(1,2)))      
        }
        names(unique) <- mapLabel
        names(mm) <- label
    } else {
        unique <- NULL
        mm <- NULL
        frag <- NULL
    }
    
    # open/close pdf stream
    if(!is.null(pdfFilename)){
        pdf(pdfFilename, paper="default", onefile=TRUE, width=0, height=0)
        on.exit(dev.off())
    }

    # Plot
    message("creating QC plots")
    plotdata <- list(raw=list(qa=qa, mm=mm, frag=frag, unique=unique, mapdata=NULL))
    if(!is.null(qa)){
        if(filetype == "fastq"){
            if(is.null(pdfFilename))
                dev.new()
            plotdata[['qualByCycle']] <- plotQualByCycle(qa, ...)
        }
    
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['nuclByCycle']] <- plotNuclByCycle(qa, ...)
    
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['duplicated']] <- plotDuplicated(qa, ...)

    } else if(!is.null(mm)){
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['nuclByCycle']] <- plotNuclByCycle(mm, ...)        
    }

    if(!is.null(alnFilename)){
        # get bamfile index statistics
        mapdata <- lapply(alnFilename, function(file){
            colSums(as.data.frame(.Call("idxstatsBam", file, PACKAGE="QuasR")[c("mapped","unmapped")]))
        })
        mapdata <- do.call(rbind, mapdata)
        rownames(mapdata) <- mapLabel
        plotdata[['raw']][['mapdata']] <- mapdata
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['mappings']] <- plotMappings(mapdata, ...)
    }

    if(!is.null(unique)){
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['uniqueness']] <- plotUniqueness(unique, ...)
    }
        
    if(!is.null(mm)){
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['errorsByCycle']] <- plotErrorsByCycle(mm, ...)

        if(is.null(pdfFilename))
            dev.new()
        plotdata[['mismatchTypes']] <- plotMismatchTypes(mm, ...)
    }

    if(!is.null(frag)){        
        if(is.null(pdfFilename))
            dev.new()
        plotdata[['fragDistribution']] <- plotFragmentDistribution(frag, ...)
    }

    invisible(plotdata)
}

truncStringToPlotWidth <- function(s, plotwidth) {
    sw <- strwidth(s)
    if(any(sw > plotwidth)) {
        w <- 10 # number of character to replace with ".."
        l <- nchar(s)
        news <- s
        i <- sw > plotwidth
        while(w < l && any(i)) {
            news <- ifelse(i, paste(substr(s, 1, ceiling((l-w)/2)+5), substr(s, floor((l+w)/2)+5, l), sep="..."), news)
            sw <- strwidth(news)
            i <- sw > plotwidth
            w <- w + 2
        }
        return(news)
    } else {
        return(s)
    }
}

plotQualByCycle <- function(qcdata, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    data <- qcdata[['perCycle']][['quality']]
    qtiles <- by(list(data$Score, data$Count), list(data$lane, data$Cycle), function(x) {
        coef <- 1.5
        scoreRle <- Rle(x[[1]], x[[2]])
        n <- length(scoreRle)
        nna <- !is.na(scoreRle)
        stats <- c(min(scoreRle), quantile(scoreRle, c(0.25, 0.5, 0.75)), max(scoreRle))
        iqr <- diff(stats[c(2,4)])
        out <- if (!is.na(iqr)) {
            scoreRle < (stats[2L] - coef * iqr) | scoreRle > (stats[4L] + coef * iqr)
        } else !is.finite(scoreRle)
        if (any(out[nna], na.rm = TRUE))
            stats[c(1, 5)] <- range(scoreRle[!out], na.rm = TRUE)
        conf <- stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
        list(stats=stats, n=n, conf=conf, out=which(out))
    }, simplify=FALSE)
    ns <- nrow(qtiles)

    qtilesL <- lapply(1:ns, function(i) {
        tmpconf <- do.call(cbind,lapply(qtiles[i,], "[[", 'conf'))
        tmpxn <- ncol(tmpconf)
        list(stats=do.call(cbind,lapply(qtiles[i,][1:tmpxn], "[[", 'stats')),
             n=sapply(qtiles[i,][1:tmpxn], "[[", 'n'),
             conf=tmpconf,
             out=numeric(0), #sapply(qtiles[i,][1:tmpxn], "[[", 'out'),
             group=numeric(0), #rep(1:ncol(qtiles), sapply(qtiles[i,], function(x) length(x$out))),
             names=colnames(tmpconf))
    })
    names(qtilesL) <- rownames(qtiles)

    layout(lmat)
    for(i in 1:ns) {
        xn <- length(qtilesL[[i]]$names)
        ym <- max(35,max(qtilesL[[i]]$stats))
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(0:1,0:1,type="n", xlab="Position in read (bp)", ylab="Quality score", xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym))
        rect(xleft=seq.int(xn)-0.5, ybottom=-10, xright=seq.int(xn)+0.5, ytop=20,    col=c("#e6afaf","#e6c3c3"), border=NA)
        rect(xleft=seq.int(xn)-0.5, ybottom=20,  xright=seq.int(xn)+0.5, ytop=28,    col=c("#e6d7af","#e6dcc3"), border=NA)
        rect(xleft=seq.int(xn)-0.5, ybottom=28,  xright=seq.int(xn)+0.5, ytop=ym+10, col=c("#afe6af","#c3e6c3"), border=NA)
        do.call("bxp", c(list(qtilesL[[i]], notch=FALSE, width=NULL, varwidth=FALSE, log="", border=par('fg'),
                              pars=list(boxwex=0.8, staplewex=0.5,  outwex=0.5, boxfill="#99999944"),
                              outline=FALSE, horizontal=FALSE, add=TRUE, at=1:xn, axes=FALSE)))
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[3]+cxy[2]/4, adj=c(0,0),
             label=truncStringToPlotWidth(rownames(qtiles)[i], diff(par("usr")[1:2]) - cxy[1]/2))
        box()
    }

    invisible(qtilesL)
}

plotNuclByCycle <- function(qcdata, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    if(!is.null(qcdata[['perCycle']])){
        data <- qcdata[['perCycle']][['baseCall']]
        nfreq <- by(list(data$Base, data$Count), list(data$lane, data$Cycle), function(x) {
            y <- x[[2]] /sum(x[[2]]) *100
            names(y) <- x[[1]]
            y
        }, simplify=TRUE)
        ns <- nrow(nfreq)
    
        nfreqL <- lapply(1:ns, function(i) {
            tmp <- do.call(cbind, lapply(nfreq[i,], function(x) x[c('A','C','G','T','N')]))
            tmp[is.na(tmp)] <- 0
            rownames(tmp) <- c('A','C','G','T','N')
            tmp
        })
        names(nfreqL) <- rownames(nfreq)
    } else {
        nfreqL <- lapply(qcdata, function(x){
            nfreq <- do.call(rbind, lapply(1:dim(x)[3], function(j){
                  c(A=sum(x[,"A",j]),
                       C=sum(x[,"C",j]),
                       G=sum(x[,"G",j]),
                       T=sum(x[,"T",j]),
                       N=sum(x[,"N",j]))
            }))
            rownames(nfreq) <- 1:dim(x)[3]
            t(nfreq / rowSums(nfreq) * 100)
        })
        ns <- length(nfreqL)
#         names(nfreqL) <- names(qcdata)
    }
    nfreqL[is.na.data.frame(nfreqL)] <- 0

    palette(c("#5050ff","#e00000","#00c000","#e6e600","darkgray"))
    layout(lmat)
    for(i in 1:ns) {
        xn <- ncol(nfreqL[[i]])
        ym <- max(50,ceiling(max(nfreqL[[i]]) /5) *5 +5)
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        matplot(1:xn, t(nfreqL[[i]]), type="o", xlab="Position in read (bp)", ylab="Nucleotide frequency (%)",
                xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym), lwd=2, lty=1, pch=20, cex=0.6)
        abline(h=0,lty=2,col="gray")
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1),
             label=truncStringToPlotWidth(names(nfreqL)[i], diff(par('usr')[1:2])-1.5*cxy[1]-strwidth(paste(rownames(nfreqL[[i]]), collapse=" "))))
        text(x=par('usr')[2]-cxy[1]/4-(4:0)*cxy[1]*0.8, y=par('usr')[4]-cxy[2]/4, adj=c(1,1), col=1:5, label=rownames(nfreqL[[i]]))
    }

    invisible(nfreqL)
}

plotDuplicated <- function(qcdata, breaks=c(1:10), lmat=matrix(1:6, nrow=3, byrow=TRUE)) {
    if(breaks[length(breaks)]<Inf)
        breaks <- c(breaks,breaks[length(breaks)]+1,Inf)
    breakNames <- c(as.character(breaks[1:(length(breaks)-2)]),paste(">",breaks[length(breaks)-2],sep=""))
    data <- qcdata[['sequenceDistribution']]
    nocc <- by(list(data$nOccurrences, data$nReads), list(data$lane),
               function(x) Rle(values=x[[1]], lengths=x[[2]]), simplify=FALSE)
    ns <- nrow(nocc)

    nbin <- length(breaks)-1
    bocc <- do.call(rbind, lapply(1:ns, function(i) {
        bin <- findInterval(runValue(nocc[[i]]), breaks)
        tmp <- tapply(runLength(nocc[[i]]),bin,sum) /length(nocc[[i]]) *100
        tmp2 <- numeric(nbin)
        tmp2[ as.integer(names(tmp)) ] <- tmp
        tmp2
    }))
    dimnames(bocc) <- list(sampleName=rownames(nocc), duplicationLevel=breakNames)

    layout(lmat)
    for(i in 1:ns) {
        nm <- rownames(nocc)[i]
        fs <- qcdata[['frequentSequences']]
        fs <- fs[ fs$lane==nm, ]
        xn <- length(breaks)-1
        ym <- max(50,ceiling(max(bocc[i,]) /5) *5)
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(1:xn, bocc[i,], type="o", xlab="Sequence duplication level", ylab="Percent of unique sequences",
             xlim=c(0,xn)+0.5, xaxs="i", ylim=c(0,ym), lwd=2, lty=1, pch=20, cex=0.6, axes=FALSE,
             panel.first=abline(h=0, lty=2, col="gray"))
        axis(1, at=1:xn, labels=breakNames)
        axis(2)
        box()
        frqseqcex <- 0.8
        frqseqS <- sprintf("%-*s",max(nchar(as.character(fs[,'sequence']))),fs[,'sequence'])
        frqseqF <- sprintf("(%6.0f)",fs[,'count']/sum(runValue(nocc[[i]])*as.numeric(runLength(nocc[[i]])))*1e6)
        frqseqJ <- " "
        frqseqW <- max(nchar(as.character(fs[,'sequence'])))
        xleft <- par('usr')[2] - max(strwidth(paste(frqseqS," ",frqseqF,sep=""), cex=frqseqcex, family="mono"))
        while(xleft < 1 && frqseqW > 7) {
            frqseqJ <- ".. "
            frqseqW <- frqseqW - 2
            frqseqS <- strtrim(frqseqS,frqseqW)
            xleft <- par('usr')[2] - max(strwidth(paste(frqseqS,frqseqJ,frqseqF,sep=""), cex=frqseqcex, family="mono"))
        }
        if(xleft >= 1 && frqseqW > 5) {
            cxy <- par('cxy')
            ytop <- par('usr')[4] - 2.0*cxy[2]
            yoff <- ytop - 1.8*cumsum(strheight(frqseqS, cex=frqseqcex, family="mono"))
            ii <- yoff+diff(yoff[1:2]) > max(bocc[i, 1:xn > floor(xleft)])
            if(any(ii)) {
                text(x=xleft, y=ytop,     adj=c(0,0),
                     label=paste(truncStringToPlotWidth(nm,par("usr")[2]-xleft-cxy[1]/2),"frequent sequences (per Mio.):",sep="\n"))
                text(x=xleft, y=yoff[ii], adj=c(0,1),
                     label=paste(frqseqS,frqseqJ,frqseqF,sep="")[ii], family="mono", cex=frqseqcex)
            } else {
                text(x=xleft, y=ytop,     adj=c(0,0),
                     label=truncStringToPlotWidth(nm,par("usr")[2]-xleft-cxy[1]/2))
            }
        }
    }

    invisible(bocc)
}

plotMappings <- function(mapdata, cols=c("#006D2C","#E41A1C"), a4layout=TRUE) {
    nr <- nrow(mapdata)
    
    # set page layout
    if(a4layout)
        layout(rbind(c(0,1),c(0,2),c(0,0)), widths=c(2,3), heights=c(2, nr+2, max(28-nr,0.5)))
    else
        layout(rbind(c(0,1),c(0,2)), widths=c(2,3), heights=c(2, nr+2))

    lapply(seq(1, nr, by=32), function(i) {
        mapdataChunk <- mapdata[min(nr,i+32-1):i, , drop=FALSE]

        if(a4layout && nr>32 && nrow(mapdataChunk)<32)
            mapdataChunk <- rbind(mapdataChunk, matrix(NA, ncol=2, nrow=32-nrow(mapdataChunk)))
        
        # draw legend
        par(mar=c(0,1,0,3)+.1)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
        legend(x="center", xjust=.5, yjust=.5, bty='n', x.intersp=0.25,
               fill=cols, ncol=length(cols), legend=colnames(mapdataChunk), xpd=NA)

        # draw bars
        par(mar=c(5,1,0,3)+.1)
        ymax <- nrow(mapdataChunk)*1.25
        mp <- barplot(t(mapdataChunk/rowSums(mapdataChunk))*100, horiz=TRUE, beside=FALSE, col=cols, border=NA,
                      ylim=c(0,ymax), names.arg=rep("",nrow(mapdataChunk)),
                      main='', xlab='Percent of sequences', ylab='', xpd=NA)

        # draw bar annotation
        cxy <- par('cxy')
        text(x=rep(par('usr')[1]+cxy[1]/3, nrow(mapdataChunk)), y=mp, col="white", adj=c(0,0.5),
             label=sprintf("%.1f%%",mapdataChunk[,'mapped']/rowSums(mapdataChunk)*100), xpd=NA)
        text(x=rep(mean(par('usr')[1:2]), nrow(mapdataChunk)), y=mp, col="white", adj=c(0.5,0.5),
             label=sprintf("total=%.3g",mapdataChunk[,'mapped']+mapdataChunk[,'unmapped']), xpd=NA)
        text(x=rep(par('usr')[2]-cxy[1]/5, nrow(mapdataChunk)), y=mp, col="white", adj=c(1,0.5),
             label=sprintf("%.1f%%",mapdataChunk[,'unmapped']/rowSums(mapdataChunk)*100), xpd=NA)

        # draw sample names
        text(x=par('usr')[1] - 1.0*cxy[1], y=mp, col="black", adj=c(1,0.5),
             label=truncStringToPlotWidth(rownames(mapdataChunk), ((diff(par("usr")[1:2]) + 4*par("cxy")[1]) /3 *2) - 3*par("cxy")[1]), xpd=NA)
    })

    invisible(mapdata)
}

plotUniqueness <- function(data, cols=c("#ff8c00","#4682b4"), a4layout=TRUE) {
    data <- do.call(rbind, data)
    nr <- nrow(data)
 
    data[,2] <- data[,2] - data[,1]
    colnames(data) <- c("unique","non-unique")
    
    # set page layout
    if(a4layout)
        layout(rbind(c(0,1),c(0,2),c(0,0)), widths=c(2,3), heights=c(2, nr+2, max(28-nr,0.5)))
    else
        layout(rbind(c(0,1),c(0,2)), widths=c(2,3), heights=c(2, nr+2))
    
    lapply(seq(1, nr, by=32), function(i) {
        dataChunk <- data[min(nr,i+32-1):i, , drop=FALSE]

        if(a4layout && nr>32 && nrow(dataChunk)<32)
            dataChunk <- rbind(dataChunk, matrix(NA, ncol=2, nrow=32-nrow(dataChunk)))

        # draw legend
        par(mar=c(0,1,0,3)+.1)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
        legend(x="center", xjust=.5, yjust=.5, bty='n', x.intersp=0.25,
               fill=cols, ncol=length(cols), legend=colnames(dataChunk), xpd=NA)
    
        # draw bars
        par(mar=c(5,1,0,3)+.1)
        ymax <- nrow(dataChunk)*1.25
        mp <- barplot(t(dataChunk/rowSums(dataChunk))*100, horiz=TRUE, beside=FALSE, col=cols, border=NA,
                      ylim=c(0,ymax), names.arg=rep("",nrow(dataChunk)),
                      main='', xlab='Percent of unique alignment positions', ylab='', xpd=NA)

        # draw bar annotation
        cxy <- par('cxy')
        text(x=rep(par('usr')[1]+cxy[1]/3, nrow(dataChunk)), y=mp, col="white", adj=c(0,0.5),
             label=sprintf("%.1f%%",dataChunk[,'unique']/rowSums(dataChunk)*100), xpd=NA)
        text(x=rep(mean(par('usr')[1:2]), nrow(dataChunk)), y=mp, col="white", adj=c(0.5,0.5),
             label=sprintf("total=%.3g",dataChunk[,'unique']+dataChunk[,'non-unique']), xpd=NA)
        text(x=rep(par('usr')[2]-cxy[1]/5, nrow(dataChunk)), y=mp, col="white", adj=c(1,0.5),
             label=sprintf("%.1f%%",dataChunk[,'non-unique']/rowSums(dataChunk)*100), xpd=NA)
    
        # draw sample names
        text(x=par('usr')[1] - 1.0*cxy[1], y=mp, col="black", adj=c(1,0.5),
             label=truncStringToPlotWidth(rownames(dataChunk), ((diff(par("usr")[1:2]) + 4*par("cxy")[1]) /3 *2) - 3*par("cxy")[1]), xpd=NA)
    })
    
    invisible(data)
}

plotErrorsByCycle <- function(data, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {

    ns <- length(data)
    layout(lmat)
    for(i in 1:ns) {
        xn <- dim(data[[i]])[3]
        cumcvg <- unlist(lapply(1:xn, function(j){ sum(data[[i]][,,j])}))
        nErr <- cumcvg - unlist(lapply(1:xn, function(j){ sum(diag(data[[i]][,,j]))}))
        frq <- nErr / cumcvg * 100
        ym <- max(5, ceiling(max(frq)), na.rm=T)
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(1:xn, frq, type='l', lwd=2, col="#E41A1C", ylim=c(0,ym), xaxs="i", xlim=c(0,xn)+.5, lty=1, pch=20, cex=0.6,
             main="", xlab='Position in read (bp)', ylab='Mismatche bases (%)')
        abline(h=0, lty=2, col='gray')
        #abline(v=c(12,25), lty=3, col='red')
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1),
             label=truncStringToPlotWidth(names(data)[i], diff(par("usr")[1:2]) - cxy[1]/2))
        if(all(is.na(data[[i]])))
            text(x=mean(par('usr')[1:2]), y=mean(par('usr')[3:4]), adj=c(0.5,0.5), label="no data")
    }

    invisible(data)
}

plotMismatchTypes <- function(data, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    ns <- length(data)
    layout(lmat)
    palette(c("#5050ff","#e00000","#00c000","#e6e600","darkgray"))
    for(i in 1:ns) {
        mtypes <- apply(data[[i]],1,rowSums)
        pmtypes <- mtypes *100 /sum(mtypes)
        diag(pmtypes) <- 0 # don't plot matches
        pmtypes <- pmtypes[,colnames(pmtypes)!="N"] # don't plot genomic N

        barplot(pmtypes, ylab="% of aligned bases", xlab="Genome", col=1:5,
                ylim=c(0,max(colSums(pmtypes),0.1,na.rm=TRUE)*1.16))
        if(all(is.na(mtypes)))
            text(x=mean(par('usr')[1:2]), y=mean(par('usr')[3:4]), adj=c(0.5,0.5), label="no data")
        box()
        cxy <- par("cxy")
        text(x=par('usr')[2]-cxy[1]/4-(4:0)*cxy[1]*0.8, y=par('usr')[4]-cxy[2]/4,
             adj=c(1,1), col=1:5, label=rownames(pmtypes))
        text(x=par('usr')[2]-cxy[1]/4-5*cxy[1]*0.8, y=par('usr')[4]-cxy[2]/4,
             col="black", label="Read:", adj=c(1,1))
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1), col="black",
             label=truncStringToPlotWidth(names(data)[i], par('usr')[2]-cxy[1]/4-5*cxy[1]*0.8-strwidth("Read:")-cxy[1]*1.0))
    }

    invisible(data)
}

plotFragmentDistribution <- function(data, lmat=matrix(1:12, nrow=6, byrow=TRUE)) {
    ns <- length(data)
    frag <- do.call(cbind, data)
    xn <- dim(frag)[1]
    # trim vector
    if(sum(frag[xn,], na.rm=T) == 0L){
        xn <- max(as.integer(rownames(frag)[rowSums(frag, na.rm=T) > 0]), 100)
        frag <- frag[1:xn,, drop=F]
    } 
      
    layout(lmat)
    for(i in 1:ns) {
        dens <- frag[,i]/sum(frag[,i])
        ym <- max(dens, 0.01, na.rm=T)
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))

        #plot(dens, type='l', lwd=2, col="#377EB8", ylim=c(0,ym), xaxs="i", xlim=c(0,xn)+.5, lty=1, pch=20, cex=0.6,
        #     main="", xlab='Fragment size (nt)', ylab='Density')
        plot(dens, type='n', ylim=c(0,ym), xaxs="i", xlim=c(0,xn)+.5, lty=1, pch=20, cex=0.6,
             main="", xlab='Fragment size (nt)', ylab='Density', xpd=NA)
        rect(xleft=1:xn-.5, ybottom=0, xright=1:xn+.5, ytop=dens, col="#377EB8", border=NA)
        abline(h=0, lty=2, col='gray')
        if(any(ii <- grepl("^>",names(dens))))
            axis(side=1, at=xn, labels=names(dens)[xn])
        cxy <- par('cxy')
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1),
             label=truncStringToPlotWidth(names(data)[i], diff(par("usr")[1:2]) - cxy[1]/2))
        medlen <- median(Rle(values=1:xn, lengths=frag[,i]))
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-5*cxy[2]/4, adj=c(0,1), col="#377EB8", label=sprintf("median = %1.f",medlen))
        abline(v=medlen, lty=3, col="black")
        if(all(is.na(data[[i]])))
            text(x=mean(par('usr')[1:2]), y=mean(par('usr')[3:4]), adj=c(0.5,0.5), label="no data")
    }
    
    invisible(frag)
}

.qa_ShortRead <-
    function(dirPath, lane, ..., verbose=FALSE)
{
    if (missing(lane))
        stop("Paramteter 'lane' is missing.")
    obj <- dirPath
    alf <- alphabetFrequency(sread(obj), baseOnly=TRUE, collapse=TRUE)
#     bqtbl <- alphabetFrequency(quality(obj), collapse=TRUE)
#     rqs <- .qa_qdensity(quality(obj))
    freqtbl <- tables(sread(obj))
    abc <- alphabetByCycle(obj)
    names(dimnames(abc)) <- c("base", "cycle")
    dimnames(abc)$cycle <- as.character(1:dim(abc)[2])
    ac <- ShortRead:::.qa_adapterContamination(obj, lane, ...)
    perCycleBaseCall <- data.frame(Cycle = as.integer(colnames(abc)[col(abc)]), 
        Base = factor(rownames(abc)[row(abc)]), Count = as.vector(abc), 
        lane = lane, row.names = NULL)
    perCycleBaseCall <- perCycleBaseCall[perCycleBaseCall$Count != 0, ]
#     perCycleBaseCall <- ShortRead:::.qa_perCycleBaseCall(abc, lane)
#     perCycleQuality <- .qa_perCycleQuality(abc, quality(obj), lane)
    lst <- list(readCounts=data.frame(
           read=length(obj), filter=NA, aligned=NA,
           row.names=lane),
         baseCalls=data.frame(
           A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
           N=alf[["other"]], row.names=lane),
          readQualityScore=data.frame(
           quality=NULL, #rqs$x,
           density=NULL, #rqs$y,
           lane=NULL, #lane,
           type=NULL #"read"
           ),
          baseQuality=data.frame(
           score=NULL, #names(bqtbl),
           count=NULL, #as.vector(bqtbl),
           lane=NULL #lane
           ),
         alignQuality=data.frame(
           score=as.numeric(NA),
           count=as.numeric(NA),
           lane=lane, row.names=NULL),
         frequentSequences=data.frame(
           sequence=names(freqtbl$top),
           count=as.integer(freqtbl$top),
           type="read",
           lane=lane),
         sequenceDistribution=cbind(
           freqtbl$distribution,
           type="read",
           lane=lane),
          perCycle=list(
            baseCall=perCycleBaseCall,
            quality=NULL #perCycleQuality
              ),
         perTile=list(
           readCounts=data.frame(
             count=integer(0), type=character(0),
             tile=integer(0), lane=character(0)),
           medianReadQualityScore=data.frame(
             score=integer(), type=character(), tile=integer(),
             lane=integer(), row.names=NULL)),
         adapterContamination=ac

         )

    ShortRead:::.ShortReadQQA(lst)
}

setMethod(qa, "ShortRead", .qa_ShortRead)
