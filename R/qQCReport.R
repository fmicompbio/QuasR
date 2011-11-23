qQCReport <- function(qproject, pdfFilename=NULL, ...)
{
    if(!is.null(pdfFilename)){
        pdf(pdfFilename)
        on.exit(dev.off())
    }

#     if(require(parallel))
#         cl <- makeCluster(2)
#     #stopCluster(cl)
#     clusterCall(cl, function() library("ShortRead"))

    # FASTQ quality control
    seqChunkSize <- 1e6
    fnames.fastq <- qproject@samples$filepath[qproject@samples$filetype == "fastq"] # TODO allow also non fastq files
    names(fnames.fastq) <- basename(fnames.fastq)

    if(is.null(qproject@qc$qa)){
#         qc1L <- parLapply(cl, seq_along(fnames.fastq),
        qc1L <- lapply(seq_along(fnames.fastq),
                          function(i, sChunkSize, fnames){
                              f <- FastqSampler(fnames[i], sChunkSize)
                              qa(yield(f), names(fnames)[i], type="fastq")
                          },
                          seqChunkSize, fnames.fastq)
        qproject@qc$qa <- do.call(rbind, qc1L)
    }
  
    if(is.null(pdfFilename))
        dev.new()
    qL <- .plotQualByCycle(qproject@qc$qa, ...)
    
    if(is.null(pdfFilename))
        dev.new()
    nL <- .plotNuclByCycle(qproject@qc$qa, ...)
    
    if(is.null(pdfFilename))
        dev.new()
    dL <- .plotDuplicated(qproject@qc$qa, ...)

    if(!is.null(qproject@qc$mappingStats)){
        if(is.null(pdfFilename))
            dev.new()
        mL <- .plotMappings(qproject@qc$mappingStats$genome, ...)
    }

    if(!is.null(qproject@alignments$genome)){
        fnames.bam <- qproject@alignments$genome
        names(fnames.bam) <- basename(qproject@alignments$genome)
        if(is.null(pdfFilename))
            dev.new()
        eL <- .plotErrorsByCycle(fnames.bam, ...)
    }
}

.plotQualByCycle <- function(qcdata, lmat=matrix(1:18, nrow=6, byrow=TRUE)) {
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
        text(x=par('usr')[1]+par('cxy')[1]/4, y=par('usr')[3]+par('cxy')[2]/4, adj=c(0,0), label=rownames(qtiles)[i])
        box()
    }

    invisible(qtilesL)
}

.plotNuclByCycle <- function(qcdata, lmat=matrix(1:18, nrow=6, byrow=TRUE)) {
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
        tmp
    })
    names(nfreqL) <- rownames(nfreq)

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
        text(x=par('usr')[1]+cxy[1]/4, y=par('usr')[4]-cxy[2]/4, adj=c(0,1), label=rownames(nfreq)[i])
        text(x=par('usr')[2]-cxy[1]/4-(4:0)*cxy[1]*0.8, y=par('usr')[4]-cxy[2]/4, adj=c(1,1), col=1:5, label=rownames(nfreqL[[i]]))
    }

    invisible(nfreqL)
}

.plotDuplicated <- function(qcdata, breaks=c(1:10), lmat=matrix(1:6, nrow=3, byrow=TRUE)) {
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
        axis(1, at=1:xn, label=breakNames)
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
            ytop <- par('usr')[4] - 2.0*par('cxy')[2]
            yoff <- ytop - 1.8*cumsum(strheight(frqseqS, cex=frqseqcex, family="mono"))
            ii <- yoff+diff(yoff[1:2]) > max(bocc[i, 1:xn > floor(xleft)])
            if(any(ii)) {
                text(x=xleft, y=ytop,     adj=c(0,0), label=paste(nm,"frequent sequences (ppm):",sep="\n"))
                text(x=xleft, y=yoff[ii], adj=c(0,1), label=paste(frqseqS,frqseqJ,frqseqF,sep="")[ii], family="mono", cex=frqseqcex)
            } else {
                text(x=xleft, y=ytop,     adj=c(0,0), label=nm)
            }
        }
    }

    invisible(bocc)
}

.plotMappings <- function(mapdata, breaks=c(0,1,5,10,50,100), cols=c("#E41A1C","#006D2C","#31A354","#74C476","#BAE4B3","#EDF8E9","#BDBDBD")) {
    mapdata <- mapdata[nrow(mapdata):1,]
    if(breaks[length(breaks)]<Inf)
        breaks <- c(breaks,breaks[length(breaks)]+1,Inf)
    breakNames <- c(as.character(breaks[1:(length(breaks)-2)]),paste(">",breaks[length(breaks)-2],sep=""))

    if(length(cols) != length(breaks)-1) {
        if(breaks[1]<=0)
            cols <- c("#E41A1C",colorRampPalette(c("#006D2C","#31A354","#74C476","#BAE4B3","#EDF8E9"))(length(breaks)-3),"#BDBDBD")
        else
            cols <- c(colorRampPalette(c("#006D2C","#31A354","#74C476","#BAE4B3","#EDF8E9"))(length(breaks)-2),"#BDBDBD")
    }

    obs <- suppressWarnings(as.numeric(colnames(mapdata)))
    if(is.na(obs[length(obs)])) obs[length(obs)] <- Inf

    bin <- findInterval(obs, breaks, all.inside=TRUE)
    mapdataBinned <- matrix(NA, nrow=nrow(mapdata), ncol=length(breaks)-1, dimnames=list(sampleName=rownames(mapdata), hits=breakNames))
    for(i in 1:nrow(mapdata)) {
        tmp <- tapply(mapdata[i,], bin, sum)
        mapdataBinned[i,as.integer(names(tmp))] <- tmp
    }

    par(las=1, mar=c(5,4+7,4,1)+.1, mfrow=c(2,1))
    mp <- barplot(t(mapdataBinned)/1e6, horiz=TRUE, beside=FALSE, col=cols, border=NA, ylim=c(0,nrow(mapdata)+1),
                  main='', xlab='Number of sequences (Mio.)', ylab='', names.arg=rownames(mapdata))
    legend(x=sum(par('usr')[1:2])/2, y=par('usr')[4], xjust=.5, yjust=0, bty='n', x.intersp=0.25,
           title='No. of hits:', fill=cols, ncol=length(cols), legend=breakNames, xpd=NA)
    mp <- barplot(t(mapdataBinned/rowSums(mapdata))*100, horiz=TRUE, beside=FALSE, col=cols, border=NA, ylim=c(0,nrow(mapdata)+1),
                  main='', xlab='Percent of sequences', ylab='', names.arg=rownames(mapdata))
    legend(x=sum(par('usr')[1:2])/2, y=par('usr')[4], xjust=.5, yjust=0, bty='n', x.intersp=0.25,
           title='No. of hits:', fill=cols, ncol=length(cols), legend=breakNames, xpd=NA)

    invisible(mapdataBinned)
}

.plotErrorsByCycle <- function(fnames, N=1e6, lmat=matrix(1:18, nrow=6, byrow=TRUE)) {
    bh <- scanBamHeader(fnames)
    names(bh) <- fnames
    ns <- length(fnames)
    # MD field regexp: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    # description: The MD field aims to achieve SNP/indel calling without looking at the reference. For example, a string `10A5^AC6'
    # means from the leftmost reference base in the alignment, there are 10 matches followed by an A on the reference which
    # is different from the aligned read base; the next 5 reference bases are matches followed by a 2bp deletion from the
    # reference; the deleted sequence is AC; the last 6 bases are matches.
    # remark: two consequtive mismatches are separated by a '0', e.g. '16T0C2'
    .localMDtoErrPos <- function(x) {
        x <- x[grep("[ACGT]",x)]                    # select mismatch alignments
        x <- x[grep("^",x,fixed=TRUE,invert=TRUE)]  # remove alignments with insertions in the read
        mtch <- lapply(strsplit(x,"[ACGT]+"), as.integer)
        errPos <- unlist(lapply(mtch, function(z) cumsum((z+1)[-length(z)])))
        errPos
    }
    #.localMDtoErrPos(c('16T0C2','0A3A20'))

    data <- lapply(fnames, function(fname) {
        # get target sizes
        trg <- bh[[fname]]$targets
        trg <- trg[ order(trg,decreasing=TRUE) ]
        # estimate average number of alignments per base, assuming 50 bytes per alignment and uniform coverage
        ac <- file.info(fname)$size /50 /sum(as.numeric(trg))
        # select chromosomes to load
        chunksize <- round(N/ac)
        selChr <- c(names(trg)[1], names(trg[-1])[cumsum(as.numeric(trg[-1]))+trg[1]<chunksize])
        reg <- GRanges(selChr, IRanges(start=1, end=trg[selChr]))
        # select subset of the last chromosome
        end(reg)[length(reg)] <- min(trg[selChr[length(selChr)]], chunksize - sum(as.numeric(width(reg[-length(reg)]))))
        # read alignments
        param <- ScanBamParam(tag="MD", what=c("cigar"), which=reg, flag=scanBamFlag(isUnmappedQuery=FALSE))
        aln <- scanBam(fname, param=param)
        cvg <- table(unlist(lapply(aln, function(x) cigarToWidth(x$cigar))))
        cumcvg <- rep(sum(cvg),as.integer(names(cvg)[length(cvg)])) - c(rep(0,as.integer(names(cvg)[1])),cumsum(as.numeric(cvg))[-length(cvg)])
        tmp <- table(unlist(lapply(aln, function(x) .localMDtoErrPos(x$tag$MD))))
        nErr <- rep(0,length(cumcvg))
        nErr[as.integer(names(tmp))] <- tmp
        list(cumcvg=cumcvg, nErr=nErr)
    })
    names(data) <- names(fnames)

    layout(lmat)
    for(i in 1:ns) {
        xn <- length(data[[i]]$cumcvg)
        frq <- data[[i]]$nErr /data[[i]]$cumcvg *100
        ym <- max(5,ceiling(max(frq)))
        par(mar=c(5-1,4-1,4-4,2-1)+.1, mgp=c(3-1,1-0.25,0))
        plot(1:xn, frq, type='l', lwd=2, col="#E41A1C", ylim=c(0,ym), xaxs="i", xlim=c(0,xn)+.5, lty=1, pch=20, cex=0.6,
             main="", xlab='Position in read (bp)', ylab='Mismatche bases (% of aligned)')
        abline(h=0, lty=2, col='gray')
        #abline(v=c(12,25), lty=3, col='red')
        legend(x="topleft", bty="n", legend=names(data)[i])
    }

    invisible(data)
}
