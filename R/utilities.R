displayNames <- function(proj) { # create unique names for each sequence file
    if(!inherits(proj, "qProject", which=FALSE))
        stop("'proj' must be an object of type 'qProject' (returned by 'qAlign')")
    samples <- proj@alignments$SampleName
    ndigits <- nchar(as.character(length(samples)))
    sprintf(paste("s%0",ndigits,"i_%s",sep=""),seq_along(samples),samples)
}

alignmentStats <-
    function(x, collapseBySample=TRUE) {
        # check argument and extract named vector of 'bamfiles'
        if(is.character(x)) {
            bamfiles <- x
            names(bamfiles) <- basename(x)
        } else if (inherits(x, "qProject", which = FALSE)) {
            aln <- alignments(x)
            bamfiles <- c(aln$genome$FileName, unlist(aln$aux, use.names=F))
            iGenome <- seq_along(aln$genome$SampleName)
            names(bamfiles) <- if(collapseBySample)
                c(sprintf("%s:genome",aln$genome$SampleName),
                  sprintf("%s:%s",rep(colnames(aln$aux), each=nrow(aln$aux)),rownames(aln$aux)))
                               else
                c(sprintf("%s:genome",displayNames(x)),
                  sprintf("%s:%s",rep(displayNames(x), each=nrow(aln$aux)),rownames(aln$aux)))
        } else {
            stop("'x' must be of type 'character' or 'qProject'")
        }
        # files exist?
        if(any(i <- ! file.exists(bamfiles)))
            stop(sprintf("cannot access bam files: %s",paste(bamfiles[i],collapse=", ")))
        # call idxstats_bam and collapse counts (seqlength, mapped and unmapped)
        res <- do.call(rbind, lapply(bamfiles, function(bf) {
            tmp <- .Call("idxstats_bam", bf, PACKAGE="QuasR")
            im <- tmp$seqname != "*"
            c(seqlength=sum(as.numeric(tmp$seqlength[im])),
              mapped=sum(as.numeric(tmp$mapped[im])),
              unmapped=tmp$unmapped[!im])
        }))
        # fix unmapped counts for auxiliaries
        if (inherits(x, "qProject", which = FALSE))
            res[-iGenome,'unmapped'] <- rep(res[iGenome,'unmapped'],each=nrow(aln$aux)) - res[-iGenome,'mapped']
        # collapse by sample
        if(collapseBySample) {
            res.old <- res
            tmp <- aggregate(res, list(factor(rownames(res),levels=unique(rownames(res)))), sum)
            res <- as.matrix(tmp[,-1])
            rownames(res) <- as.character(tmp[,1])
            res[,'seqlength'] <- res.old[rownames(res),'seqlength']
        }
        return(res)
    }

freeDiskSpace <- function(path) {
    # return the number of free bytes at 'path' as reported by operating system utilities
    path <- sub("(\\\\|/)$","",path)
    res <- NA
    if(file.exists(path)) {
        path <- tools::file_path_as_absolute(path)
        tryCatch({ # fail silently, just return NA
            if(.Platform$OS.type == "unix") {
                ret <- system2("df", shQuote(path), stdout=TRUE, stderr=TRUE)
                if(grepl("^.*[0-9]+ +[0-9]+ +[0-9]+.*$", perl=TRUE, ret[length(ret)]))
                    res <- as.numeric(sub("^.*[0-9]+ +[0-9]+ +([0-9]+) .*$", "\\1", perl=TRUE, ret[length(ret)])) * 1024
            } else if (.Platform$OS.type == "windows") {
                ret <- shell(paste("dir", shQuote(path)), intern=TRUE, translate=TRUE, mustWork=NA)
                ret <- ret[nchar(ret)>0]
                if(grepl("^.* [0-9,]+ bytes free.*$",ret[length(ret)]))
                    res <- as.numeric(gsub(",","",sub("^.* ([0-9,]+) bytes free.*$","\\1",ret[length(ret)])))
            }}, error=function(ex) NULL, warning=function(ex) NULL)
    } else {
        warning("'", path, "' does not exist")
    }
    return(res)
}

truncPath <-
    function(path, w=getOption('width')) {
        # always print full basename, but shorten path if longer than w
        truncpath <- lapply(path, function(p) {
            if(nchar(p[1])>w)
                file.path(paste(substr(p,1,w-nchar(basename(p))-4),"...",sep=""),basename(p))
            else
                p
        })
        return(unlist(truncpath, use.names=FALSE))
    }

truncString <-
    function(s, w=getOption('width')) {
        # truncate 's' in the middle
        truncs <- lapply(s, function(ss) {
            if((l <- nchar(ss[1])) > w)
                paste(substr(ss, 1, floor(w/2)-2), substr(ss, floor(l-w/2)+2, l), sep="...")
            else
                ss
        })
        return(unlist(truncs, use.names=FALSE))
    }

consolidateFileExtensions <-
    function(filename, compressed=FALSE) {
        # convert various sequence file extension to one representative
        if(compressed)
            fileExtension <- tolower(tools::file_ext(sub("[.](gz|bz2|xz)$","",filename)))
        else
            fileExtension <- tolower(tools::file_ext(filename))
        fileExtension[ fileExtension %in% c("fa", "fna", "fasta") ] <- "fasta"
        fileExtension[ fileExtension %in% c("fq", "fastq") ] <- "fastq"
        return(fileExtension)
    }

compressedFileFormat <-
    function(filename) {
        ifelse(grepl("[.](gz|bz2|xz)$",filename),
               c("gz"="gzip", "bz2"="bzip2", "xz"="xz")[tools::file_ext(filename)],
               "none")
    }

compressFile <-
    function(filename, destname, level=6,
             overwrite=FALSE,  remove=TRUE, BFR.SIZE=1e+07) {
        # gzip/gunzip a file - based on R.utils::gzip.default and R.utils::tar
        if (filename == destname) 
            stop(sprintf("'filename' and 'destname' are identical: %s", filename))
        if (!overwrite && file.exists(destname)) 
            stop(sprintf("'destname' already exists: %s", destname))
        #inn <- file(filename, "rb")
        inn <- switch(compressedFileFormat(filename),
                      "none" =    file(filename, "rb"),
                      "gzip" =  gzfile(filename, "rb"),
                      "bzip2" = bzfile(filename, "rb"),
                      "xz" =    xzfile(filename, "rb"))
        on.exit(if (!is.null(inn)) close(inn))
        outComplete <- FALSE
        out <- switch(compressedFileFormat(destname),
                      "none" =    file(destname, "wb"),
                      "gzip" =  gzfile(destname, "wb", compression = level),
                      "bzip2" = bzfile(destname, "wb", compression = level),
                      "xz" =    xzfile(destname, "wb", compression = level))
        on.exit({
            close(out)
            if (!outComplete) {
                file.remove(destname)
            }
        }, add = TRUE)
        nbytes <- 0
        repeat {
            bfr <- readBin(inn, what = raw(0), size = 1, n = BFR.SIZE)
            n <- length(bfr)
            if (n == 0) 
                break
            nbytes <- nbytes + n
            writeBin(bfr, con = out, size = 1)
        }
        outComplete <- TRUE
        if (remove) {
            close(inn)
            inn <- NULL
            file.remove(filename)
        }
        invisible(nbytes)
    }

grangesFromGff <-
    function(con, version="2", split=NULL) {
        # create a GRangesList object from regions defined in a GFF/GTA file, split by 'type'
        if (requireNamespace("rtracklayer", quietly=TRUE)) {
            gr <- rtracklayer::import.gff(con, version=version, asRangedData=FALSE,
                                          colnames=c("strand", "type", "source", "gene_id", "transcript_id", "exon_number"))

            if(is.null(split) || !(split %in% colnames(mcols(gr))))
                return(gr)
            else
                #return(lapply(split(seq.int(length(gr)), mcols(gr)[[split]]), function(i) gr[i])) # creates a list() of GRanges
                return(split(gr, mcols(gr)[[split]])) # creates a GRangesList (compound elements, e.g. in findOverlaps)
        } else {
            stop("rtracklayer package not found - aborting grangesFromGff")
        }
    }

concatenateFiles <-
    function(filenames, destname, BFR.SIZE=1e+07) {
        # concatenate 'filenames' to 'destname'
        if(any(f <- !(file.exists(filenames))))
            stop("'filenames' don't exist: ", paste(filenames[f], collapse=", "))
        if(!is.character(destname) || length(destname)!=1)
            stop("'destname' needs to be a single output file name")
        if(file.exists(destname))
            stop(sprintf("'destname' already exists: %s", destname))
        outComplete <- FALSE
        out <- file(destname, "wb")
        on.exit({
            close(out)
            if (!outComplete) {
                file.remove(destname)
            }
        })
        nbytes <- 0
        for(filename in filenames) {
            inn <- file(filename, "rb")
            repeat {
                bfr <- readBin(inn, what = raw(0), size = 1, n = BFR.SIZE)
                n <- length(bfr)
                if (n == 0) 
                    break
                nbytes <- nbytes + n
                writeBin(bfr, con = out, size = 1)
            }
            close(inn)
        }
        outComplete <- TRUE
        invisible(nbytes)
    }

md5subsum <-
    function(filenames) {
        # calculate md5sum() on a repoducible random subset of the file's content
        unlist(lapply(filenames, function(fname) {
	    funit <- 1e6
            fs <- floor(file.info(fname)$size /funit)
	    if(fs==0) {
	        funit <- 1
		fs <- floor(file.info(fname)$size)
            }
            if(!is.na(fs)) {
                inn <- file(fname, "rb")
                outname <- tempfile()
                out <- file(outname, "wb")
                set.seed(1)
                for(pos in c(0, funit * sort(sample.int(n=fs, size=20, replace=TRUE)))) {
                    seek(inn, where=pos)
                    bfr <- readBin(inn, what = raw(0), size = 1, n = 10000)
                    if (length(bfr) == 0) 
                        break
                    writeBin(bfr, con = out, size = 1)
                }
                close(inn)
                close(out)
                res <- tools::md5sum(outname)
                unlink(outname)
                return(res)
                
            } else {
                warning(sprintf("could not stat file '%s'; returning 'NA' as md5subsum",fname))
                return(NA)
            }
        }), use.names=FALSE)
    }

# return a list(2) of BiocParallel::BiocParallelParam for nested parallelization
#  [[1]] is used to parallelize across files
#  [[2]] is used to parallelize across threads (accessing a single file) -> has to be a MulticoreParam or SerialParam object
# for downwards compatibility, a parallel::cluster object can be passed that will be used to create the BiocParallel parameter objects
getListOfBiocParallelParam <- function(clObj=NULL) {
    if(is.null(clObj)) { # no 'clObj' argument
        if("package:BiocParallel" %in% search()) { # BiocParallel is loaded
            bppl <- registered() # use registered backend(s)
        } else {                                   # BiocParallel not loaded
            bppl <- list(BiocParallel::SerialParam(), BiocParallel::SerialParam())
        }
    } else {             # have 'clObj' argument
        if(inherits(clObj, "SOCKcluster")) {
            # get node names
            tryCatch({
                nodeNames <- unlist(parallel::clusterEvalQ(clObj, Sys.info()['nodename']))
            }, error = function(ex) {
                message("FAILED")
                stop("The cluster object does not work properly on this system. Please consult the manual of the package 'parallel'\n",call.=FALSE)
            })
            coresPerNode <- table(nodeNames)
            # subset cluster object (represent each node just a single time)
            clObjSub <- clObj[!duplicated(nodeNames)]
            bppl <- if(min(coresPerNode)==1)
                list(as(clObj, "SnowParam"), SerialParam())
            else
                list(as(clObjSub, "SnowParam"), MulticoreParam(workers=min(coresPerNode)))
        } else if(is.list(clObj) && all(sapply(clObj, inherits, "BiocParallelParam"))) {
            bppl <- clObj
        }
    }
    if(length(bppl)==1)
        bppl[[2]] <- BiocParallel::SerialParam()
    if(!inherits(bppl[[2]], c("MulticoreParam","SerialParam")))
        stop('Error configuring the parallel backend. The second registered backend (registered()[[2]] or clObj[[2]]) has to be of class "MulticoreParam" or "SerialParam"')
    return(bppl[1:2])
}
