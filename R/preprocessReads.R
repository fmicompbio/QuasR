preprocessReads <-
    function(filename, outputFilename=NULL,
             filenameMate=NULL, outputFilenameMate=NULL,
             truncateStartBases=NULL, truncateEndBases=NULL, 
             Lpattern="", Rpattern="",
             max.Lmismatch=rep(0:2, c(6,3,100)), max.Rmismatch=rep(0:2, c(6,3,100)),
             with.Lindels=FALSE, with.Rindels=FALSE,
             minLength=14L, nBases=2L, complexity=NULL,
             nrec=1000000L, clObj=NULL) {

        ## check parameters
        if(!is.character(filename))
            stop("'filename' must be of type character.")  
        if(is.null(outputFilename)) {
            outputFilename <- paste(tools::file_path_sans_ext(filename), "_prep.", tools::file_ext(filename), sep="")
            message(paste(c("using default 'outputFilename':",truncPath(outputFilename)),collapse="\n "),"\n")
        }
        if(!is.character(outputFilename))
            stop("'outputFilename' must be of type character.")    
        if(length(filename) != length(outputFilename))
            stop("'filename' and 'outputFilename' must have equal length.")
        if(any(f <- file.exists(outputFilename)))
            stop("existing output files: ", paste(outputFilename[f], collapse=", "))
        fileformat <- consolidateFileExtensions(filename, compressed=TRUE)
        if(any(f <- !(fileformat %in% c("fasta","fastq"))))
            stop("unsupported file format (must be one of 'fasta' or 'fastq'): ", paste(filename[f], collapse=", "))
        if(any(fileformat != consolidateFileExtensions(outputFilename, compressed=TRUE)))
            stop("format of 'filename' and 'outputFilename' must be identical")
        filecompr <- compressedFileFormat(outputFilename)
        if(any(fileformat == "fasta") && any(compressedFileFormat(filename) != "none"))
            stop("compressed 'fasta' input is not yet supported")
        
        paired <- !is.null(filenameMate)
        if(paired) {
            paired <- TRUE
            if(Lpattern != "" || Rpattern != "")
                stop("Removing adapters from paired-end samples is not yet supported ('Lpattern' and 'Rpattern' must be set to \"\").")
            if(!is.character(filenameMate))
                stop("'filenameMate' must be of type character.")
            if(is.null(outputFilenameMate)) {
                outputFilenameMate <- paste(tools::file_path_sans_ext(filenameMate), "_prep.", tools::file_ext(filenameMate), sep="")
                message(paste(c("using default 'outputFilenameMate':",truncPath(outputFilenameMate)),collapse="\n "),"\n")
            }
            if(!is.character(outputFilenameMate))
                stop("'outputFilenameMate' must be of type character.")
            if(length(filenameMate) != length(outputFilenameMate))
                stop("'filenameMate' and 'outputFilenameMate' must have equal length.")
            if(length(filename) != length(filenameMate))
                stop("'filename' and 'filenameMate' must have equal length.")
            if(any(f <- file.exists(outputFilenameMate)))
                stop("existing output files: ", paste(outputFilenameMate[f], collapse=", "))
            if(any(fileformat != consolidateFileExtensions(filenameMate, compressed=TRUE)))
                stop("format of 'filename' and 'filenameMate' must be identical")
            if(any(fileformat != consolidateFileExtensions(outputFilenameMate, compressed=TRUE)))
                stop("format of 'filenameMate' and 'outputFilenameMate' must be identical")
            if(any(filecompr != compressedFileFormat(outputFilenameMate)))
                stop("compression format of 'outputFilename' and 'outputFilenameMate' must be identical")
        }

        ## set 'from' and 'to' for truncation (from=1 is the first base, to=-1 is the last one)
        truncateFromBase <- ifelse(is.null(truncateStartBases), 1, truncateStartBases + 1)
        truncateToBase <- ifelse(is.null(truncateEndBases), -1, -(truncateEndBases + 1))
        
        ## initialize filters
        activeFilters <- c(!is.null(nBases),
                           !is.null(complexity),
                           !is.null(minLength))
        nFilt <- nFilter(threshold = ifelse(activeFilters[1], nBases[1], 0L))
        cFilt <- filterLowComplexity(threshold = ifelse(activeFilters[2], complexity[1], 0.5))
        lFilt <- filterLength(threshold = ifelse(activeFilters[3], minLength, 0L))
        filters <- c(nFilt, cFilt, lFilt)
        active(filters) <- activeFilters

        ## parallel execution?
        if(!is.null(clObj) & inherits(clObj, "cluster", which=FALSE)) {
            message("preparing to run on ", length(clObj), " nodes...", appendLF=FALSE)
            myapply <- function(...) do.call(cbind, parallel::clusterMap(clObj, ...))
            # load libraries on nodes
            ret <- clusterEvalQ(clObj, library("QuasR"))
            if(!all(sapply(ret, function(x) "QuasR" %in% x)))
                stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
            # avoid nested parallelization
            nthreads <- clusterEvalQ(clObj, .Call(ShortRead:::.set_omp_threads, 1L))
            on.exit(clusterMap(clObj, function(n) .Call(ShortRead:::.set_omp_threads, n), nthreads))
            message("done")
        } else {
            myapply <- mapply
        }
        
        ## do the filtering
        if(paired){
            #message("start filtering (paired-end mode)")
            filterReport <- myapply(preprocessPairedReads,
                                    filename=filename, filenameMate=filenameMate,
                                    outputFilename=outputFilename, outputFilenameMate=outputFilenameMate,
                                    fileformat=fileformat, filecompr=filecompr,
                                    MoreArgs=list(   
                                    truncateFromBase=truncateFromBase, truncateToBase=truncateToBase, 
                                    filters=filters,
                                    nrec=nrec))
            colnames(filterReport) <- paste(basename(filename), basename(filenameMate), sep=":")
    
        }else{
            #message("start filtering (single read mode)")
            filterReport <- myapply(preprocessSingleReads,
                                    filename=filename,
                                    outputFilename=outputFilename,
                                    fileformat=fileformat, filecompr=filecompr,
                                    MoreArgs=list(
                                    truncateFromBase=truncateFromBase, truncateToBase=truncateToBase, 
                                    Lpattern=Lpattern, Rpattern=Rpattern,
                                    max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                    with.Lindels=with.Lindels, with.Rindels=with.Rindels,
                                    filters=filters,
                                    nrec=nrec))
            colnames(filterReport) <- basename(filename)
        }
        #message("finished filtering")
        return(filterReport)
}

preprocessSingleReads <-
    function(filename, outputFilename, fileformat, filecompr,
             truncateFromBase, truncateToBase, 
             Lpattern, Rpattern,
             max.Lmismatch, max.Rmismatch,
             with.Lindels, with.Rindels,
             filters,
             nrec) {

        message("  filtering ", truncPath(filename, getOption('width')-13))

        ## create (temporary) output file name, will be compressed if filecompr!="none"
        tmpOutputFilename <- if(filecompr!="none") tempfile(pattern="preprocessReadsTemp") else outputFilename

        ## extend R/Lpattern by Ns
        numNs <- 90
        if(nchar(Lpattern) > 0){
            max.Lmismatch <- max.Lmismatch[1:nchar(Lpattern)]
            max.Lmismatch <- c(max.Lmismatch, 1:numNs+max(max.Lmismatch))
            Lpattern <- paste(c(rep("N", numNs), Lpattern), collapse="")
        }
        if(nchar(Rpattern) > 0){
            max.Rmismatch <- max.Rmismatch[1:nchar(Rpattern)]
            max.Rmismatch <- c(max.Rmismatch, 1:numNs+max(max.Rmismatch))
            Rpattern <- paste(c(Rpattern, rep("N", numNs)), collapse="")
        }
    
        ## filter chunks
        filterReport <- c(totalSequences=0, matchTo5pAdapter=0, matchTo3pAdapter=0, tooShort=0, tooManyN=0, lowComplexity=0, totalPassed=0)

        if(fileformat=="fasta"){
            mode <- 'w'
            cycle <- 1L
            while(length(chunks <- readFasta(filename, nrec=nrec, skip=(cycle-1)*nrec)) != 0L){
                filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
        
                ## tuncate start or end bases
                chunks <- narrow(x=chunks, start=truncateFromBase, end=truncateToBase)
        
                ## trim adaptor
                ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, 
                                         max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                         with.Lindels=with.Lindels, with.Rindels=with.Rindels, 
                                         ranges=TRUE)
                filterReport['matchTo5pAdapter'] <- filterReport['matchTo5pAdapter'] + sum(start(ranges) != 1)
                filterReport['matchTo3pAdapter'] <- filterReport['matchTo3pAdapter'] + sum(end(ranges) != width(chunks))
                chunks <- narrow(x=chunks, start=start(ranges), end=end(ranges))
        
                ## filter and write short reads
                filterResults <- evalSeparately(filters, chunks)
                if(is.list(filterResults)) # return value of evalSeparately?
                    filterResults <- do.call(cbind, filterResults)
                #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, 'CleanNFilter'])
                # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
                filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, grep('CleanNFilter',colnames(filterResults))[1]])
                filterReport['tooShort'] <- filterReport['tooShort'] + sum(!filterResults[, 'LengthFilter'])
                filterReport['lowComplexity'] <- filterReport['lowComplexity'] + sum(!filterResults[, 'LowComplexityFilter'])
                filter <- apply(filterResults, 1, all)
                filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)
                if(sum(filter))
                    writeFasta(chunks[filter], tmpOutputFilename, mode=mode, compress=FALSE)

                mode <- 'a'
                cycle <- cycle + 1
            }

        } else if(fileformat=="fastq") {
            fs1 <- FastqStreamer(filename, n=nrec)
            on.exit(close(fs1))
            on.exit(rm(fs1), add=TRUE)
            mode <- 'w'
            while(length(chunks <- yield(fs1)) != 0L){
                filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
        
                ## tuncate start or end bases
                chunks <- narrow(x=chunks, start=truncateFromBase, end=truncateToBase)
        
                ## trim adaptor
                ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, 
                                         max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                         with.Lindels=with.Lindels, with.Rindels=with.Rindels, 
                                         ranges=TRUE)
                filterReport['matchTo5pAdapter'] <- filterReport['matchTo5pAdapter'] + sum(start(ranges) != 1)
                filterReport['matchTo3pAdapter'] <- filterReport['matchTo3pAdapter'] + sum(end(ranges) != width(chunks))
                chunks <- narrow(x=chunks, start=start(ranges), end=end(ranges))
        
                ## filter and write short reads
                filterResults <- evalSeparately(filters, chunks)
                #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, 'CleanNFilter'])
                # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
                filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, grep('CleanNFilter',colnames(filterResults))[1]])
                filterReport['tooShort'] <- filterReport['tooShort'] + sum(!filterResults[, 'LengthFilter'])
                filterReport['lowComplexity'] <- filterReport['lowComplexity'] + sum(!filterResults[, 'LowComplexityFilter'])
                filter <- apply(filterResults, 1, all)
                filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)
                if(sum(filter))
                    writeFastq(chunks[filter], tmpOutputFilename, mode=mode, qualityType="Auto", compress=FALSE)
            
                mode <- 'a'
            }

        } else {
            stop("unknown file format: ", fileformat)
        }

        if(filecompr != "none")
            compressFile(tmpOutputFilename, destname=outputFilename, remove=TRUE)

        return(filterReport)
    }


preprocessPairedReads <-
    function(filename, filenameMate, outputFilename, outputFilenameMate, fileformat, filecompr,
             truncateFromBase, truncateToBase, 
             filters,
             nrec) {

        message("  filtering ", truncPath(filename, getOption('width')-17), " and\n    ",
                truncPath(filenameMate, getOption('width')-16))

        ## create (temporary) output file names, will be compressed if filecompr!="none"
        if(filecompr!="none") {
            tmpOutputFilename <- tempfile(pattern="preprocessReadsTemp")
            tmpOutputFilenameMate <- tempfile(pattern="preprocessReadsTemp")
        } else {
            tmpOutputFilename <- outputFilename
            tmpOutputFilenameMate <- outputFilenameMate
        }

        ## filter chunks
        nrec <- round(nrec/2) # only load half the number of pairs
        filterReport <- c(totalSequences=0, matchTo5pAdapter=NA, matchTo3pAdapter=NA, tooShort=0, tooManyN=0, lowComplexity=0, totalPassed=0)
    
        if(fileformat=="fasta"){
            mode <- 'w'
            cycle <- 1L
            while(length(chunks <- readFasta(filename, nrec=nrec, skip=(cycle-1)*nrec)) != 0L){
                chunksMate <- readFasta(filenameMate, nrec=nrec, skip=(cycle-1)*nrec)
                filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
            
                chunks <- narrow(x=chunks, start=truncateFromBase, end=truncateToBase)
                chunksMate <- narrow(x=chunksMate, start=truncateFromBase, end=truncateToBase)

                ## filter and write short reads
                filterResults <- evalSeparately(filters, chunks)
                if(is.list(filterResults)) # return value of evalSeparately?
                    filterResults <- do.call(cbind, filterResults)
                filterResultsMate <- evalSeparately(filters, chunksMate)
                if(is.list(filterResultsMate)) # return value of evalSeparately?
                    filterResults <- do.call(cbind, filterResultsMate)
                #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, 'CleanNFilter'] & filterResultsMate[, 'CleanNFilter']))
                # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
                filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, grep('CleanNFilter',colnames(filterResults))[1]] &
                                                           filterResultsMate[, grep('CleanNFilter',colnames(filterResultsMate))[1]]))
                filterReport['tooShort'] <- filterReport['tooShort'] + sum(!(filterResults[, 'LengthFilter'] &
                                                                       filterResultsMate[, 'LengthFilter']))
                filterReport['lowComplexity'] <- filterReport['lowComplexity'] + sum(!(filterResults[, 'LowComplexityFilter'] &
                                                                           filterResultsMate[, 'LowComplexityFilter']))
                filter <- apply(cbind(filterResults, filterResultsMate), 1, all)
                filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)         
                if(sum(filter)) {
                    writeFasta(chunks[filter], tmpOutputFilename, mode=mode, compress=FALSE)
                    writeFasta(chunksMate[filter], tmpOutputFilenameMate, mode=mode, compress=FALSE)
                }
            
                mode <- 'a'
                cycle <- cycle + 1
            }

        } else if(fileformat=="fastq") {
            fs1 <- FastqStreamer(filename, n=nrec)
            fs2 <- FastqStreamer(filenameMate, n=nrec)
            on.exit(close(fs1))
            on.exit(close(fs2), add=TRUE)
            on.exit(rm(fs1, fs2), add=TRUE)
            mode <- 'w'
            while(length(chunks <- yield(fs1)) != 0L){
                chunksMate <- yield(fs2)
                filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
                
                chunks <- narrow(x=chunks, start=truncateFromBase, end=truncateToBase)
                chunksMate <- narrow(x=chunksMate, start=truncateFromBase, end=truncateToBase)
                
                ## filter and write short reads
                filterResults <- evalSeparately(filters, chunks)
                filterResultsMate <- evalSeparately(filters, chunksMate)
                #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, 'CleanNFilter'] & filterResultsMate[, 'CleanNFilter']))
                # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
                filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, grep('CleanNFilter',colnames(filterResults))[1]] &
                                                           filterResultsMate[, grep('CleanNFilter',colnames(filterResultsMate))[1]]))
                filterReport['tooShort'] <- filterReport['tooShort'] + sum(!(filterResults[, 'LengthFilter'] &
                                                                       filterResultsMate[, 'LengthFilter']))
                filterReport['lowComplexity'] <- filterReport['lowComplexity'] + sum(!(filterResults[, 'LowComplexityFilter'] &
                                                                           filterResultsMate[, 'LowComplexityFilter']))
                filter <- apply(cbind(filterResults, filterResultsMate), 1, all)
                filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)            
                if(sum(filter)) {
                    writeFastq(chunks[filter], tmpOutputFilename, mode=mode, qualityType="Auto", compress=FALSE)
                    writeFastq(chunksMate[filter], tmpOutputFilenameMate, mode=mode, qualityType="Auto", compress=FALSE)
                }
            
                mode <- 'a'
            }

        } else {
            stop("unknown file format: ", fileformat)
        }
        
        if(filecompr != "none") {
            compressFile(tmpOutputFilename, destname=outputFilename, remove=TRUE)
            compressFile(tmpOutputFilenameMate, destname=outputFilenameMate, remove=TRUE)
        }

        return(filterReport)
    }

filterLength <-
    function(threshold=0L, .name = "LengthFilter"){
        #.check_type_and_length(threshold, "numeric", 1)
        srFilter(function(x) {
            width(x) >= threshold
        }, name = .name)
    }

filterLowComplexity <-
    function(threshold=0.5, referenceEntropy=3.908135, .name = "LowComplexityFilter"){
        srFilter(function(x) {
            # less than half the average entropy per dinucleotide of non-random chromosomes in hg18 (3.908135 bits):
            #   entropy(x)/3.908135 >= threshold
            diNucFreq <- dinucleotideFrequency(sread(x))
            if(is.null(dim(diNucFreq))) {
                diNucFreq <- diNucFreq/sum(diNucFreq)
                H <- -sum(diNucFreq * log2(diNucFreq), na.rm = TRUE)
            } else {
                diNucFreq <- diNucFreq/rowSums(diNucFreq)
                H <- -rowSums(diNucFreq * log2(diNucFreq), na.rm = TRUE)
            }
            H/referenceEntropy >= threshold
        }, name = .name)
    }

