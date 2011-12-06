
preprocessReads <- function(filename, filenameMate=NULL,
                            outputFilename=NULL, outputFilenameMate=NULL,
                            Lpattern="", Rpattern="",
                            max.Lmismatch=rep(0:2, c(6,3,100)), max.Rmismatch=rep(0:2, c(6,3,100)),
                            with.Lindels=FALSE, with.Rindels=FALSE,
                            truncateStartBases=NULL, truncateEndBases=NULL, 
                            minLength=14L, nBases=2L, complexity=NULL,
                            nrec=1000000L
                            )
{
    .progressReport("Start filtering", phase=-1)

    # init variables
    pairedSamples <- FALSE
    
    # check parameters
    if(!is.character(filename))
        stop("The Parameter 'filename' must be of type character.")  
    if(is.null(outputFilename))
        outputFilename <- sprintf("%s/%s_prep.%s",
                             dirname(filename),
                             .baseFileName(filename),
                             .fileExtension(filename)
                             )
    if(!is.character(outputFilename))
        stop("The Parameter 'filename' must be of type character.")    
    if(length(filename) != length(outputFilename))
        stop("The parameters 'filename' and 'outputFilename' must habe equal length.")
    if(any(f <- file.exists(outputFilename)))
        stop("The following files exists already: ", paste(outputFilename[f], collapse=", "), ". You have do delete them.")
        
    if(!is.null(filenameMate)){
        pairedSamples <- TRUE
        if(Lpattern != "" || Rpattern != "")
            stop("Removing adapter isn't supported for paired-end sample. The parameters 'Lpattern' and 'Rpattern' must be \"\".")
        if(!is.character(filenameMate))
            stop("Parameter 'filenameMate' must be of type character.")
        if(is.null(outputFilenameMate))
            outputFilenameMate <- sprintf("%s/%s_prep.%s",
                                 dirname(filenameMate),
                                 .baseFileName(filenameMate),
                                 .fileExtension(filenameMate)
                                 )
        if(!is.character(outputFilenameMate))
            stop("Parameter 'outputFilenameMate' must be of type character.")
        if(length(filenameMate) != length(outputFilenameMate))
            stop("The parameters 'filenameMate' and 'outputFilenameMate' must habe equal length.")
        if(length(filename) != length(filenameMate))
            stop("The parameters 'filename' and 'filenameMate' must habe equal length.")
        if(!all(.fileType(filename) == .fileType(filenameMate)))
            stop("The 'filename' and 'filenameMate' must have equal filetype.")
        if(any(f <- file.exists(outputFilenameMate)))
            stop("The following files exists already: ", paste(outputFilenameMate[f], collapse=", "), ". You have do delete them.")
    }

    ## set the ranges for the truncation
    truncateStartBases <- ifelse(is.null(truncateStartBases), 1, truncateStartBases + 1)
    truncateEndBases <- ifelse(is.null(truncateEndBases), -1, -(truncateEndBases + 1))
        
    # TODO configure filters here and not in subfunction

    if(pairedSamples){
        report <- mapply(.preprocessPairedReads, filename=filename, filenameMate=filenameMate,
                            outputFilename=outputFilename, outputFilenameMate=outputFilenameMate,
                            MoreArgs=list(   
                            truncateStartBases=truncateStartBases, truncateEndBases=truncateEndBases, 
                            minLength=minLength, nBases=nBases, complexity=complexity,
                            nrec=nrec))
        colnames(report) <- paste(basename(filename), basename(filenameMate), sep=":")
    
    }else{
        report <- mapply(.preprocessSingleReads, filename=filename, outputFilename=outputFilename,
                            MoreArgs=list(Lpattern=Lpattern, Rpattern=Rpattern,
                            max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                            with.Lindels=with.Lindels, with.Rindels=with.Rindels,
                            truncateStartBases=truncateStartBases, truncateEndBases=truncateEndBases, 
                            minLength=minLength, nBases=nBases, complexity=complexity,
                            nrec=nrec))
        colnames(report) <- basename(filename)
    }
    .progressReport("Successfully finished the filtering.", phase=1)
    return(as.data.frame(report))
}

.preprocessSingleReads <- function(filename, outputFilename,
                            Lpattern, Rpattern,
                            max.Lmismatch, max.Rmismatch,
                            with.Lindels, with.Rindels,
                            truncateStartBases, truncateEndBases, 
                            minLength, nBases, complexity,
                            nrec
                            )
{
    .progressReport(sprintf("Filter file %s", filename))
    
    ## extend R/Lpattern by Ns
    numNs <- 20
    if(nchar(Lpattern) > 0){
        max.Lmismatch <- max.Lmismatch[1:nchar(Lpattern)]
        max.Lmismatch <- c(max.Lmismatch, 1:numNs+max(max.Lmismatch))
        Lpattern <- paste(c(Lpattern, rep("N", numNs)), collapse="")
    }
    if(nchar(Rpattern) > 0){
        max.Rmismatch <- max.Rmismatch[1:nchar(Rpattern)]
        max.Rmismatch <- c(max.Rmismatch, 1:numNs+max(max.Rmismatch))
        Rpattern <- paste(c(Rpattern, rep("N", numNs)), collapse="")
    }
    
    format <- .fileType(filename)

    ## define filters
    ## number of N bases allowed
    if(is.null(nBases))
        nFilt <- nFilter()
    else
        nFilt <- nFilter(nBases)
    ## complexity entropy
    if(is.null(complexity))
        cFilt <- filterLowComplexity()
    else
        cFilt <- filterLowComplexity(complexity)
    ## short read length
    if(is.null(minLength))
        lFilt <- filterLength()
    else
        lFilt <- filterLength(minLength)
    filters <- c(nFilt, cFilt, lFilt)
    ## deactivate filters if NULL
    active(filters) <- c(!is.null(nBases), !is.null(complexity), !is.null(minLength))
    
    report <- data.frame(totalSequence=0, matchTo5pAdaptor=0, matchTo3pAdaptor=0, tooShort=0, tooManyN=0, lowEntropy=0, totalPassed=0)

    if(format=="fasta"){
        mode <- 'w'
        cycle <- 1L
        while(length(chunks <- readFasta(filename, nrec=nrec, skip=(cycle-1)*nrec)) != 0L){
            report$totalSequence <- report$totalSequence + length(chunks)
        
            ## tuncate start or end bases
            chunks <- narrow(chunks, truncateStartBases, truncateEndBases)
        
            ## trim adaptor
            ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, 
                                     max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                     with.Lindels=with.Lindels, with.Rindels=with.Rindels, 
                                     ranges=TRUE)
            report$matchTo5pAdaptor <- report$matchTo5pAdaptor + sum(start(ranges) != 1)
            report$matchTo3pAdaptor <- report$matchTo3pAdaptor + sum(end(ranges) != width(chunks))
            chunks <- narrow(chunks, start(ranges), end(ranges))
        
            ## filter and write short reads
            filterResults <- as.data.frame(evalSeparately(filters, chunks))
            report$tooManyN <- report$tooManyN + sum(!filterResults$CleanNFilter)
            report$tooShort <- report$tooShort + sum(!filterResults$LengthFilter)
            report$lowEntropy <- report$lowEntropy + sum(!filterResults$LowComplexityFilter)
            filter <- as.logical(do.call(pmin.int, filterResults))
            report$totalPassed <- report$totalPassed + sum(filter)
            writeFasta(chunks[filter], outputFilename, mode=mode)
            
            mode <- 'a'
            cycle <- cycle + 1
        }
    }else{
        fs1 <- FastqStreamer(filename, n=nrec)
        on.exit(rm(fs1))
        mode <- 'w'
        while(length(chunks <- yield(fs1)) != 0L){
            report$totalSequence <- report$totalSequence + length(chunks)
        
            ## tuncate start or end bases
            chunks <- narrow(chunks, truncateStartBases, truncateEndBases)
        
            ## trim adaptor
            ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, 
                                     max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                     with.Lindels=with.Lindels, with.Rindels=with.Rindels, 
                                     ranges=TRUE)
            report$matchTo5pAdaptor <- report$matchTo5pAdaptor + sum(start(ranges) != 1)
            report$matchTo3pAdaptor <- report$matchTo3pAdaptor + sum(end(ranges) != width(chunks))
            chunks <- narrow(chunks, start(ranges), end(ranges))
        
            ## filter and write short reads
            filterResults <- evalSeparately(filters, chunks)
            report$tooManyN <- report$tooManyN + sum(!filterResults$CleanNFilter)
            report$tooShort <- report$tooShort + sum(!filterResults$LengthFilter)
            report$lowEntropy <- report$lowEntropy + sum(!filterResults$LowComplexityFilter)
            filter <- as.logical(do.call(pmin.int, filterResults))
            report$totalPassed <- report$totalPassed + sum(filter)
            writeFastq(chunks[filter], outputFilename, mode=mode, qualityType="Auto")
            
            mode <- 'a'
        }

    }
    return(report)
}

.preprocessPairedReads <- function(filename, filenameMate,
                            outputFilename, outputFilenameMate,
                            truncateStartBases, truncateEndBases, 
                            minLength, nBases, complexity,
                            nrec
                            )
{
    .progressReport(sprintf("Filter paired files %s and %s", filename, filenameMate))

    format <- .fileType(filename)
                                   
    ## define filters
    ## number of N bases allowed
    if(is.null(nBases))
        nFilt <- nFilter()
    else
        nFilt <- nFilter(nBases)
    ## complexity entropy
    if(is.null(complexity))
        cFilt <- filterLowComplexity()
    else
        cFilt <- filterLowComplexity(complexity)
    ## short read length
    if(is.null(minLength))
        lFilt <- filterLength()
    else
        lFilt <- filterLength(minLength)
    filters <- c(nFilt, cFilt, lFilt)
    active(filters) <- c(!is.null(nBases), !is.null(complexity), !is.null(minLength))
     
    report <- data.frame(totalSequence=0, matchTo5pAdaptor=NA, matchTo3pAdaptor=NA, tooShort=0, tooManyN=0, lowEntropy=0, totalPassed=0)
    
    if(format=="fasta"){
        mode <- 'w'
        cycle <- 1L
        while(length(chunks <- readFasta(filename, nrec=nrec, skip=(cycle-1)*nrec)) != 0L){
            chunksMate <- readFasta(filenameMate, nrec=nrec, skip=(cycle-1)*nrec)
            report$totalSequence <- report$totalSequence + length(chunks)
            
            chunks <- narrow(chunks, truncateStartBases, truncateEndBases)
            chunksMate <- narrow(chunksMate, truncateStartBases, truncateEndBases)

            ## filter and write short reads
            filterResults <- evalSeparately(filters, chunks)
            filterResultsMate <- evalSeparately(filters, chunksMate)
            report$tooManyN <- report$tooManyN + sum(!(filterResults$CleanNFilter & filterResultsMate$CleanNFilter))
            report$tooShort <- report$tooShort + sum(!(filterResults$LengthFilter & filterResultsMate$LengthFilter))
            report$lowEntropy <- report$lowEntropy + sum(!(filterResults$LowComplexityFilter & filterResultsMate$LowComplexityFilter))
            filter <- as.logical(do.call(pmin.int, c(filterResults, filterResultsMate)))
            report$totalPassed <- report$totalPassed + sum(filter)         
            writeFasta(chunks[filter], outputFilename, mode=mode)
            writeFasta(chunksMate[filter], outputFilenameMate, mode=mode)
            
            mode <- 'a'
            cycle <- cycle + 1
        }
    }else{
        fs1 <- FastqStreamer(filename, n=nrec)
        fs2 <- FastqStreamer(filenameMate, n=nrec)
        on.exit(rm(fs1, fs2))
        mode <- 'w'
        while(length(chunks <- yield(fs1)) != 0L){
            chunksMate <- yield(fs2)
            report$totalSequence <- report$totalSequence + length(chunks)
                
            chunks <- narrow(chunks, truncateStartBases, truncateEndBases)
            chunksMate <- narrow(chunksMate, truncateStartBases, truncateEndBases)
                
            ## filter and write short reads
            filterResults <- evalSeparately(filters, chunks)
            filterResultsMate <- evalSeparately(filters, chunksMate)
            report$tooManyN <- report$tooManyN + sum(!(filterResults$CleanNFilter & filterResultsMate$CleanNFilter))
            report$tooShort <- report$tooShort + sum(!(filterResults$LengthFilter & filterResultsMate$LengthFilter))
            report$lowEntropy <- report$lowEntropy + sum(!(filterResults$LowComplexityFilter & filterResultsMate$LowComplexityFilter))
            filter <- as.logical(do.call(pmin.int, c(filterResults, filterResultsMate)))
            report$totalPassed <- report$totalPassed + sum(filter)            
            writeFastq(chunks[filter], outputFilename, mode=mode, qualityType="Auto")
            writeFastq(chunksMate[filter], outputFilenameMate, mode=mode, qualityType="Auto")
            
            mode <- 'a'
        }
    }
    return(report)
}                
                
                
                
filterLength <- function(threshold=0L, .name = "LengthFilter"){
    #.check_type_and_length(threshold, "numeric", 1)
    srFilter(function(x) {
        width(x) >= threshold
    }, name = .name)
}

filterLowComplexity <- function(threshold=0.5, referenceEntropy=0.9770337, .name = "LowComplexityFilter"){
    srFilter(function(x) {
        # less than half the average entropy per dinucleotide of hg18:
        # H/0.9770337 < 0.5
        # 0.9770337*log(16)/log(2) = 3.908135 bit
        #entropy(x)/0.9770337 >= threshold
        diNucFreq <- dinucleotideFrequency(sread(x))
        if(is.null(dim(diNucFreq)[1] == 1)){
            diNucFreq <- diNucFreq/sum(diNucFreq)
            H <- -sum(diNucFreq * log(diNucFreq, base=16), na.rm = TRUE)
        }else{
            diNucFreq <- diNucFreq/rowSums(diNucFreq)
            H <- -rowSums(diNucFreq * log(diNucFreq, base=16), na.rm = TRUE)
        }
        H/referenceEntropy >= threshold
    }, name = .name)
}
