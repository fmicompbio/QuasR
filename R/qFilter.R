
setGeneric("qFilter", function(subject,
                               Lpattern="", Rpattern="",
                               max.Lmismatch = rep(0:2, c(6,3,100)), max.Rmismatch = rep(0:2, c(6,3,100)),
                               with.Lindels = FALSE, with.Rindels = FALSE,
                               minLength=14L, nBases=2L, complexity=0.5,
                               nrec=1000000L, ...
                               ) standardGeneric("qFilter"))

setMethod("qFilter",
          signature(subject="qProject"),
          function(subject, ...){
              if(subject@paired){
                  files <- subject@samples[c("filepath1", "filepath2")]
                  files <- lapply(seq(nrow(files)), function(i) {
                      qFilter(as.character(files[i,]))
                  })
                  files <- do.call(rbind,files)
                  colnames(files) <- c("filtered1","filtered2")
                  subject@samples <- cbind.data.frame(subject@samples, files, stringsAsFactors=F)
              }else{
                  subject@samples$filtered <- unlist(lapply(subject@samples$filepath, qFilter))
              }
              return(subject)
          })

setMethod("qFilter",
          signature(subject="character"),
          function(subject,
                   Lpattern, Rpattern, max.Lmismatch, max.Rmismatch, with.Lindels, with.Rindels,
                   minLength, nBases, complexity,
                   nrec, ...){

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

              .progressReport("Start filtering", phase=-1)

              if(length(subject) > 2)
                  stop("More then two 'subject' filename.")

              pairedSample <- length(subject) == 2

              format <- unique(.fileType(subject))
              if(length(format) != 1)
                  stop("'subject' must have equal filetype.")

              #if(missing(outputFilenames))
                outputFilenames <- sprintf("%s/%s_filtered.%s",
                                         dirname(subject),
                                         .baseFileName(subject),
                                         .fileExtension(subject)
                                         )
#               if(length(outputFilenames) != length(subject))
#                 error("'outputFilenames' must have equal length as 'subject'")

              ## define filters
              nFilt <- nFilter(nBases) ## number of N bases allowed
              lFilt <- filterLength(minLength) ## short read length
              cFilt <- filterLowComplexity(complexity) ## complexity
              filt <- compose(nFilt, lFilt, cFilt)

              eof <- FALSE
              if(format=="fasta"){
                  mode <- 'w'
                  cycle <- 1L
                  while(length(chunks <- readFasta(subject, nrec=nrec, skip=(cycle-1)*nrec)) != 0L){
                      ## trim adaptor, filter and write short reads
                      if(pairedSample == TRUE){
                          chunksMate <- readFasta(subject[2], nrec=nrec, skip=(cycle-1)*nrec)
                          ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=TRUE)
                          rangesMate <- trimLRPatterns(subject=chunksMate, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=TRUE)
                          ranges[start(ranges) == 0] <- IRanges(1,0)
                          rangesMate[start(rangesMate) == 0] <- IRanges(1,0)
                          chunks <- narrow(chunks, start(ranges), end(ranges))
                          chunksMate <- narrow(chunksMate,start(rangesMate), end(rangesMate))
                          filter <- filt(chunks) & filt(chunksMate)
                          writeFasta(chunks[filter], outputFilenames[1], mode=mode)
                          writeFasta(chunksMate[filter], outputFilenames[2], mode=mode)
                      }else{
                          ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=T)
                          ranges[start(ranges)==0] <- IRanges(1,0)
                          chunks <- narrow(chunks, start(ranges), end(ranges))
                          filter <- filt(chunks)
                          writeFasta(chunks[filter], outputFilenames, mode=mode)
                      }
                      mode <- 'a'
                      cycle <- cycle + 1
                  }
              }else{
                  fs1 <- FastqStreamer(subject[1], n=nrec)
                  if(pairedSample == TRUE)
                      fs2 <- FastqStreamer(subject[2], n=nrec)
                  mode <- 'w'
                  while(length(chunks <- yield(fs1)) != 0L){
                      if(pairedSample == TRUE){
                          chunksMate <- yield(fs2)
                          ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=TRUE)
                          rangesMate <- trimLRPatterns(subject=chunksMate, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=TRUE)
                          ranges[start(ranges) == 0] <- IRanges(1,0)
                          rangesMate[start(rangesMate) == 0] <- IRanges(1,0)
                          chunks <- narrow(chunks, start(ranges), end(ranges))
                          chunksMate <- narrow(chunksMate, start(rangesMate), end(rangesMate))
                          filter <- filt(chunks) & filt(chunksMate)
                          writeFastq(chunks[filter], outputFilenames[1], mode=mode, qualityType=Auto)
                          writeFastq(chunksMate[filter], outputFilenames[2], mode=mode, qualityType=Auto)
                      }else{
                          ranges <- trimLRPatterns(subject=chunks, Lpattern=Lpattern, Rpattern=Rpattern, max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch, with.Lindels=with.Lindels, with.Rindels=with.Rindels, ranges=TRUE)
                          ranges[start(ranges) == 0] <- IRanges(1,0)
                          chunks <- narrow(chunks, start(ranges), end(ranges))
                          filter <- filt(chunks)
                          writeFastq(chunks[filter], outputFilenames, mode=mode, qualityType=Auto)
                      }
                      mode <- 'a'
                  }
                  rm(chunks, ranges, filter, fs1)
                  if(pairedSample == TRUE){
                      rm(chunksMate, rangesMate, fs2)
                  }
              }
              .progressReport("Successfully terminated the filtering.", phase=1)
              return(outputFilenames)
          })

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
