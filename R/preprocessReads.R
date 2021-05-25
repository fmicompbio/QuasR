#' Preprocess Short Reads
#' 
#' Truncate sequences, remove parts matching to adapters and filter out low 
#' quality or low complexity sequences from (compressed) 'fasta' or 'fastq' files.
#' 
#' Sequence files can be in fasta or fastq format, and can be compressed by
#' either gzip, bzip2 or xz (extensions .gz, .bz2 or .xz). Multiple files
#' can be processed by a single call to \code{preprocessReads}; in that
#' case all sequence file vectors must have identical lengths.
#' 
#' \code{nrec} can be used to limit the memory usage when processing
#' large input files. \code{preprocessReads} iteratively loads chunks of
#' \code{nrec} sequences from the input until all data been processed.
#' 
#' Sequence pairs from paired-end experiments can be processed by
#' specifying pairs of input and output files (\code{filenameMate} and
#' \code{outputFilenameMate} arguments). In that case, it is assumed that
#' pairs appear in the same order in the two input files, and only pairs
#' in which both reads pass all filtering criteria are written to the
#' output files, maintaining the consistent ordering.
#' 
#' If output files are compressed, the processed sequences are first
#' written to temporary files (created in the same directory as the final
#' output file), and the output files are generated at the end by compressing
#' the temporary files.
#' 
#' For the trimming of left and/or right flanking sequences (adapters) from
#' sequence reads, the \code{\link[Biostrings]{trimLRPatterns}} function
#' from package \pkg{Biostrings} is used, and the arguments \code{Lpattern},
#' \code{Rpattern}, \code{max.Lmismatch}, \code{max.Rmismatch},
#' \code{with.Lindels} and \code{with.Rindels} are used in the call to
#' \code{trimLRPatterns}. \code{Lfixed} and \code{Rfixed} arguments
#' of \code{trimLRPatterns} are set to \code{TRUE}, thus only fixed
#' patterns (without IUPAC codes for ambigous bases) can be
#' used. Currently, trimming of adapters is only supported for single read
#' experiments.
#' 
#' Sequence complexity (\eqn{H}) is calculated based on the dinucleotide
#' composition using the formula (Shannon entropy): \deqn{H = -\sum_i {f_i \log_2 f_i},}
#' where \eqn{f_i} is the fraction of dinucleotide \eqn{i} from all
#' dinucleotides in the sequence. Sequence reads that fulfill the condition
#' \eqn{H/H_r \ge c} are retained (not filtered out), where \eqn{H_r =
#' 3.908} is the reference complexity in bits obtained from the human
#' genome, and \eqn{c} is the value given to the argument \code{complexity}.
#' 
#' If an object that inherits from class \code{cluster} is provided to
#' the \code{clObj} argument, for example an object returned by
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel},
#' multiple files will be processed in parallel using
#' \code{\link[parallel:clusterApply]{clusterMap}} from package \pkg{parallel}.
#' 
#' @param filename the name(s) of the input sequence file(s).
#' @param outputFilename the name(s) of the output sequence file(s).
#' @param filenameMate for paired-end experiments, the name(s) of the
#'   input sequence file(s) containing the second read (mate) of each pair.
#' @param outputFilenameMate for paired-end experiments, the name(s) of the
#'   output sequence file(s) containing the second read (mate) of each pair.
#' @param truncateStartBases integer(1): the number of bases to be truncated 
#'   (removed) from the beginning of each sequence.
#' @param truncateEndBases integer(1): the number of bases to be truncated 
#'   (removed) from the end of each sequence.
#' @param Lpattern character(1): the left (5'-end) adapter sequence.
#' @param Rpattern character(1): the right (3'-end) adapter sequence.
#' @param max.Lmismatch mismatch tolerance when searching for matches of 
#'   \code{Lpattern} (see \sQuote{Details}).
#' @param max.Rmismatch mismatch tolerance when searching for matches of 
#'   \code{Rpattern} (see \sQuote{Details}).
#' @param with.Lindels if \code{TRUE}, indels are allowed in the alignments of 
#'   the suffixes of \code{Lpattern} with the subject, at its beginning 
#'   (see \sQuote{Details}).
#' @param with.Rindels same as \code{with.Lindels} but for alignments of the 
#'   prefixes of \code{Rpattern} with the subject, at its end (see 
#'   \sQuote{Details}).
#' @param minLength integer(1): the minimal allowed sequence length.
#' @param nBases integer(1): the maximal number of Ns allowed per sequence.
#' @param complexity \code{NULL} (default) or numeric(1): If not \code{NULL},
#'   the minimal sequence complexity, as a fraction of the average complexity
#'   in the human genome (~3.9bits). For example, \code{complexity = 0.5} will
#'   filter out sequences that do not have at least half the complexity of the
#'   human genome. See \sQuote{Details} on how the complexity is calculated.
#' @param nrec integer(1): the number of sequence records to read at a time.
#' @param clObj a cluster object to be used for parallel processing of multiple 
#' files (see \sQuote{Details}).
#' 
#' @name preprocessReads
#' @aliases preprocessReads
#' 
#' @returns 
#' A matrix with summary statistics on the processed sequences, containing:
#' \itemize{
#'   \item One column per input file (or pair of input files for paired-end
#'                                    experiments).
#'   \item The number of sequences or sequence pairs in rows:
#'   \itemize{
#'     \item \code{totalSequences} - the total number in the input
#'     \item \code{matchTo5pAdaptor} - matching to \code{Lpattern}
#'     \item \code{matchTo3pAdaptor} - matching to \code{Rpattern}
#'     \item \code{tooShort} - shorter than \code{minLength}
#'     \item \code{tooManyN} - more than \code{nBases} Ns
#'     \item \code{lowComplexity} - relative complexity below \code{complexity}
#'     \item \code{totalPassed} - the number of sequences/sequence pairs
#'       that pass all filtering criteria and were written to the output file(s).
#'   }
#' }
#' 
#' @author Anita Lerch, Dimos Gaidatzis and Michael Stadler
#' 
#' @keywords utilities misc 
#' 
#' @export
#' 
#' @seealso \code{\link[Biostrings]{trimLRPatterns}} from package \pkg{Biostrings},
#' \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
#' 
#' @importFrom parallel clusterEvalQ clusterMap
#' @importFrom ShortRead nFilter SRFilterResult
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom S4Vectors active
#' 
#' @examples 
#' # sample files
#' infiles <- system.file(package="QuasR", "extdata",
#'                        c("rna_1_1.fq.bz2","rna_1_2.fq.bz2"))
#' outfiles <- paste(tempfile(pattern=c("output_1_","output_2_")),".fastq",sep="")
#' # single read example
#' preprocessReads(infiles, outfiles, nBases=0, complexity=0.6)
#' unlink(outfiles)
#' # paired-end example
#' preprocessReads(filename=infiles[1],
#'                 outputFilename=outfiles[1],
#'                 filenameMate=infiles[2],
#'                 outputFilenameMate=outfiles[2],
#'                 nBases=0, complexity=0.6)
#' unlink(outfiles)
#'
preprocessReads <- function(filename, outputFilename = NULL,
                            filenameMate = NULL, outputFilenameMate = NULL,
                            truncateStartBases = NULL, truncateEndBases = NULL, 
                            Lpattern = "", Rpattern = "",
                            max.Lmismatch = rep(0:2, c(6, 3, 100)), 
                            max.Rmismatch = rep(0:2, c(6, 3, 100)),
                            with.Lindels = FALSE, with.Rindels = FALSE,
                            minLength = 14L, nBases = 2L, complexity = NULL,
                            nrec = 1000000L, clObj = NULL) {
    
    ## check parameters
    if (!is.character(filename))
        stop("'filename' must be of type character.")  
    if (is.null(outputFilename)) {
        outputFilename <- paste(tools::file_path_sans_ext(filename), "_prep.", 
                                tools::file_ext(filename), sep = "")
        message(paste(c("using default 'outputFilename':",
                        truncPath(outputFilename)), collapse = "\n "), "\n")
    }
    if (!is.character(outputFilename))
        stop("'outputFilename' must be of type character.")    
    if (length(filename) != length(outputFilename))
        stop("'filename' and 'outputFilename' must have equal length.")
    if (any(f <- !file.exists(filename)))
        stop("non-existing input file(s): ", 
             paste(filename[f], collapse = ", "))
    if (any(f <- file.exists(outputFilename)))
        stop("existing output file(s): ", 
             paste(outputFilename[f], collapse = ", "))
    fileformat <- consolidateFileExtensions(filename, compressed = TRUE)
    if (any(f <- !(fileformat %in% c("fasta","fastq"))))
        stop("unsupported file format (must be one of 'fasta' or 'fastq'): ", 
             paste(filename[f], collapse = ", "))
    if (any(fileformat != consolidateFileExtensions(outputFilename, compressed = TRUE)))
        stop("format of 'filename' and 'outputFilename' must be identical")
    filecompr <- compressedFileFormat(outputFilename)
    if (any(fileformat == "fasta") && any(compressedFileFormat(filename) != "none"))
        stop("compressed 'fasta' input is not yet supported")
    
    paired <- !is.null(filenameMate)
    if (paired) {
        paired <- TRUE
        if (Lpattern != "" || Rpattern != "")
            stop("Removing adapters from paired-end samples is not yet ", 
                 "supported ('Lpattern' and 'Rpattern' must be set to \"\").")
        if (!is.character(filenameMate))
            stop("'filenameMate' must be of type character.")
        if (is.null(outputFilenameMate)) {
            outputFilenameMate <- paste(tools::file_path_sans_ext(filenameMate), 
                                        "_prep.", tools::file_ext(filenameMate), sep = "")
            message(paste(c("using default 'outputFilenameMate':",
                            truncPath(outputFilenameMate)), collapse = "\n "), "\n")
        }
        if (!is.character(outputFilenameMate))
            stop("'outputFilenameMate' must be of type character.")
        if (length(filenameMate) != length(outputFilenameMate))
            stop("'filenameMate' and 'outputFilenameMate' must have equal length.")
        if (length(filename) != length(filenameMate))
            stop("'filename' and 'filenameMate' must have equal length.")
        if (any(f <- file.exists(outputFilenameMate)))
            stop("existing output files: ", paste(outputFilenameMate[f], collapse = ", "))
        if (any(fileformat != consolidateFileExtensions(filenameMate, compressed = TRUE)))
            stop("format of 'filename' and 'filenameMate' must be identical")
        if (any(fileformat != consolidateFileExtensions(outputFilenameMate, compressed = TRUE)))
            stop("format of 'filenameMate' and 'outputFilenameMate' must be identical")
        if (any(filecompr != compressedFileFormat(outputFilenameMate)))
            stop("compression format of 'outputFilename' and 'outputFilenameMate' must be identical")
    }
    
    ## set 'from' and 'to' for truncation (from=1 is the first base, to=-1 is the last one)
    truncateFromBase <- ifelse(is.null(truncateStartBases), 1, truncateStartBases + 1)
    truncateToBase <- ifelse(is.null(truncateEndBases), -1, -(truncateEndBases + 1))
    
    ## initialize filters
    activeFilters <- c(!is.null(nBases),
                       !is.null(complexity),
                       !is.null(minLength))
    nFilt <- ShortRead::nFilter(threshold = ifelse(activeFilters[1], nBases[1], 0L))
    cFilt <- filterLowComplexity(threshold = ifelse(activeFilters[2], complexity[1], 0.5))
    lFilt <- filterLength(threshold = ifelse(activeFilters[3], minLength, 0L))
    filters <- c(nFilt, cFilt, lFilt)
    S4Vectors::active(filters) <- activeFilters
    
    ## parallel execution?
    if (!is.null(clObj) & inherits(clObj, "cluster", which = FALSE)) {
        message("preparing to run on ", length(clObj), " nodes...", appendLF = FALSE)
        myapply <- function(...) do.call(cbind, parallel::clusterMap(clObj, ...))
        # load libraries on nodes
        ret <- parallel::clusterEvalQ(clObj, library("QuasR"))
        if (!all(sapply(ret, function(x) "QuasR" %in% x)))
            stop("'QuasR' package could not be loaded on all nodes in 'clObj'")
        # avoid nested parallelization
        nthreads <- parallel::clusterEvalQ(clObj, .Call(ShortRead:::.set_omp_threads, 1L))
        on.exit(parallel::clusterMap(clObj, function(n) .Call(ShortRead:::.set_omp_threads, n), nthreads))
        message("done")
    } else {
        myapply <- mapply
    }
        
    ## do the filtering
    if (paired) {
        #message("start filtering (paired-end mode)")
        filterReport <- myapply(preprocessPairedReads,
                                filename = filename, filenameMate = filenameMate,
                                outputFilename = outputFilename, 
                                outputFilenameMate = outputFilenameMate,
                                fileformat = fileformat, filecompr = filecompr,
                                MoreArgs = list(   
                                    truncateFromBase = truncateFromBase, 
                                    truncateToBase = truncateToBase, 
                                    filters = filters,
                                    nrec = nrec))
        colnames(filterReport) <- paste(basename(filename), 
                                        basename(filenameMate), sep = ":")
    } else {
        #message("start filtering (single read mode)")
        filterReport <- myapply(preprocessSingleReads,
                                filename = filename,
                                outputFilename = outputFilename,
                                fileformat = fileformat, filecompr = filecompr,
                                MoreArgs = list(
                                    truncateFromBase = truncateFromBase, 
                                    truncateToBase = truncateToBase, 
                                    Lpattern = Lpattern, Rpattern = Rpattern,
                                    max.Lmismatch = max.Lmismatch, 
                                    max.Rmismatch = max.Rmismatch,
                                    with.Lindels = with.Lindels, 
                                    with.Rindels = with.Rindels,
                                    filters = filters,
                                    nrec = nrec))
        colnames(filterReport) <- basename(filename)
    }
    #message("finished filtering")
    return(filterReport)
}

#' @keywords internal
#' @importFrom ShortRead readFasta writeFastq FastqStreamer yield narrow 
#'   SRFilterResult
#' @importFrom Biostrings trimLRPatterns
#' @importFrom S4Vectors evalSeparately
#' @importFrom BiocGenerics width start end
preprocessSingleReads <- function(filename, outputFilename, fileformat, filecompr,
                                  truncateFromBase, truncateToBase, 
                                  Lpattern, Rpattern,
                                  max.Lmismatch, max.Rmismatch,
                                  with.Lindels, with.Rindels,
                                  filters,
                                  nrec) {
    
    message("  filtering ", truncPath(filename, getOption('width') - 13))
    
    ## create (temporary) output file name, will be compressed if filecompr!="none"
    tmpOutputFilename <- if (filecompr != "none") tempfile(pattern = "preprocessReadsTemp") 
    else outputFilename
    
    ## extend R/Lpattern by Ns
    numNs <- 90
    if (nchar(Lpattern) > 0) {
        max.Lmismatch <- max.Lmismatch[seq_len(nchar(Lpattern))]
        max.Lmismatch <- c(max.Lmismatch, seq_len(numNs) + max(max.Lmismatch))
        Lpattern <- paste(c(rep("N", numNs), Lpattern), collapse = "")
    }
    if (nchar(Rpattern) > 0) {
        max.Rmismatch <- max.Rmismatch[seq_len(nchar(Rpattern))]
        max.Rmismatch <- c(max.Rmismatch, seq_len(numNs) + max(max.Rmismatch))
        Rpattern <- paste(c(Rpattern, rep("N", numNs)), collapse = "")
    }
    
    ## filter chunks
    filterReport <- c(totalSequences = 0, matchTo5pAdapter = 0, 
                      matchTo3pAdapter = 0, tooShort = 0, tooManyN = 0,
                      lowComplexity = 0, totalPassed = 0)

    if (fileformat == "fasta") {
        mode <- 'w'
        cycle <- 1L
        while (length(chunks <- ShortRead::readFasta(filename, nrec = nrec, 
                                                     skip = (cycle - 1) * nrec)) != 0L) {
            filterReport['totalSequences'] <- filterReport['totalSequences'] + 
                length(chunks)
            
            ## truncate start or end bases
            chunks <- ShortRead::narrow(x = chunks, start = truncateFromBase, 
                                        end = truncateToBase)
            
            ## trim adaptor
            ranges <- Biostrings::trimLRPatterns(
                subject = chunks, Lpattern = Lpattern, Rpattern = Rpattern, 
                max.Lmismatch = max.Lmismatch, max.Rmismatch = max.Rmismatch,
                with.Lindels = with.Lindels, with.Rindels = with.Rindels, 
                ranges = TRUE
            )
            filterReport['matchTo5pAdapter'] <- filterReport['matchTo5pAdapter'] + 
                sum(BiocGenerics::start(ranges) != 1)
            filterReport['matchTo3pAdapter'] <- filterReport['matchTo3pAdapter'] + 
                sum(BiocGenerics::end(ranges) != BiocGenerics::width(chunks))
            chunks <- ShortRead::narrow(x = chunks, 
                                        start = BiocGenerics::start(ranges), 
                                        end = BiocGenerics::end(ranges))
        
            ## filter and write short reads
            filterResults <- S4Vectors::evalSeparately(filters, chunks)
            if (is.list(filterResults)) # return value of evalSeparately?
                filterResults <- do.call(cbind, filterResults)
            #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, 'CleanNFilter'])
            # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
            filterReport['tooManyN'] <- filterReport['tooManyN'] + 
                sum(!filterResults[, grep('CleanNFilter',colnames(filterResults))[1]])
            filterReport['tooShort'] <- filterReport['tooShort'] + sum(!filterResults[, 'LengthFilter'])
            filterReport['lowComplexity'] <- filterReport['lowComplexity'] + 
                sum(!filterResults[, 'LowComplexityFilter'])
            filter <- apply(filterResults, 1, all)
            filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)
            if (sum(filter))
                ShortRead::writeFasta(chunks[filter], tmpOutputFilename, 
                                      mode = mode, compress = FALSE)
            mode <- 'a'
            cycle <- cycle + 1
        }
        
    } else if (fileformat == "fastq") {
        fs1 <- ShortRead::FastqStreamer(filename, n = nrec)
        on.exit(close(fs1))
        on.exit(rm(fs1), add = TRUE)
        mode <- 'w'
        while (length(chunks <- ShortRead::yield(fs1)) != 0L) {
            filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
            
            ## truncate start or end bases
            chunks <- ShortRead::narrow(x = chunks, start = truncateFromBase, 
                                        end = truncateToBase)
            
            ## trim adaptor
            ranges <- Biostrings::trimLRPatterns(
                subject = chunks, Lpattern = Lpattern, Rpattern = Rpattern, 
                max.Lmismatch = max.Lmismatch, max.Rmismatch = max.Rmismatch,
                with.Lindels = with.Lindels, with.Rindels = with.Rindels, 
                ranges = TRUE
            )
            filterReport['matchTo5pAdapter'] <- filterReport['matchTo5pAdapter'] + 
                sum(BiocGenerics::start(ranges) != 1)
            filterReport['matchTo3pAdapter'] <- filterReport['matchTo3pAdapter'] + 
                sum(BiocGenerics::end(ranges) != BiocGenerics::width(chunks))
            chunks <- ShortRead::narrow(x = chunks, 
                                        start = BiocGenerics::start(ranges), 
                                        end = BiocGenerics::end(ranges))
        
            ## filter and write short reads
            filterResults <- S4Vectors::evalSeparately(filters, chunks)
            #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!filterResults[, 'CleanNFilter'])
            # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
            filterReport['tooManyN'] <- filterReport['tooManyN'] + 
                sum(!filterResults[, grep('CleanNFilter',colnames(filterResults))[1]])
            filterReport['tooShort'] <- filterReport['tooShort'] + sum(!filterResults[, 'LengthFilter'])
            filterReport['lowComplexity'] <- filterReport['lowComplexity'] + 
                sum(!filterResults[, 'LowComplexityFilter'])
            filter <- apply(filterResults, 1, all)
            filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)
            if (sum(filter))
                ShortRead::writeFastq(chunks[filter], tmpOutputFilename, 
                                      mode = mode, qualityType = "Auto", compress = FALSE)
            mode <- 'a'
        }
        
    } else {
        stop("unknown file format: ", fileformat)
    }
    
    if(filecompr != "none")
        compressFile(tmpOutputFilename, destname = outputFilename, remove = TRUE)
    
    return(filterReport)
}

#' @keywords internal
#' @importFrom ShortRead readFasta writeFastq FastqStreamer yield narrow
#' @importFrom Biostrings trimLRPatterns
#' @importFrom S4Vectors evalSeparately
preprocessPairedReads <- function(filename, filenameMate, outputFilename, 
                                  outputFilenameMate, fileformat, filecompr,
                                  truncateFromBase, truncateToBase, 
                                  filters,
                                  nrec) {

    message("  filtering ", truncPath(filename, getOption('width') - 17), " and\n    ",
            truncPath(filenameMate, getOption('width') - 16))
    
    ## create (temporary) output file names, will be compressed if filecompr!="none"
    if (filecompr != "none") {
        tmpOutputFilename <- tempfile(pattern = "preprocessReadsTemp")
        tmpOutputFilenameMate <- tempfile(pattern = "preprocessReadsTemp")
    } else {
        tmpOutputFilename <- outputFilename
        tmpOutputFilenameMate <- outputFilenameMate
    }
    
    ## filter chunks
    nrec <- round(nrec/2) # only load half the number of pairs
    filterReport <- c(totalSequences = 0, matchTo5pAdapter = NA, 
                      matchTo3pAdapter = NA, tooShort = 0, tooManyN = 0, 
                      lowComplexity = 0, totalPassed = 0)
    
    if (fileformat == "fasta") {
        mode <- 'w'
        cycle <- 1L
        while(length(chunks <- ShortRead::readFasta(filename, nrec = nrec, 
                                                    skip = (cycle - 1)*nrec)) != 0L) {
            chunksMate <- ShortRead::readFasta(filenameMate, nrec = nrec, 
                                               skip = (cycle - 1)*nrec)
            filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
            
            chunks <- ShortRead::narrow(x = chunks, start = truncateFromBase, 
                                        end = truncateToBase)
            chunksMate <- ShortRead::narrow(x = chunksMate, 
                                            start = truncateFromBase, 
                                            end = truncateToBase)
            
            ## filter and write short reads
            filterResults <- S4Vectors::evalSeparately(filters, chunks)
            if (is.list(filterResults)) # return value of evalSeparately?
                filterResults <- do.call(cbind, filterResults)
            filterResultsMate <- S4Vectors::evalSeparately(filters, chunksMate)
            if (is.list(filterResultsMate)) # return value of evalSeparately?
                filterResults <- do.call(cbind, filterResultsMate)
            #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, 'CleanNFilter'] & filterResultsMate[, 'CleanNFilter']))
            # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
            filterReport['tooManyN'] <- filterReport['tooManyN'] + 
                sum(!(filterResults[, grep('CleanNFilter',colnames(filterResults))[1]] &
                          filterResultsMate[, grep('CleanNFilter', colnames(filterResultsMate))[1]]))
            filterReport['tooShort'] <- filterReport['tooShort'] + 
                sum(!(filterResults[, 'LengthFilter'] &
                          filterResultsMate[, 'LengthFilter']))
            filterReport['lowComplexity'] <- filterReport['lowComplexity'] + 
                sum(!(filterResults[, 'LowComplexityFilter'] &
                          filterResultsMate[, 'LowComplexityFilter']))
            filter <- apply(cbind(filterResults, filterResultsMate), 1, all)
            filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)         
            if (sum(filter)) {
                ShortRead::writeFasta(chunks[filter], tmpOutputFilename, 
                                      mode = mode, compress = FALSE)
                ShortRead::writeFasta(chunksMate[filter], tmpOutputFilenameMate, 
                                      mode = mode, compress = FALSE)
                }
            
            mode <- 'a'
            cycle <- cycle + 1
        }
        
    } else if (fileformat == "fastq") {
        fs1 <- ShortRead::FastqStreamer(filename, n = nrec)
        fs2 <- ShortRead::FastqStreamer(filenameMate, n = nrec)
        on.exit(close(fs1))
        on.exit(close(fs2), add = TRUE)
        on.exit(rm(fs1, fs2), add = TRUE)
        mode <- 'w'
        while (length(chunks <- ShortRead::yield(fs1)) != 0L) {
            chunksMate <- ShortRead::yield(fs2)
            filterReport['totalSequences'] <- filterReport['totalSequences'] + length(chunks)
            
            chunks <- ShortRead::narrow(x = chunks, start = truncateFromBase, 
                                        end = truncateToBase)
            chunksMate <- ShortRead::narrow(x = chunksMate, 
                                            start = truncateFromBase, 
                                            end = truncateToBase)
                
            ## filter and write short reads
            filterResults <- S4Vectors::evalSeparately(filters, chunks)
            filterResultsMate <- S4Vectors::evalSeparately(filters, chunksMate)
            #filterReport['tooManyN'] <- filterReport['tooManyN'] + sum(!(filterResults[, 'CleanNFilter'] & filterResultsMate[, 'CleanNFilter']))
            # workaround for filter being renamed to 'CleanNFilter.other' if nrec=1:
            filterReport['tooManyN'] <- filterReport['tooManyN'] + 
                sum(!(filterResults[, grep('CleanNFilter',colnames(filterResults))[1]] &
                          filterResultsMate[, grep('CleanNFilter',colnames(filterResultsMate))[1]]))
            filterReport['tooShort'] <- filterReport['tooShort'] + 
                sum(!(filterResults[, 'LengthFilter'] &
                          filterResultsMate[, 'LengthFilter']))
            filterReport['lowComplexity'] <- filterReport['lowComplexity'] + 
                sum(!(filterResults[, 'LowComplexityFilter'] &
                          filterResultsMate[, 'LowComplexityFilter']))
            filter <- apply(cbind(filterResults, filterResultsMate), 1, all)
            filterReport['totalPassed'] <- filterReport['totalPassed'] + sum(filter)            
            if (sum(filter)) {
                ShortRead::writeFastq(chunks[filter], tmpOutputFilename, 
                                      mode = mode, qualityType = "Auto", compress = FALSE)
                ShortRead::writeFastq(chunksMate[filter], tmpOutputFilenameMate, 
                                      mode = mode, qualityType = "Auto", compress = FALSE)
            }
            
            mode <- 'a'
        }
        
    } else {
        stop("unknown file format: ", fileformat)
    }
    
    if (filecompr != "none") {
        compressFile(tmpOutputFilename, destname = outputFilename, remove = TRUE)
        compressFile(tmpOutputFilenameMate, destname = outputFilenameMate, remove = TRUE)
    }
    
    return(filterReport)
}

#' @keywords internal
#' @importFrom ShortRead srFilter
#' @importFrom BiocGenerics width
filterLength <- function(threshold = 0L, .name = "LengthFilter") {
    #.check_type_and_length(threshold, "numeric", 1)
    ShortRead::srFilter(function(x) {
        BiocGenerics::width(x) >= threshold
    }, name = .name)
}

#' @keywords internal
#' @importFrom ShortRead srFilter sread
#' @importFrom Biostrings dinucleotideFrequency
filterLowComplexity <- function(threshold = 0.5, referenceEntropy = 3.908135, 
                                .name = "LowComplexityFilter") {
    ShortRead::srFilter(function(x) {
        # less than half the average entropy per dinucleotide of non-random 
        # chromosomes in hg18 (3.908135 bits):
        #   entropy(x)/3.908135 >= threshold
        diNucFreq <- Biostrings::dinucleotideFrequency(ShortRead::sread(x))
        if (is.null(dim(diNucFreq))) {
            diNucFreq <- diNucFreq/sum(diNucFreq)
            H <- -sum(diNucFreq * log2(diNucFreq), na.rm = TRUE)
        } else {
            diNucFreq <- diNucFreq/rowSums(diNucFreq)
            H <- -rowSums(diNucFreq * log2(diNucFreq), na.rm = TRUE)
        }
        H/referenceEntropy >= threshold
    }, name = .name)
}

