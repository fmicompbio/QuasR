#' @keywords internal
#' @importFrom parallel makeCluster stopCluster clusterEvalQ parLapply
#' @importFrom tools file_path_sans_ext
createGenomicAlignments <- function(proj, clObj) {

    # check if a cluster object is provided. if not, use only one core.
    if (is.null(clObj)) {
        clObj <- parallel::makeCluster(1)
        on.exit(parallel::stopCluster(clObj))
    }
    # retrieve information about the cluster
    message("Testing the compute nodes...", appendLF = FALSE)
    tryCatch({
        nodeNamesList <- parallel::clusterEvalQ(clObj, Sys.info()['nodename'])
    }, error = function(ex) {
        message("FAILED")
        stop("The cluster object does not work properly on this system. ",
             "Please consult the manual of the package 'parallel'\n",
             call. = FALSE)
    })
    message("OK")

    message("Loading QuasR on the compute nodes...", appendLF = FALSE)
    # load QuasR package on all the nodes
    loadQuasR(clObj)

    nodeNames <- unlist(nodeNamesList)
    coresPerNode <- table(nodeNames)
    message("Available cores:")
    for (i in seq_along(coresPerNode))
        message(names(coresPerNode)[i], ": ", as.character(coresPerNode[i]))

    # log file that is passed to all the nodes. each node can write into
    # that file and serves as a monitor of their progress
    logFile <- tempfile(tmpdir = getwd(), pattern = "QuasR_log_", fileext = ".txt")

    #create the parameters necessary for the individual processes
    # performing the genomic alignments
    paramsListGenomic = NULL
    for (i in seq_len(nrow(proj@reads))) {
        if (is.na(proj@alignments$FileName)[i]) {
            # create a filename for the current bam file and update the
            # qProject. ensure uniqueness by adding a random suffix
            if (is.na(proj@alignmentsDir)) {
                bamDir <- dirname(proj@reads[i, 1])
            } else {
                bamDir <- proj@alignmentsDir
            }
            samplePrefix <- basename(tools::file_path_sans_ext(proj@reads[i, 1],
                                                               compression = TRUE))
            proj@alignments$FileName[i] <- tempfile(tmpdir = bamDir,
                                                    pattern = paste(samplePrefix,
                                                                    "_", sep = ""),
                                                    fileext = ".bam")

            paramsListGenomic[[length(paramsListGenomic) + 1]] <- list(
                "sampleNr" = i,
                "qProject" = proj,
                "coresPerNode" = coresPerNode,
                "logFile" = logFile
            )
        }
    }

    # perform all necessary genomic alignments
    n_tasks <- length(paramsListGenomic)
    if (n_tasks > 0) {
        message(paste("Performing genomic alignments for", n_tasks,
                      "samples. See progress in the log file:"))
        message(logFile)

        # create a cluster object that has only one entry per machine.
        # multithreading (on a lower level) results
        # in a fully occupied node. This new clObj does not have to be stopped.
        # it closes once the original one closes
        clObjNR <- clObj[!duplicated(nodeNames)]
        
        parallel::clusterExport(clObjNR, "n_tasks", envir=environment())
        parallel::parLapply(clObjNR, paramsListGenomic,
                            createGenomicAlignmentsController)
        message("Genomic alignments have been created successfully")
        message("")
    }

    return(proj)
}

#' @keywords internal
#' @importFrom parallel makeCluster stopCluster clusterEvalQ parLapply
#' @importFrom tools file_path_sans_ext
createAuxAlignments <- function(proj, clObj) {

    # check if a cluster object is provided. if not, use only one core.
    if (is.null(clObj)) {
        clObj <- parallel::makeCluster(1)
        on.exit(parallel::stopCluster(clObj))
    }
    # retrieve information about the cluster
    message("Testing the compute nodes...", appendLF = FALSE)
    tryCatch({
        nodeNamesList <- parallel::clusterEvalQ(clObj, Sys.info()['nodename'])
    }, error = function(ex) {
        message("FAILED")
        stop("The cluster object does not work properly on this system. ",
             "Please consult the manual of the package 'parallel'\n",
             call. = FALSE)
    })
    message("OK")

    message("Loading QuasR on the compute nodes...", appendLF = FALSE)
    # load QuasR package on all the nodes
    loadQuasR(clObj)

    nodeNames <- unlist(nodeNamesList)
    coresPerNode <- table(nodeNames)
    message("Available cores:")
    print(coresPerNode)

    # log file that is passed to all the nodes. each node can write into
    # that file and serves as a monitor of their progress
    logFile <- tempfile(tmpdir = getwd(), pattern = "QuasR_log_", fileext = ".txt")

    # create the parameters necessary for the individual processes
    # performing the aux alignments
    paramsListAux = NULL
    for (i in seq_len(ncol(proj@auxAlignments))) {
        if (any(is.na(proj@auxAlignments[, i]))) {
            if (is.na(proj@alignmentsDir)) {
                bamDir <- dirname(proj@reads[i, 1])
            } else {
                bamDir <- proj@alignmentsDir
            }
            samplePrefix <- basename(tools::file_path_sans_ext(proj@reads[i, 1],
                                                               compression = TRUE))
            auxNrs <- NULL # the aux files to map in the children
            for (j in seq_len(nrow(proj@auxAlignments))) {
                if (is.na(proj@auxAlignments[j,i])) {
                    auxNrs <- c(auxNrs, j)
                    proj@auxAlignments[j, i] <- tempfile(
                        tmpdir = bamDir,
                        pattern = paste(samplePrefix, "_", sep = ""),
                        fileext = ".bam"
                    )
                }
            }
            paramsListAux[[length(paramsListAux) + 1]] <- list(
                "sampleNr" = i,
                "auxNrs" = auxNrs,
                "qProject" = proj,
                "coresPerNode" = coresPerNode,
                "logFile" = logFile
            )
        }
    }

    n_tasks <- length(paramsListAux)
    if (n_tasks > 0) {
        message(paste("Performing auxiliary alignments for",
                      n_tasks, "samples. See progress in the log file:"))
        message(logFile)

        # create a cluster object that has only one entry per machine.
        # multithreading (on a lower level) results
        # in a fully occupied node. This new clObj does not have to be stopped.
        # it closes once the original one closes
        clObjNR <- clObj[!duplicated(nodeNames)]
        
        parallel::clusterExport(clObjNR, "n_tasks", envir=environment())
        parallel::parLapply(clObjNR, paramsListAux, createAuxAlignmentsController)
        message("Auxiliary alignments have been created successfully")
        message("")
    }

    return(proj)
}

#' @keywords internal
#' @importFrom tools file_path_sans_ext
#' @importFrom utils installed.packages write.table
#' @importFrom stats na.omit
#' @importFrom Rsamtools asBam
createGenomicAlignmentsController <- function(params) {
    tryCatch({# try catch block goes through the whole function

        # extract the parameters from params
        sampleNr <- params$sampleNr
        proj <- params$qProject
        coresPerNode <- params$coresPerNode
        logFile <- params$logFile
        sink(logFile, append = TRUE) # redirect all the print output to a logfile
        cacheDir <- resolveCacheDir(proj@cacheDir) # tmp dir can change for each machine

        # find out how many threads are available on this node
        # (for running the alignments)
        thisNodesName <- Sys.info()['nodename']
        if (!(thisNodesName %in% names(coresPerNode))) {
            stop("Fatal error 2394793")
        }
        coresThisNode <- coresPerNode[names(coresPerNode) %in% thisNodesName]

        # try to load all the required libraries on the compute node.
        # these are the aligner package
        # and in the case of a BSgenome, the BSgenome itself needs to be
        # loaded and the associated index
        if (!require(proj@aligner, character.only = TRUE, quietly = TRUE)) {
            stop("Could not load the aligner package ", proj@aligner)
        }

        if (proj@genomeFormat == "file") {
            indexDir <- paste(proj@genome, proj@alnModeID, sep = ".")
        } else if (proj@genomeFormat == "BSgenome") {
            if (!(proj@genome %in% utils::installed.packages()[, 'Package'])) {
                stop("The genome package ", proj@genome, " is not installed")
            }
            if (is.na(proj@snpFile)) {
                if (!(paste(proj@genome, proj@alnModeID, sep = ".") %in%
                      utils::installed.packages()[, 'Package'])) {
                    stop("The genome index package ", paste(proj@genome,
                                                            proj@alnModeID, sep = "."),
                         " is not installed")
                }
            }

            indexDir <- system.file("alignmentIndex", package = paste(proj@genome,
                                                                      proj@alnModeID,
                                                                      sep = "."))
        } else {
            stop("Fatal error 4778493")
        }

        # create the info file for the bam file
        bamInfo <- qProjectBamInfo(proj, sampleNr)

        # uncompress the files containing the reads (twice for paired end samples)
        if (proj@paired == "no") {
            # SINGLE END SAMPLES
            # decompress the reads if necessary
            reads <- proj@reads$FileName[sampleNr]
            if (compressedFileFormat(reads) != "none") {
                worker_message("Decompressing SE read file on [", Sys.info()['nodename'], "]:", reads)
                readsUNC <- tempfile(tmpdir = cacheDir, pattern = basename(reads),
                                     fileext = paste(".", proj@samplesFormat, sep = ""))
                compressFile(reads, readsUNC, remove = FALSE)
                proj@reads$FileName[sampleNr] <- readsUNC
                on.exit(file.remove(readsUNC)) # make sure that the temp file is deleted
            }
        } else {
            # PAIRED END SAMPLES
            # decompress the first read pair if necessary
            reads1 <- proj@reads$FileName1[sampleNr]
            if (compressedFileFormat(reads1) != "none") {
                worker_message("Decompressing PE R1 file on [", Sys.info()['nodename'], "]:", reads1)
                readsUNC1 <- tempfile(tmpdir = cacheDir, pattern = basename(reads1),
                                      fileext = paste(".", proj@samplesFormat, sep = ""))
                compressFile(reads1, readsUNC1, remove = FALSE)
                proj@reads$FileName1[sampleNr] <- readsUNC1
                # make sure that the temp file is deleted at the end
                on.exit(file.remove(readsUNC1))
            }
            # decompress the second read pair if necessary
            reads2 <- proj@reads$FileName2[sampleNr]
            if (compressedFileFormat(reads2) != "none") {
                worker_message("Decompressing PE R2 file on [", Sys.info()['nodename'], "]:", reads2)
                readsUNC2 <- tempfile(tmpdir = cacheDir, pattern = basename(reads2),
                                      fileext = paste(".", proj@samplesFormat, sep = ""))
                compressFile(reads2, readsUNC2, remove = FALSE)
                proj@reads$FileName2[sampleNr] <- readsUNC2
                # make sure that the temp file is deleted at the end
                on.exit(file.remove(readsUNC2), add = TRUE)
            }
        }

        samFile <- tempfile(tmpdir = cacheDir,
                            pattern = basename(proj@reads[sampleNr, 1]),
                            fileext = ".sam")
        # make sure that the temp file is deleted at the end
        on.exit(file.remove(samFile), add = TRUE)

        # perform the required alignments based on the information in qProject
        if (proj@alnModeID == "Rbowtie") {
            if (!proj@splicedAlignment) {
                if (is.na(proj@snpFile)) {
                    align_Rbowtie(indexDir, proj@reads[sampleNr, ],
                                  proj@samplesFormat, proj@paired,
                                  proj@alignmentParameter, coresThisNode,
                                  samFile, cacheDir)
                } else {
                    # add numeric id to the reads, this is required for the correct
                    # operation of mergeReorderSam in allelic mode
                    proj@reads[sampleNr, ] <- addNumericToID(proj@reads[sampleNr, ],
                                                             proj@paired, cacheDir)
                    # make sure that the temp file(s) are deleted at the end
                    on.exit(file.remove(unlist(
                        proj@reads[sampleNr, stats::na.omit(match(c("FileName",
                                                                    "FileName1",
                                                                    "FileName2"),
                                                                  colnames(proj@reads)))]
                    )), add = TRUE)

                    samFileR <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    samFileA <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    on.exit(file.remove(samFileR),add = TRUE)
                    on.exit(file.remove(samFileA),add = TRUE)
                    align_Rbowtie(paste(proj@snpFile, basename(proj@genome),
                                        "R", "fa", proj@alnModeID, sep = "."),
                                  proj@reads[sampleNr, ], proj@samplesFormat,
                                  proj@paired, proj@alignmentParameter,
                                  coresThisNode, samFileR, cacheDir)
                    align_Rbowtie(paste(proj@snpFile, basename(proj@genome),
                                        "A", "fa", proj@alnModeID, sep = "."),
                                  proj@reads[sampleNr, ], proj@samplesFormat,
                                  proj@paired, proj@alignmentParameter,
                                  coresThisNode, samFileA, cacheDir)
                    mrQuSize <- .Call(mergeReorderSam, c(samFileR, samFileA),
                                      samFile, as.integer(2), as.integer(proj@maxHits))
                    print(paste("mergeReorderMaxQueueSize", mrQuSize))
                }
            } else {
                if (is.na(proj@snpFile)) {
                    # in the case of BSgenome, flush it because it is needed by SpliceMap
                    if (proj@genomeFormat == "BSgenome") {
                        if (!require(paste(proj@genome, proj@alnModeID, sep = "."),
                                     character.only = TRUE, quietly = TRUE)) {
                            stop("Could not load the genome index package ",
                                 paste(proj@genome, proj@alnModeID, sep = "."))
                        }
                        genomeObj <- get(proj@genome) # access the BSgenome
                        fastaFilepath <- tempfile(tmpdir = cacheDir, fileext = ".fa")
                        print(paste("Writing BSgenome to disk on",
                                    Sys.info()['nodename'], ":", fastaFilepath))
                        # flush the BSgenome to disk
                        BSgenomeSeqToFasta(genomeObj, fastaFilepath)
                        on.exit(unlink(fastaFilepath), add = TRUE)
                    } else {
                        fastaFilepath <- proj@genome
                    }
                    align_RbowtieSpliced(fastaFilepath, indexDir, proj@reads[sampleNr, ],
                                         proj@samplesFormat, proj@paired,
                                         proj@alignmentParameter, coresThisNode,
                                         samFile, cacheDir)
                } else {
                    # add numeric id to the reads, this is required for the correct
                    # operation of mergeReorderSam in allelic mode
                    proj@reads[sampleNr, ] <- addNumericToID(proj@reads[sampleNr, ],
                                                             proj@paired, cacheDir)
                    # make sure that the temp file(s) are deleted at the end
                    on.exit(file.remove(unlist(
                        proj@reads[sampleNr, stats::na.omit(match(c("FileName",
                                                                    "FileName1",
                                                                    "FileName2"),
                                                                  colnames(proj@reads)))]
                    )), add = TRUE)

                    samFileR <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    samFileA <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    on.exit(file.remove(samFileR), add = TRUE)
                    on.exit(file.remove(samFileA), add = TRUE)
                    align_RbowtieSpliced(paste(proj@snpFile, basename(proj@genome),
                                               "R", "fa", sep = "."),
                                         paste(proj@snpFile, basename(proj@genome),
                                               "R", "fa", proj@alnModeID,
                                               sep = "."),
                                         proj@reads[sampleNr, ], proj@samplesFormat,
                                         proj@paired, proj@alignmentParameter,
                                         coresThisNode, samFileR, cacheDir)
                    align_RbowtieSpliced(paste(proj@snpFile, basename(proj@genome),
                                               "A", "fa", sep = "."),
                                         paste(proj@snpFile, basename(proj@genome),
                                               "A", "fa", proj@alnModeID, sep = "."),
                                         proj@reads[sampleNr, ], proj@samplesFormat,
                                         proj@paired, proj@alignmentParameter,
                                         coresThisNode, samFileA, cacheDir)
                    mrQuSize <- .Call(mergeReorderSam, c(samFileR, samFileA),
                                      samFile, as.integer(2), as.integer(proj@maxHits))
                    print(paste("mergeReorderMaxQueueSize", mrQuSize))
                }
            }
        } else if (proj@alnModeID == "RbowtieCtoT") {
            if (proj@bisulfite == "dir") {
                if (is.na(proj@snpFile)) {
                    align_RbowtieCtoT_dir(indexDir, proj@reads[sampleNr, ],
                                          proj@samplesFormat, proj@paired,
                                          proj@alignmentParameter, !is.na(proj@snpFile),
                                          proj@maxHits, coresThisNode, samFile, cacheDir)
                } else {
                    samFileR <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    samFileA <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    on.exit(file.remove(samFileR), add = TRUE)
                    on.exit(file.remove(samFileA), add = TRUE)

                    align_RbowtieCtoT_dir(paste(proj@snpFile, basename(proj@genome),
                                                "R", "fa", proj@alnModeID, sep = "."),
                                          proj@reads[sampleNr, ], proj@samplesFormat,
                                          proj@paired, proj@alignmentParameter,
                                          !is.na(proj@snpFile), proj@maxHits,
                                          coresThisNode, samFileR, cacheDir)
                    align_RbowtieCtoT_dir(paste(proj@snpFile, basename(proj@genome),
                                                "A", "fa", proj@alnModeID, sep = "."),
                                          proj@reads[sampleNr, ], proj@samplesFormat,
                                          proj@paired, proj@alignmentParameter,
                                          !is.na(proj@snpFile), proj@maxHits,
                                          coresThisNode, samFileA, cacheDir)
                    mrQuSize <- .Call(mergeReorderSam, c(samFileR, samFileA),
                                      samFile, as.integer(2), as.integer(proj@maxHits))
                    print(paste("mergeReorderMaxQueueSize", mrQuSize))
                }
            } else {
                if (is.na(proj@snpFile)) {
                    align_RbowtieCtoT_undir(indexDir, proj@reads[sampleNr, ],
                                            proj@samplesFormat, proj@paired,
                                            proj@alignmentParameter,
                                            !is.na(proj@snpFile),
                                            proj@maxHits, coresThisNode, samFile, cacheDir)
                } else {
                    samFileR <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    samFileA <- tempfile(tmpdir = cacheDir,
                                         pattern = basename(proj@reads[sampleNr, 1]),
                                         fileext = ".sam")
                    on.exit(file.remove(samFileR), add = TRUE)
                    on.exit(file.remove(samFileA), add = TRUE)

                    align_RbowtieCtoT_undir(paste(proj@snpFile, basename(proj@genome),
                                                  "R", "fa", proj@alnModeID, sep = "."),
                                            proj@reads[sampleNr, ], proj@samplesFormat,
                                            proj@paired, proj@alignmentParameter,
                                            !is.na(proj@snpFile), proj@maxHits,
                                            coresThisNode, samFileR, cacheDir)
                    align_RbowtieCtoT_undir(paste(proj@snpFile, basename(proj@genome),
                                                  "A", "fa", proj@alnModeID, sep = "."),
                                            proj@reads[sampleNr, ], proj@samplesFormat,
                                            proj@paired, proj@alignmentParameter,
                                            !is.na(proj@snpFile), proj@maxHits,
                                            coresThisNode, samFileA, cacheDir)
                    mrQuSize <- .Call(mergeReorderSam, c(samFileR, samFileA),
                                      samFile, as.integer(2), as.integer(proj@maxHits))
                    print(paste("mergeReorderMaxQueueSize", mrQuSize))
                }
            }
        } else if (proj@alnModeID == "Rhisat2") {
            if (is.na(proj@snpFile)) {
                align_Rhisat2(indexDir, proj@reads[sampleNr, ],
                              proj@samplesFormat, proj@paired, proj@alignmentParameter,
                              coresThisNode, samFile, cacheDir,
                              proj@splicedAlignment, proj@maxHits)
            } else {
                # add numeric id to the reads, this is required for the correct
                # operation of mergeReorderSam in allelic mode
                proj@reads[sampleNr,] <- addNumericToID(proj@reads[sampleNr, ],
                                                        proj@paired, cacheDir)
                # make sure that the temp file(s) are deleted at the end
                on.exit(file.remove(unlist(
                    proj@reads[sampleNr, stats::na.omit(match(c("FileName", "FileName1",
                                                                "FileName2"),
                                                              colnames(proj@reads)))]
                )), add = TRUE)
                samFileR <- tempfile(tmpdir = cacheDir,
                                     pattern = basename(proj@reads[sampleNr, 1]),
                                     fileext = ".sam")
                samFileA <- tempfile(tmpdir = cacheDir,
                                     pattern = basename(proj@reads[sampleNr, 1]),
                                     fileext = ".sam")
                on.exit(file.remove(samFileR), add = TRUE)
                on.exit(file.remove(samFileA), add = TRUE)
                align_Rhisat2(paste(proj@snpFile, basename(proj@genome), "R",
                                    "fa", proj@alnModeID, sep = "."),
                              proj@reads[sampleNr, ], proj@samplesFormat,
                              proj@paired, proj@alignmentParameter,
                              coresThisNode, samFileR, cacheDir,
                              proj@splicedAlignment, proj@maxHits)
                align_Rhisat2(paste(proj@snpFile, basename(proj@genome), "A",
                                    "fa", proj@alnModeID, sep = "."),
                              proj@reads[sampleNr, ], proj@samplesFormat,
                              proj@paired, proj@alignmentParameter,
                              coresThisNode, samFileA, cacheDir,
                              proj@splicedAlignment, proj@maxHits)
                mrQuSize <- .Call(mergeReorderSam, c(samFileR, samFileA),
                                  samFile, as.integer(2), as.integer(proj@maxHits))
                print(paste("mergeReorderMaxQueueSize", mrQuSize))
            }
        } else {
            stop("Fatal error 23484303")
        }

        worker_message(
          "Converting sam file to sorted bam file on [",
          Sys.info()['nodename'], "] :", samFile)
        if (coresThisNode > 3) {
            # sort sam and convert to bam parallel
            samToSortedBamParallel(
                samFile, tools::file_path_sans_ext(proj@alignments$FileName[sampleNr]),
                coresThisNode, cacheDir
            )
        } else {
            # sort sam and convert to bam
            Rsamtools::asBam(
                samFile,
                tools::file_path_sans_ext(proj@alignments$FileName[sampleNr])
            )
        }

        # save the info file for the bam file
        utils::write.table(bamInfo, paste(proj@alignments$FileName[sampleNr],
                                          "txt", sep = "."),
                           sep = "\t", quote = FALSE, col.names = FALSE)

        worker_message(
          "Genomic alignments for sample ", sampleNr, " (",
          proj@reads$SampleName[sampleNr],
          ") have been successfully created on [",
          Sys.info()['nodename'], "]")

        # if one process stops due to an error, catch it, concatenate the
        # message with specific information about
        # the compute node and then print it to the log file. this procedure
        # provides the information about which
        # sample was processed and on which machine when the error occured.
    }, error = function(ex) {
        emsg <- paste("Error on", Sys.info()['nodename'],
                      "processing sample", proj@reads[sampleNr, 1], ":", ex$message)
        print(emsg)
        stop(emsg)
    }, finally = {
        sink() # close the redirection of the print statements
    })
}

#' @keywords internal
#' @importFrom tools file_path_sans_ext
#' @importFrom Rsamtools sortBam indexBam
#' @importFrom utils write.table
createAuxAlignmentsController <- function(params) {
    tryCatch({ # try catch block goes through the whole function

        # extract the parameters from params
        sampleNr <- params$sampleNr
        auxNrs <- params$auxNrs
        proj <- params$qProject
        coresPerNode <- params$coresPerNode
        logFile <- params$logFile
        sink(logFile, append = TRUE) # redirect all the print output to a logfile
        cacheDir <- resolveCacheDir(proj@cacheDir) # tmp dir can change for each machine

        # load the aligner package on the compute node
        if (!require(proj@aligner, character.only = TRUE, quietly = TRUE)) {
            stop("Could not load the aligner package ", proj@aligner)
        }

        # find out how many threads are available on this node
        # (for running the alignments)
        thisNodesName <- Sys.info()['nodename']
        if (!(thisNodesName %in% names(coresPerNode))) {
            stop("Fatal error 23594793")
        }
        coresThisNode <- coresPerNode[names(coresPerNode) %in% thisNodesName]

        # extract the unmapped reads from bam file. in the case of a bisulfite
        # sample, this requires the conversion from ff (that is present in the bam file)
        # to fr to make the rest of the execution compatible (see last parameter
        # of extractUnmappedReads)
        unmappedReadsInfo <- proj@reads[sampleNr, ]
        if (proj@paired == "no") {
            unmappedReadsFile <- tempfile(
                tmpdir = cacheDir,
                pattern = basename(proj@alignments$FileName[sampleNr]),
                fileext = proj@samplesFormat
            )
            .Call(extractUnmappedReads, proj@alignments$FileName[sampleNr],
                  unmappedReadsFile, proj@samplesFormat == "fastq",
                  proj@bisulfite != "no")
            unmappedReadsInfo$FileName <- unmappedReadsFile
            # make sure that the temp file is deleted
            on.exit(file.remove(unmappedReadsFile))
        } else {
            unmappedReadsFile1 <- tempfile(
                tmpdir = cacheDir,
                pattern = basename(proj@alignments$FileName[sampleNr]),
                fileext = proj@samplesFormat
            )
            unmappedReadsFile2 <- tempfile(
                tmpdir = cacheDir,
                pattern = basename(proj@alignments$FileName[sampleNr]),
                fileext = proj@samplesFormat
            )
            .Call(extractUnmappedReads, proj@alignments$FileName[sampleNr],
                  c(unmappedReadsFile1, unmappedReadsFile2),
                  proj@samplesFormat == "fastq", proj@bisulfite != "no")
            unmappedReadsInfo$FileName1 <- unmappedReadsFile1
            unmappedReadsInfo$FileName2 <- unmappedReadsFile2
            # make sure that the temp files are deleted
            on.exit(file.remove(unmappedReadsFile1))
            on.exit(file.remove(unmappedReadsFile2), add = TRUE)
        }

        # for the current sample, create alignments against all the aux files
        for (j in auxNrs) {
            samFile <- tempfile(tmpdir = cacheDir,
                                pattern = basename(proj@reads[sampleNr, 1]),
                                fileext = ".sam")
            bamFileNoUnmapped <- tempfile(tmpdir = cacheDir,
                                          pattern = basename(proj@reads[sampleNr, 1]),
                                          fileext = ".bam")

            if (proj@alnModeID == "Rbowtie") {
                if (!proj@splicedAlignment) {
                    align_Rbowtie(paste(proj@aux$FileName[j], proj@alnModeID, sep = "."),
                                  unmappedReadsInfo, proj@samplesFormat,
                                  proj@paired, proj@alignmentParameter,
                                  coresThisNode, samFile, cacheDir)
                } else {
                    align_RbowtieSpliced(proj@aux$FileName[j],
                                         paste(proj@aux$FileName[j], proj@alnModeID,
                                               sep = "."),
                                         unmappedReadsInfo, proj@samplesFormat,
                                         proj@paired, proj@alignmentParameter,
                                         coresThisNode, samFile, cacheDir)
                }
            } else if (proj@alnModeID == "RbowtieCtoT") {
                if (proj@bisulfite == "dir") {
                    align_RbowtieCtoT_dir(paste(proj@aux$FileName[j], proj@alnModeID,
                                                sep = "."),
                                          unmappedReadsInfo, proj@samplesFormat,
                                          proj@paired, proj@alignmentParameter,
                                          FALSE, proj@maxHits,
                                          coresThisNode, samFile, cacheDir)
                } else {
                    align_RbowtieCtoT_undir(paste(proj@aux$FileName[j], proj@alnModeID,
                                                  sep = "."),
                                            unmappedReadsInfo, proj@samplesFormat,
                                            proj@paired, proj@alignmentParameter,
                                            FALSE, proj@maxHits,
                                            coresThisNode, samFile, cacheDir)
                }
            } else if (proj@alnModeID == "Rhisat2") {
                align_Rhisat2(paste(proj@aux$FileName[j], proj@alnModeID, sep = "."),
                              unmappedReadsInfo, proj@samplesFormat, proj@paired,
                              proj@alignmentParameter, coresThisNode, samFile,
                              cacheDir, proj@splicedAlignment, proj@maxHits)
            } else {
                stop("Fatal error 23484303")
            }

            # remove the unmapped reads and convert to sorted bam
            worker_message(
              "Converting sam file to sorted bam file on [",
              Sys.info()['nodename'], "] :", samFile)
            .Call(removeUnmappedFromSamAndConvertToBam, samFile, bamFileNoUnmapped)
            file.remove(samFile)
            # sort bam
            Rsamtools::sortBam(bamFileNoUnmapped,
                               tools::file_path_sans_ext(proj@auxAlignments[j, sampleNr]))
            Rsamtools::indexBam(proj@auxAlignments[j, sampleNr])
            file.remove(bamFileNoUnmapped)

            # create the info file for the bam file
            bamInfo <- qProjectBamInfo(proj, sampleNr, j)
            utils::write.table(bamInfo,
                               paste(proj@auxAlignments[j, sampleNr],
                                     "txt", sep = "."),
                               sep = "\t", quote = FALSE, col.names = FALSE)

            worker_message(
              "Auxiliary alignments ", j, " for sample ", sampleNr,
              " (", proj@reads$SampleName[sampleNr],
              ") have been successfully created on [",
              Sys.info()['nodename'])

  }

        # if one process stops due to an error, catch it, concatenate the
        # message with specific information about
        # the compute node and then print it to the log file. this procedure
        # provides the information about which
        # sample was processed and on which machine when the error occured.
    }, error = function(ex) {
        emsg <- paste("Error on", Sys.info()['nodename'],"processing sample",
                      proj@reads[sampleNr, 1], ":", ex$message)
        print(emsg)
        stop(emsg)
    }, finally = {
        sink() # close the redirection of the print statements
    })
}

#' @keywords internal
#' @import Rbowtie
align_Rbowtie <- function(indexDir, reads, samplesFormat, paired,
                          alignmentParameter, threads, outFile, cacheDir) {

    # add some variable parameters based on the input format
    if (samplesFormat == "fasta") {
        alignmentParameterAdded <- "-f"
    } else {
        alignmentParameterAdded <- paste("--phred", reads$phred, "-quals", sep = "")
    }

    worker_message(
      "Executing bowtie on [", Sys.info()['nodename'],"] ",
      "using ", threads, " cores. Parameters:")
    if (paired == "no") {
        args <- paste(shQuote(file.path(indexDir, "bowtieIndex")),
                      shQuote(reads$FileName), alignmentParameter,
                      alignmentParameterAdded, "-S", "-p", threads, shQuote(outFile))
        worker_message(args)
        ret <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                       args, stdout = TRUE, stderr = TRUE)
    } else {
        args <- paste(shQuote(file.path(indexDir, "bowtieIndex")),
                      "-1", shQuote(reads$FileName1),
                      "-2", shQuote(reads$FileName2),
                      paste("--", paired, sep = ""), alignmentParameter,
                      alignmentParameterAdded, "-S", "-p", threads, shQuote(outFile))
        worker_message(args)
        ret <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                       args, stdout = TRUE, stderr = TRUE)
    }
    if (!(grepl(" alignments", ret[length(ret)]))) {
        stop("bowtie failed to perform the alignments")
    }
}

#' @keywords internal
align_Rhisat2 <- function(indexDir, reads, samplesFormat, paired,
                          alignmentParameter, threads, outFile, cacheDir,
                          splicedAlignment, maxHits) {
    # add some variable parameters based on the input format
    if (samplesFormat == "fasta") {
        alignmentParameterAdded <- "-f"
    } else {
        alignmentParameterAdded <- paste("--phred", reads$phred, sep = "")
    }
    worker_message(
      "Executing hisat2 on [", Sys.info()['nodename'], "] ",
      "using ", threads, " cores. Parameters:")
    if (paired == "no") {
        args <- paste(shQuote(file.path(indexDir, "hisat2Index")),
                      shQuote(reads$FileName), alignmentParameter,
                      alignmentParameterAdded, "-p", threads,
                      "-S", shQuote(paste(outFile, "tmp", sep = ".")))
    } else {
        args <- paste(shQuote(file.path(indexDir, "hisat2Index")),
                      "-1", shQuote(reads$FileName1), "-2", shQuote(reads$FileName2),
                      paste0("--", paired), alignmentParameter, alignmentParameterAdded,
                      "-p", threads, "-S", shQuote(paste(outFile, "tmp", sep = ".")))
    }
    if (!splicedAlignment) {
        args <- paste(args, "--no-spliced-alignment")
    }
    worker_message(args)
    ret <- system2(file.path(system.file(package = "Rhisat2"), "hisat2"),
                   args, stdout = TRUE, stderr = TRUE)
    if (!(grepl(" reads", ret[1]))) {
        stop("hisat2 failed to perform the alignments")
    }

    ## Filter alignments (keep only reads with at most maxHits alignments, only 1
    ## hit for multimapping reads)
    fhs <- .Call(filterHisat2, paste(outFile, "tmp", sep = "."),
                 outFile, as.integer(maxHits))
    print(paste("Number of filtered secondary alignments:", fhs["n_secondary"]))
    print(paste("Number of filtered overmapped alignments:", fhs["n_overmapped"]))
}

#' @keywords internal
#' @importFrom Rbowtie SpliceMap
#' @importFrom Rsamtools scanFaIndex
align_RbowtieSpliced <- function(genomeFilepath, indexDir, reads, samplesFormat,
                                 paired, alignmentParameter, threads,
                                 outFile, cacheDir) {
    # determine the number of chromosomes (or sequences in case of aux). this is needed to ensure that no more threads
    # than chromosomes are used in the spliced alignment step of SpliceMap
    nrChr <- length(Rsamtools::scanFaIndex(genomeFilepath))
    nrChrTogether <- min(threads, nrChr)

    # parse alignment parameters
    keyValueMatrix <- matrix(unlist(strsplit(alignmentParameter,
                                             "\\s+", perl = TRUE)),
                             ncol = 2, byrow = TRUE)
    if (!all(substr(keyValueMatrix[, 1], 1, 1) == "-")) {
        stop(paste("The parameters for spliced-aligment are in an incorrect format:",
                   alignmentParameter))
    }

    rownames(keyValueMatrix) <- substring(keyValueMatrix[, 1], 2)
    sm_cfg <- as.list(keyValueMatrix[, 2])

    # convert integer strings to actual integers
    all.numbers <- grep("^[[:digit:]]*$", keyValueMatrix[, 2])
    sm_cfg[all.numbers] <- lapply(sm_cfg[all.numbers], as.integer)

    # convert boolean to actual logicals
    all.logicals <- grep("^(TRUE|FALSE)$", keyValueMatrix[, 2])
    sm_cfg[all.logicals] <- lapply(sm_cfg[all.logicals], as.logical)

    sm_cfg[["genome_dir"]] <- genomeFilepath
    if (paired == "no") {
        sm_cfg[["reads_list1"]] <- reads$FileName
    } else {
        sm_cfg[["reads_list1"]] <- reads$FileName1
        sm_cfg[["reads_list2"]] <- reads$FileName2
    }
    sm_cfg[["read_format"]] <- toupper(samplesFormat)
    if (samplesFormat == "fastq") {
        sm_cfg[["quality_format"]] <- paste("phred", reads$phred, sep = "-")
    }
    sm_cfg[["temp_path"]] <- cacheDir
    sm_cfg[["num_chromosome_together"]] <- nrChrTogether
    sm_cfg[["num_threads"]] <- threads
    sm_cfg[["bowtie_base_dir"]] <- file.path(indexDir, "bowtieIndex")
    sm_cfg[["outfile"]] <- outFile

    sm_cfgParString <- cbind(paste("-", as.character(names(sm_cfg)), sep = ""),
                             as.character(sm_cfg))
    print(paste("Executing Splicemap on", Sys.info()['nodename'],
                "using", threads, "cores. Parameters:"))
    print(paste(paste(sm_cfgParString[, 1], sm_cfgParString[, 2]), collapse = " "))

    Rbowtie::SpliceMap(sm_cfg)
}

#' @keywords internal
#' @import Rbowtie
align_RbowtieCtoT_dir <- function(indexDir, reads, samplesFormat, paired,
                                  alignmentParameter, allelic, maxHits,
                                  threads, outFile, cacheDir) {
    # add some variable parameters based on the input format
    if (samplesFormat == "fasta") {
        alignmentParameterAdded <- "-f"
    } else {
        alignmentParameterAdded <- paste("--phred", reads$phred, "-quals", sep = "")
    }

    # set the format of the read identifier. this is necessary for mergeReorderSam outside of this function in allelic mode
    if (!allelic) {
        idMode <- 1
    } else {
        idMode <- 3
    }

    print(paste("Executing bowtie (CtoT and GtoA) on", Sys.info()['nodename'],
                "using", threads, "cores. Parameters:"))
    if (paired == "no") {
        # CtoT convert the reads. include the original sequence in the identifier.
        readsCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName), "readsCtoT_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        .Call(convertReadsIdBisRc, reads$FileName, readsCtoT, c("C", "T"), FALSE)
        # make sure that the temp file is deleted
        on.exit(file.remove(readsCtoT), add = TRUE)

        # create temp filenames for the aligment results
        readsCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName), "readsCtoT_genomeCtoT_", sep = "_"),
            fileext = ".sam"
        )
        readsCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName), "readsCtoT_genomeGtoA_", sep = "_"),
            fileext = ".sam"
        )
        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsCtoT_genomeGtoA), add = TRUE)

        argsCtoT <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                          shQuote(readsCtoT), alignmentParameter,
                          alignmentParameterAdded, "-S", "-p", threads,
                          "--norc", shQuote(readsCtoT_genomeCtoT))
        argsGtoA <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                          shQuote(readsCtoT), alignmentParameter,
                          alignmentParameterAdded, "-S", "-p", threads,
                          "--nofw", shQuote(readsCtoT_genomeGtoA))
        print(argsCtoT)
        print(argsGtoA)
        # perform two alignments, readsCtoT agains genomeCtoT (plus strand)
        # and readsCtoT agains genomeGtoA (minus strand)
        ret1 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        argsCtoT, stdout = TRUE, stderr = TRUE)
        ret2 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        argsGtoA, stdout = TRUE, stderr = TRUE)

        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (GtoA) failed to perform the alignments")
        }

        mrQuSize <- .Call(mergeReorderSam, c(readsCtoT_genomeCtoT,
                                             readsCtoT_genomeGtoA),
                          outFile, as.integer(idMode), as.integer(maxHits))
        print(paste("mergeReorderMaxQueueSize", mrQuSize))

    } else if (paired == "fr") {
        # CtoT convert the reads. include the original sequence in the identifier.
        readsCtoT_1 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsCtoT_1_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        readsCtoT_2 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsCtoT_2_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        .Call(convertReadsIdBisRc, reads$FileName1, readsCtoT_1, c("C", "T"), FALSE)
        .Call(convertReadsIdBisRc, reads$FileName2, readsCtoT_2, c("C", "T"), TRUE) # reverse complement the second read because its fr
        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_1), add = TRUE)
        on.exit(file.remove(readsCtoT_2), add = TRUE)

        # create temp filenames for the aligment results
        readsCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsCtoT_genomeCtoT_", sep = "."),
            fileext = ".sam"
        )
        readsCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsCtoT_genomeGtoA_", sep = "."),
            fileext = ".sam"
        )
        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsCtoT_genomeGtoA), add = TRUE)

        argsCtoT <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                          "-1", shQuote(readsCtoT_1),
                          "-2", shQuote(readsCtoT_2),
                          "--ff", alignmentParameter, alignmentParameterAdded,
                          "-S", "-p", threads, "--norc", shQuote(readsCtoT_genomeCtoT))
        argsGtoA <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                          "-1", shQuote(readsCtoT_1),
                          "-2",shQuote(readsCtoT_2),
                          "--ff", alignmentParameter, alignmentParameterAdded,
                          "-S", "-p", threads, "--nofw", shQuote(readsCtoT_genomeGtoA))
        print(argsCtoT)
        print(argsGtoA)

        # perform two alignments, readsCtoT against genomeCtoT (plus strand) and
        # readsCtoT against genomeGtoA (minus strand)
        ret1 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        argsCtoT, stdout = TRUE, stderr = TRUE)
        ret2 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        argsGtoA, stdout = TRUE, stderr = TRUE)

        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (GtoA) failed to perform the alignments")
        }

        mrQuSize <- .Call(mergeReorderSam,
                          c(readsCtoT_genomeCtoT, readsCtoT_genomeGtoA),
                          outFile, as.integer(idMode), as.integer(maxHits))
        print(paste("mergeReorderMaxQueueSize", mrQuSize))
    }
}

# For undirected bisulfite, these are the four alignments that are being produced
#
# single end:
#    read        genome
# 0  C->T        C->T Plus
# 1  C->T        G->A Minus
# 2  rc, C->T    C->T Plus
# 3  rc, C->T    G->A Minus
#
# paired end:
#    read1       read2        genome
# 0  C->T        rc, C->T     C->T Plus
# 1  C->T        rc, C->T     G->A Minus
# 2  rc, C->T    C->T *       C->T Plus   -| for 2&3, swap first and second read
# 3  rc, C->T    C->T *       G->A Minus  -|
#
# * reverse complemented twice which cancels out

#' @keywords internal
#' @import Rbowtie
align_RbowtieCtoT_undir <- function(indexDir, reads, samplesFormat, paired,
                                    alignmentParameter, allelic, maxHits,
                                    threads, outFile, cacheDir) {
    # add some variable parameters based on the input format
    if (samplesFormat == "fasta") {
        alignmentParameterAdded <- "-f"
    } else {
        alignmentParameterAdded <- paste("--phred", reads$phred, "-quals", sep = "")
    }

    # set the format of the read identifier. this is necessary for mergeReorderSam outside of this function in allelic mode
    if (!allelic) {
        idMode <- 1
    } else {
        idMode <- 3
    }

    print(paste("Executing bowtie (CtoT and GtoA) on", Sys.info()['nodename'],
                "using", threads, "cores. Parameters:"))
    if (paired == "no") {

        # CtoT convert the reads. include the original sequence in the identifier.
        readsCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName), "readsCtoT_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        readsRcCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName), "readsRcCtoT_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        .Call(convertReadsIdBisRc, reads$FileName, readsCtoT, c("C", "T"), FALSE)
        .Call(convertReadsIdBisRc, reads$FileName, readsRcCtoT, c("C", "T"), TRUE)
        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT), add = TRUE)
        on.exit(file.remove(readsRcCtoT), add = TRUE)

        # create temp filenames for the aligment results
        readsCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName),
                            "readsCtoT_genomeCtoT_", sep = "_"),
            fileext = ".sam"
        )
        readsCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName),
                            "readsCtoT_genomeGtoA_", sep = "_"),
            fileext = ".sam"
        )
        readsRcCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName),
                            "readsRcCtoT_genomeCtoT_", sep = "_"),
            fileext = ".sam"
        )
        readsRcCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName),
                            "readsRcCtoT_genomeGtoA_", sep = "_"),
            fileext = ".sam"
        )

        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsCtoT_genomeGtoA), add = TRUE)
        on.exit(file.remove(readsRcCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsRcCtoT_genomeGtoA), add = TRUE)

        # compile bowtie parameters for the four alignments
        args0 <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                       shQuote(readsCtoT), alignmentParameter,
                       alignmentParameterAdded, "-S", "-p", threads,
                       "--norc", shQuote(readsCtoT_genomeCtoT))
        args1 <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                       shQuote(readsCtoT), alignmentParameter,
                       alignmentParameterAdded, "-S", "-p", threads,
                       "--nofw", shQuote(readsCtoT_genomeGtoA))
        args2 <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                       shQuote(readsRcCtoT), alignmentParameter,
                       alignmentParameterAdded, "-S", "-p", threads,
                       "--norc", shQuote(readsRcCtoT_genomeCtoT))
        args3 <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                       shQuote(readsRcCtoT), alignmentParameter,
                       alignmentParameterAdded, "-S", "-p", threads,
                       "--nofw", shQuote(readsRcCtoT_genomeGtoA))

        print(args0)
        print(args1)
        print(args2)
        print(args3)

        # perform four alignments, readsCtoT agains genomeCtoT (plus strand) and
        # readsCtoT agains genomeGtoA (minus strand)
        ret1 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args0, stdout = TRUE, stderr = TRUE)
        ret2 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args1, stdout = TRUE, stderr = TRUE)
        ret3 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args2, stdout = TRUE, stderr = TRUE)
        ret4 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args3, stdout = TRUE, stderr = TRUE)

        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (GtoA) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (RC, CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (RC, GtoA) failed to perform the alignments")
        }

        mrQuSize <- .Call(mergeReorderSam,
                          c(readsCtoT_genomeCtoT, readsCtoT_genomeGtoA,
                            readsRcCtoT_genomeCtoT, readsRcCtoT_genomeGtoA),
                          outFile, as.integer(idMode), as.integer(maxHits))
        print(paste("mergeReorderMaxQueueSize", mrQuSize))

    } else if (paired == "fr") {

        # CtoT convert the reads. include the original sequence in the identifier.
        readsCtoT_1 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsCtoT_1_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        readsCtoT_2 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsCtoT_2_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        readsRcCtoT_1 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsRcCtoT_1_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )
        readsRcCtoT_2 <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1), "readsRcCtoT_2_", sep = "_"),
            fileext = paste(".", samplesFormat, sep = "")
        )

        .Call(convertReadsIdBisRc, reads$FileName1, readsCtoT_1, c("C", "T"), FALSE)
        .Call(convertReadsIdBisRc, reads$FileName2, readsCtoT_2, c("C", "T"), TRUE) # reverse complement the second read because its fr
        .Call(convertReadsIdBisRc, reads$FileName1, readsRcCtoT_1, c("C", "T"), TRUE)
        .Call(convertReadsIdBisRc, reads$FileName2, readsRcCtoT_2, c("C", "T"), FALSE) # don't reverse complement the second read because its like reverse complementing twice

        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_1), add = TRUE)
        on.exit(file.remove(readsCtoT_2), add = TRUE)
        on.exit(file.remove(readsRcCtoT_1), add = TRUE)
        on.exit(file.remove(readsRcCtoT_2), add = TRUE)

        # create temp filenames for the aligment results
        readsCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsCtoT_genomeCtoT_", sep = "."),
            fileext = ".sam"
        )
        readsCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsCtoT_genomeGtoA_", sep = "."),
            fileext = ".sam"
        )
        readsRcCtoT_genomeCtoT <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsRcCtoT_genomeCtoT_", sep = "."),
            fileext = ".sam"
        )
        readsRcCtoT_genomeGtoA <- tempfile(
            tmpdir = cacheDir,
            pattern = paste(basename(reads$FileName1),
                            "readsRcCtoT_genomeGtoA_", sep = "."),
            fileext = ".sam"
        )

        # make sure that the temp files are deleted
        on.exit(file.remove(readsCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsCtoT_genomeGtoA), add = TRUE)
        on.exit(file.remove(readsRcCtoT_genomeCtoT), add = TRUE)
        on.exit(file.remove(readsRcCtoT_genomeGtoA), add = TRUE)

        args0 <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                       "-1", shQuote(readsCtoT_1),
                       "-2", shQuote(readsCtoT_2),
                       "--ff", alignmentParameter, alignmentParameterAdded,
                       "-S", "-p", threads, "--norc", shQuote(readsCtoT_genomeCtoT))
        args1 <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                       "-1", shQuote(readsCtoT_1),
                       "-2", shQuote(readsCtoT_2),
                       "--ff", alignmentParameter, alignmentParameterAdded,
                       "-S", "-p", threads, "--nofw", shQuote(readsCtoT_genomeGtoA))
        args2 <- paste(shQuote(file.path(indexDir, "bowtieIndexCtoT")),
                       "-2", shQuote(readsRcCtoT_1),
                       "-1", shQuote(readsRcCtoT_2),
                       "--ff", alignmentParameter, alignmentParameterAdded,
                       "-S", "-p", threads, "--norc", shQuote(readsRcCtoT_genomeCtoT))
        args3 <- paste(shQuote(file.path(indexDir, "bowtieIndexGtoA")),
                       "-2", shQuote(readsRcCtoT_1),
                       "-1", shQuote(readsRcCtoT_2),
                       "--ff", alignmentParameter, alignmentParameterAdded,
                       "-S", "-p", threads, "--nofw", shQuote(readsRcCtoT_genomeGtoA))

        print(args0)
        print(args1)
        print(args2)
        print(args3)

        # perform four alignments, readsCtoT agains genomeCtoT (plus strand) and
        # readsCtoT agains genomeGtoA (minus strand)
        ret1 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args0, stdout = TRUE, stderr = TRUE)
        ret2 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args1, stdout = TRUE, stderr = TRUE)
        ret3 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args2, stdout = TRUE, stderr = TRUE)
        ret4 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie"),
                        args3, stdout = TRUE, stderr = TRUE)

        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (GtoA) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret1[length(ret1)]))) {
            stop("bowtie (RC, CtoT) failed to perform the alignments")
        }
        if (!(grepl(" alignments", ret2[length(ret2)]))) {
            stop("bowtie (RC, GtoA) failed to perform the alignments")
        }

        mrQuSize <- .Call(mergeReorderSam,
                          c(readsCtoT_genomeCtoT, readsCtoT_genomeGtoA,
                            readsRcCtoT_genomeCtoT, readsRcCtoT_genomeGtoA),
                          outFile, as.integer(idMode), as.integer(maxHits))
        print(paste("mergeReorderMaxQueueSize", mrQuSize))
    }
}


worker_message <- function(..., sep="", appendLF=TRUE) {
  cat(..., sep=sep)
  if (appendLF) cat("\n")
}


# For a given sample, add an integer to the id. This is necessary for allelic analysis
# when calling mergeReorderSam
#' @keywords internal
#' @importFrom tools file_ext
addNumericToID <- function(reads, paired, cacheDir) {

    if (paired == "no") {
        readsFileNameTmp <- tempfile(
            tmpdir = cacheDir,
            pattern = basename(reads$FileName),
            paste(".", fileext = tools::file_ext(reads$FileName), sep = "")
        )
        .Call(convertReadsIdBisRc, reads$FileName, readsFileNameTmp, NULL, FALSE)
        reads$FileName <- readsFileNameTmp
    } else {
        readsFileNameTmp1 <- tempfile(
            tmpdir = cacheDir,
            pattern = basename(reads$FileName1),
            paste(".", fileext = tools::file_ext(reads$FileName1), sep = "")
        )
        readsFileNameTmp2 <- tempfile(
            tmpdir = cacheDir,
            pattern = basename(reads$FileName2),
            paste(".", fileext = tools::file_ext(reads$FileName2), sep = "")
        )
        .Call(convertReadsIdBisRc, reads$FileName1, readsFileNameTmp1, NULL, FALSE)
        .Call(convertReadsIdBisRc, reads$FileName2, readsFileNameTmp2, NULL, FALSE)
        reads$FileName1 <- readsFileNameTmp1
        reads$FileName2 <- readsFileNameTmp2
    }
    return(reads)
}


# helper function executed in parallel by samToSortedBamParallel
#' @keywords internal
#' @importFrom Rsamtools asBam sortBam
samToSortedBamCore <- function(ind, paramsL) {
    set.seed(0)
    params <- paramsL[[ind]]

    # don't sort in the case of splitChrSam_unaligned
    if (basename(params[2]) != "splitChrSam_unaligned") {
        bamtempFile <- tempfile()
        Rsamtools::asBam(params[1], bamtempFile, indexDestination = FALSE)
        on.exit(file.remove(bamtempFile))
        Rsamtools::sortBam(paste(bamtempFile, "bam", sep = "."), params[2])
    } else {
        Rsamtools::asBam(params[1], params[2], indexDestination = FALSE)
    }
}

# convert a sam file to a sorted bam file in parallel fashion using p threads.
# it splits the sam file into individual chromosomes, creates a cluster object with
# p nodes and sort and converts to bam in parallel fashion. After this step, the bam
# files are concatentated (in the order of the header of the original sam file and
# an index in built for the final bam file
#' @keywords internal
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterApplyLB
#' @importFrom Rsamtools indexBam
samToSortedBamParallel <- function(file, destination, p, cacheDir = NULL) {
    # test if the input file exists
    if (!file.exists(file)) {
        stop("Cannot read ", file, call. = FALSE)
    }

    if (is.null(cacheDir))
        cacheDir <- tempdir()

    # split the input sam file into seperate files, one per chromosome in
    # a temporary directory
    splitDir <- tempfile(tmpdir = cacheDir, pattern = "samToBam_")
    if (!dir.create(path = splitDir, showWarnings = FALSE)) {
        stop("No permissions to create a directory in the cacheDir", call. = FALSE)
    }
    # make sure that the temp dir is deleted
    on.exit(unlink(splitDir, recursive = TRUE))

    # perform the split
    chrNames <- .Call(splitSamChr, file, splitDir)

    # collect the sam files that are not empty. The empty ones would cause a
    # problem during sam to bam conversion
    splitSamAndBamNames <- NULL
    for (i in seq_len(length(chrNames))) {
        samName <- file.path(splitDir, paste(chrNames[i], "sam", sep = "."))
        bamName <- file.path(splitDir, chrNames[i])
        con <- file(samName, "r")
        fileHead <- scan(con, as.list(rep('character', 20)),
                         sep = "\t", comment.char = "@", fill = TRUE,
                         nmax = 30, quiet = TRUE)
        if (length(fileHead[[1]]) > 0) {
            splitSamAndBamNames[[length(splitSamAndBamNames) + 1]] <-
                c(samName, bamName)
        }
        close(con)
    }

    # sort the jobs based on the file sizes
    sortedOrder <- order(file.info(do.call(rbind, splitSamAndBamNames)[, 1])$size,
                         decreasing = TRUE)

    clObjS <- parallel::makeCluster(p)
    on.exit(parallel::stopCluster(clObjS), add = TRUE)

    # load QuasR package on all the nodes
    loadQuasR(clObjS)

    parallel::clusterApplyLB(clObjS, seq_len(length(splitSamAndBamNames)),
                             samToSortedBamCore,
                             splitSamAndBamNames[sortedOrder])

    # concatenate all the sorted bam files in the order of the original sam file header
    .Call(catBam, paste0(do.call(rbind, splitSamAndBamNames)[, 2], ".bam"),
          paste0(destination, ".bam"))

    # create index for the final bam file
    Rsamtools::indexBam(paste0(destination, ".bam"))

    invisible(paste0(destination, ".bam"))
}

