# build an index based on a given BSgenome, create a package with the files and install it on the system
#' @keywords internal
#' @importFrom Biobase createPackage
#' @importFrom utils installed.packages install.packages
buildIndexPackage <- function(genome, aligner, alnModeID, cacheDir, lib.loc) {
    indexPackageName <- paste(genome, alnModeID, sep = ".")
    
    lib.locTemp <- lib.loc
    if (is.na(lib.loc)){
        lib.locTemp <- NULL
    } # this is a way to convert NA to NULL needed for install.packages
    
    # Create the index and install it if it is not yet installed on the system
    if (!(indexPackageName %in% utils::installed.packages(lib.loc = lib.locTemp)[, 'Package'])) {
        genomeObj <- get(genome) # access the BSgenome
        # flush the BSgenome to disk
        fastaFilepath <- BSgenomeSeqToFasta(genomeObj, tempfile(tmpdir = cacheDir, 
                                                                fileext = ".fa"))  
        on.exit(unlink(fastaFilepath))

        # create the package (without installing it yet). if it is already in 
        # the temp dir, keep it. It might contain an index
        # calculated in a previous round where the installation was not successful 
        # (due to permissions)
        if (!file.exists(file.path(cacheDir, indexPackageName))) {
            seedList <- createSeedList(genomeObj, aligner, indexPackageName)
            templatePath <- system.file("AlignerIndexPkg-template", package = "QuasR")
            Biobase::createPackage(indexPackageName, cacheDir, templatePath, 
                                   seedList, quiet = TRUE)
        }
        # create the index files
        buildIndex(fastaFilepath,file.path(cacheDir, indexPackageName, "inst", 
                                           "alignmentIndex"), alnModeID, cacheDir)
        
        # install the package
        utils::install.packages(
            file.path(cacheDir, indexPackageName), repos = NULL, 
            dependencies = FALSE, type = "source", lib = lib.locTemp
        )
        
        if (indexPackageName %in% utils::installed.packages(lib.loc = lib.locTemp)[, 'Package']) {
            # package installation was successful, clean up
            unlink(file.path(cacheDir, indexPackageName), recursive = TRUE)
        } else {
            stop("Failed to install the index package '", indexPackageName, "'")
        }
    }
}

#' @keywords internal
#' @importFrom utils read.table
#' @importFrom Rsamtools scanFaIndex scanFa indexFa
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings masks injectHardMask replaceLetterAt DNAStringSet
#'   writeXStringSet
buildIndexSNP <- function(snpFile, indexPath, genome, genomeFormat, 
                          alnModeID, cacheDir, checkMD5 = FALSE) {

    fastaOutFileR <- paste(snpFile, basename(genome), "R", "fa", sep = ".")
    fastaOutFileA <- paste(snpFile, basename(genome), "A", "fa", sep = ".")
    
    if (!file.exists(fastaOutFileR) | !file.exists(fastaOutFileA)) {
        # read in the SNPs
        message(paste("Reading and processing the SNP file:", snpFile))
        snps <- utils::read.table(snpFile, colClasses = c("factor", "numeric", 
                                                          "character", "character"))
        colnames(snps) <- c("chrom", "pos", "ref", "alt")
 
        # convert ref nucleotides to upper case if not already the case
        snps[, 3] <- toupper(snps[, 3]) 
        # convert alt nucleotides to upper case if not already the case
        snps[, 4] <- toupper(snps[, 4]) 
        
        # check if there are only regular nucleotides (ACGT) in the snp file
        if (!all(unique(c(snps[, 3], snps[, 4])) %in% c("A", "C", "G", "T"))) {
            stop("There are non-regular nucleoides in snpFile. Only ACGT are allowed.",
                 call. = FALSE)
        }
        
        snpsL <- split(snps, snps[, 1]) # split SNPs accoring to chromosome (first column)
        
        # check for duplicate snp entries (not allowed)
        if (sum(sapply(snpsL,function(x){
            sum(duplicated(x[, 2]))
        })) > 0) {
            stop("There are duplicate SNP positions in snpFile.", call. = FALSE)
        }
        
        if (genomeFormat == "file") {
            idx <- Rsamtools::scanFaIndex(genome)
            allChrs <- as.character(GenomeInfoDb::seqnames(idx))
        } else {
            genomeObj <- get(genome) # access the BSgenome
            allChrs <- GenomeInfoDb::seqnames(genomeObj)
        }
        if (!all(names(snpsL) %in% allChrs)){
            stop("The snpFile contains chromosomes that are not present in the genome",
                 call. = FALSE)
        }
    }

    for (j in seq_len(2)) {
        fastaOutFile <- c(fastaOutFileR, fastaOutFileA)[j]
        
        if (!file.exists(fastaOutFile)) {
            append <- FALSE
            message(paste("Creating the genome fasta file containing the SNPs:",
                          fastaOutFile))
            for (i in seq_len(length(allChrs))) {
                if (genomeFormat == "file") {
                    seq <- Rsamtools::scanFa(genome, idx[i])[[1]]
                } else {
                    if (is.null(Biostrings::masks(genomeObj[[allChrs[i]]]))) {
                        seq <- genomeObj[[allChrs[i]]]
                    } else {
                        seq <- Biostrings::injectHardMask(genomeObj[[allChrs[i]]],
                                                          letter = "N")
                    }
                }
                snpsLC <- snpsL[[allChrs[i]]]
                # check if the ref nucleotide in snpFile matches the genome
                if (!is.null(snpsLC)) {
                    # inject the SNPs
                    seqSNP <- Biostrings::replaceLetterAt(seq, snpsLC[, 2], 
                                                          snpsLC[, 2 + j])
                    seqSNP_SS <- Biostrings::DNAStringSet(seqSNP)
                } else {
                    seqSNP_SS <- Biostrings::DNAStringSet(seq)
                }
                names(seqSNP_SS) <- allChrs[i]
                
                Biostrings::writeXStringSet(seqSNP_SS, filepath = fastaOutFile, 
                                            format = "fasta", append = append)
                append <- TRUE
            }
        }
    }
    
    # create .fai file for the snp genome
    for (fastaOutFile in c(fastaOutFileR, fastaOutFileA)) {
        if (!file.exists(paste(fastaOutFile, "fai", sep = "."))) {
            message(paste("Creating a .fai file for the snp genome:", fastaOutFile))
            if (class(try(Rsamtools::indexFa(fastaOutFile))) == "try-error") {
                stop("Cannot write into the directory where ", fastaOutFile, 
                     " is located. Make sure you have the right permissions",
                     call. = FALSE)
            }
        }
    }

    # create The two indices
    buildIndex(fastaOutFileR, paste(snpFile, basename(genome), "R", "fa",
                                    alnModeID, sep = "."), alnModeID, cacheDir)
    buildIndex(fastaOutFileA, paste(snpFile, basename(genome), "A", "fa",
                                    alnModeID, sep = "."), alnModeID, cacheDir)
}



# build an index based on a given fasta file with one or more sequences.
#' @keywords internal
#' @importFrom utils read.delim write.table
#' @importFrom tools md5sum
buildIndex <- function(seqFile, indexPath, alnModeID, cacheDir, checkMD5 = FALSE) {
    
    # check if the directory exists but contains no ref_md5Sum.txt file. This means that a last index builder call did not
    # finish properly. Delete the directory containing the partial index
    if (file.exists(indexPath)) {
        if (!("ref_md5Sum.txt" %in% dir(indexPath))) {
            message(paste("Deleting an incomplete index at:", indexPath))
            file.remove(dir(indexPath, full.names = TRUE))
            unlink(indexPath, recursive = TRUE)
        } else {
            # if checkMD5==TRUE, check consistency between the sequences and the index
            # delete index in the case of inconsistency
            if (checkMD5) {
                MD5_fromIndexTab <- utils::read.delim(
                    file.path(indexPath, "ref_md5Sum.txt"), 
                    header = FALSE, colClasses = "character"
                )
                if (!all(dim(MD5_fromIndexTab) == c(1, 1))){
                    stop("The md5sum file does not have the right format: ",
                         file.path(indexPath, "ref_md5Sum.txt"), call. = FALSE)
                }
                MD5_fromIndex <- MD5_fromIndexTab[1, 1]
                MD5_fromSeq <- tools::md5sum(seqFile)
                if (MD5_fromIndex != MD5_fromSeq) {
                    message(paste("The sequence file", seqFile,
                                  "was changed. Updating the index: ", indexPath)) 
                    file.remove(dir(indexPath, full.names = TRUE))
                    unlink(indexPath, recursive = TRUE)
                }
            }
        }
    }
    if (!file.exists(indexPath)) {
        # create the directory for the index
        if (!dir.create(indexPath)){
            stop("Cannot create the directory: ", indexPath, call. = FALSE)
        }
        
        # Create all the various indices
        message(paste("Creating an", alnModeID, "index for", seqFile))

        if (alnModeID == "Rbowtie") {
            ret <- buildIndex_Rbowtie(seqFile, indexPath)
        } else if (alnModeID == "RbowtieCtoT") {
            ret <- buildIndex_RbowtieCtoT(seqFile, indexPath, cacheDir)
        } else if (alnModeID == "RbowtieCs") {
            ret <- buildIndex_RbowtieCs(seqFile, indexPath)
        } else if (alnModeID == "Rhisat2") {
            ret <- buildIndex_Rhisat2(seqFile, indexPath)
        } else {
            stop("Fatal error 2374027")
        }
        
        if (ret == 0) {
            indexMD5 <- tools::md5sum(seqFile)
            utils::write.table(indexMD5, file = file.path(indexPath, "ref_md5Sum.txt"),
                               row.names = FALSE, col.names = FALSE, 
                               sep = "\t", quote = FALSE)
            message("Finished creating index")
        } else {
            # the execution of index-build failed, delete the directory.
            file.remove(dir(indexPath, full.names = TRUE))
            unlink(indexPath, recursive = TRUE)
            stop("The execution of the index builder failed. No index was created", 
                 call. = FALSE)
        }
    }
}

# build index for Rhisat2, base space, non-bisulfite converted
#' @keywords internal
buildIndex_Rhisat2 <- function(seqFile,indexPath) {
    indexFullPath <- file.path(indexPath, "hisat2Index")
    
    ret <- system2(file.path(system.file(package = "Rhisat2"), "hisat2-build"),
                   c(shQuote(seqFile), shQuote(indexFullPath)), 
                   stdout = TRUE, stderr = TRUE)
    if (!(grepl("^Total time for call to driver", ret[length(ret)]))) {
        ret <- 1
    } else {
        ret <- 0
    }
    return(ret)
}

# generate splice site file
#' @keywords internal
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom AnnotationDbi loadDb
#' @importFrom tools md5sum
#' @importFrom utils read.delim write.table
buildSpliceSiteFile <- function(geneAnnotation, geneAnnotationFormat) {
    if (file.exists(paste0(geneAnnotation, ".SpliceSites.txt"))) {
        # if the SpliceSites.txt file exists, but not the md5sum file, remove the former
        if (!file.exists(paste0(geneAnnotation, ".SpliceSites.txt.md5"))) {
            message("Deleting an incomplete SpliceSites.txt file at: ",
                    paste0(geneAnnotation, ".SpliceSites.txt"))
            unlink(paste0(geneAnnotation, ".SpliceSites.txt"))
        } else {
            # both SpliceSites.txt and md5 file exist
            md5FromFile <- utils::read.delim(
                paste0(geneAnnotation, ".SpliceSites.txt.md5"),
                header = FALSE, colClasses = "character"
            )
            if (!all(dim(md5FromFile) == c(1, 1))) {
                stop("The md5sum file does not have the right format: ",
                     paste0(geneAnnotation, ".SpliceSites.txt.md5"), call. = FALSE)
            }
            md5FromFile <- md5FromFile[1, 1]
            md5FromObj <- tools::md5sum(geneAnnotation)
            if (md5FromFile != md5FromObj) {
                message("The annotation file ", geneAnnotation,
                        " was changed. Recreating the SpliceSites file.")
                unlink(paste0(geneAnnotation, ".SpliceSites.txt"))
            }
        }
    }
    if (!file.exists(paste0(geneAnnotation, ".SpliceSites.txt"))) {
        if (geneAnnotationFormat == "TxDb") {
            txdb <- AnnotationDbi::loadDb(geneAnnotation)
        } else if (geneAnnotationFormat == "gtf") {
            txdb <- suppressWarnings(
                GenomicFeatures::makeTxDbFromGFF(geneAnnotation, format = "auto")
            )
        } else {
            stop("Fatal error 81956293")
        }
        Rhisat2::extract_splice_sites(txdb, paste0(geneAnnotation, ".SpliceSites.txt"), 
                                      min_length = 5)
        md5FromObj <- tools::md5sum(geneAnnotation)
        utils::write.table(
            md5FromObj, file = paste0(geneAnnotation, ".SpliceSites.txt.md5"),
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
        )
    }
}

# build index for Rbowtie, base space, non-bisulfite converted
#' @keywords internal
#' @import Rbowtie
buildIndex_Rbowtie <- function(seqFile, indexPath) {
    indexFullPath <- file.path(indexPath, "bowtieIndex")
    
    ret <- system2(file.path(system.file(package = "Rbowtie"), "bowtie-build"),
                   c(shQuote(seqFile), shQuote(indexFullPath)), stdout = TRUE, stderr = TRUE)
    if (!(grepl("^Total time for backward call to driver", ret[length(ret)]))) {
        ret <- 1
    } else {
        ret <- 0
    }
    
    return(ret)
}

# build index for Rbowtie, base space, bisulfite converted
#' @keywords internal
#' @importFrom Rsamtools scanFaIndex scanFa
#' @importFrom Biostrings writeXStringSet chartr
buildIndex_RbowtieCtoT <- function(seqFile, indexPath, cacheDir) {
    indexFullPath <- file.path(indexPath, "bowtieIndex")
    
    # read the reference sequences
    idx <- Rsamtools::scanFaIndex(seqFile)
    
    # create two temporary sequences files, C->T the plus strand and G->A for the minus strand
    outFilePlus <- tempfile(tmpdir = cacheDir, fileext = ".fa")
    outFileMinus <- tempfile(tmpdir = cacheDir, fileext = ".fa")
    
    on.exit(unlink(c(outFilePlus, outFileMinus)))
    
    # read one chromosome after the other, convert and write to disk
    append <- FALSE
    for (i in seq_len(length(idx))) {
        seq <- Rsamtools::scanFa(seqFile, idx[i])
        plus_strand <- Biostrings::chartr("C", "T", seq)
        minus_strand <- Biostrings::chartr("G", "A", seq)
        Biostrings::writeXStringSet(plus_strand, filepath = outFilePlus, 
                                    format = "fasta", append = append)
        Biostrings::writeXStringSet(minus_strand, filepath = outFileMinus, 
                                    format = "fasta", append = append)
        append <- TRUE
    }

    # execute bowtie twice to create the two indices
    ret1 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie-build"),
                    c(shQuote(outFilePlus), 
                      shQuote(file.path(indexPath, "bowtieIndexCtoT"))), 
                    stdout = TRUE, stderr = TRUE)
    if (!(grepl("^Total time for backward call to driver", ret1[length(ret1)]))) {
        ret1 <- 1
    } else {
        ret1 <- 0
    }
    
    ret2 <- system2(file.path(system.file(package = "Rbowtie"), "bowtie-build"),
                    c(shQuote(outFileMinus), 
                      shQuote(file.path(indexPath, "bowtieIndexGtoA"))), 
                    stdout = TRUE, stderr = TRUE)
    if (!(grepl("^Total time for backward call to driver", ret2[length(ret2)]))) {
        ret2 <- 1
    } else {
        ret2 <- 0
    }
 
    if ((ret1 == 0) & (ret2 == 0)) {
        return(0)
    } else {
        return(1)
    }
}

# build index for Rbowtie, color space, non-bisulfite converted
#' @keywords internal
#' @import Rbowtie
buildIndex_RbowtieCs <- function(seqFile, indexPath) {
    indexFullPath <- file.path(indexPath, "bowtieIndexCs")

    ret <- system2(file.path(system.file(package = "Rbowtie"), "bowtie-build"),
                   c(shQuote(seqFile), shQuote(indexFullPath), "-C"), 
                   stdout = TRUE, stderr = TRUE)
    if (!(grepl("^Total time for backward call to driver", ret[length(ret)]))) {
        ret <- 1
    } else {
        ret <- 0
    }
    
    return(ret)
}


# flush all chromosomes of a BSgenome into a file
#' @keywords internal
#' @importFrom methods is
#' @importFrom Biostrings masks DNAStringSet writeXStringSet injectHardMask
#' @importFrom GenomeInfoDb seqnames
BSgenomeSeqToFasta <- function(bsgenome, outFile) {
    if (!methods::is(bsgenome, "BSgenome")) {
        stop("The variable 'bsgenome' is not a BSgenome")
    }
    append <- FALSE
    for (chrT in GenomeInfoDb::seqnames(bsgenome)) {
        if (is.null(Biostrings::masks(bsgenome[[chrT]])))
            chrSeq <- Biostrings::DNAStringSet(bsgenome[[chrT]])
        else
            chrSeq <- Biostrings::DNAStringSet(
                Biostrings::injectHardMask(bsgenome[[chrT]], 
                                           letter = "N"))
        names(chrSeq) <- chrT
        Biostrings::writeXStringSet(chrSeq, filepath = outFile, 
                                    format = "fasta", append = append)
        append <- TRUE
    }
    return(outFile)
}

#' @keywords internal
#' @importFrom S4Vectors metadata
#' @importFrom utils packageVersion
#' @importFrom GenomeInfoDb bsgenomeName provider releaseDate organism 
#'   commonName
createSeedList <- function(genome, aligner, indexPackageName) {
    pv <- S4Vectors::metadata(genome)$genome
    seed <- list(##package seeds
        PKGTITLE = paste(aligner, "Index of", GenomeInfoDb::bsgenomeName(genome)),
        PKGDESCRIPTION = paste(aligner, "Index of", GenomeInfoDb::bsgenomeName(genome)),
        PKGVERSION = "1.0",
        DATE = format(Sys.time(), "%Y-%M-%d"),
        AUTHOR = "QuasR",
        MAINTAINER = "This package was automatically created <yourfault@somewhere.net>",
        LIC = paste("see", GenomeInfoDb::bsgenomeName(genome)),
        PKGDETAILS = "Storing genome index for BSgenome which is needed for alignments.",
        PKGEXAMPLES = "No examples",
        
        ##genome seeds
        GENOMENAME = GenomeInfoDb::bsgenomeName(genome),
        PROVIDER = GenomeInfoDb::provider(genome),
        PROVIDERVERSION = ifelse(!is.null(pv) && is.character(pv) && 
                                     length(pv) == 1L, pv, "not_available"),
        RELEASEDATE = GenomeInfoDb::releaseDate(genome),
        RELEASENAME = "not_available",
        ORGANISM = GenomeInfoDb::organism(genome),
        SPECIES = GenomeInfoDb::commonName(genome),
        SRCDATAFILES = GenomeInfoDb::bsgenomeName(genome),
        ORGANISMBIOCVIEW = gsub(" ", "_", GenomeInfoDb::organism(genome)),
        
        #aligner seeds
        ALIGNER = aligner,
        ALIGNERVERSION = as.character(utils::packageVersion(aligner))
    )
    
    return(seed)
}

