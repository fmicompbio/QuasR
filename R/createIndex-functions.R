# build an index based on a given BSgenome, create a package with the files and install it on the system
buildIndexPackage <- function(genome,aligner,alnModeID,cacheDir,lib.loc){
  indexPackageName <- paste(genome,alnModeID,sep=".")

  # Create the index and install it if it is not yet installed on the system
  if(!(indexPackageName %in% installed.packages())){
    genomeObj <- get(genome) # access the BSgenome
    fastaFilepath <- BSgenomeSeqToFasta(genomeObj, tempfile(tmpdir=cacheDir, fileext=".fa"))  # flush the BSgenome to disk
    on.exit(unlink(fastaFilepath))

    # create the package (without installing it yet). if it is already in the temp dir, keep it. It might contain an index
    # calculated in a previous round where the intallation was not successful (due to permissions)
    if(!file.exists(file.path(cacheDir,indexPackageName))){
      seedList <- createSeedList(genomeObj, aligner, indexPackageName)
      templatePath <- system.file("AlignerIndexPkg-template", package="QuasR")
      Biobase::createPackage(indexPackageName, cacheDir, templatePath, seedList, quiet=TRUE)
    }
    # create the index files
    buildIndex(fastaFilepath,file.path(cacheDir,indexPackageName,"inst","alignmentIndex"),alnModeID,cacheDir)

    # install the package
    lib.locTemp <- lib.loc
    if(is.na(lib.loc)){lib.locTemp<-NULL;} # this is a way to convert NA to NULL needed for install.packages
    install.packages(file.path(cacheDir,indexPackageName), repos=NULL, dependencies=FALSE, type="source", lib=lib.locTemp)

    if(indexPackageName %in% installed.packages()[,'Package']){
      # package installation was successful, clean up
      unlink(file.path(cacheDir,indexPackageName),recursive=TRUE)
    }else{"Fatal error 2309420"}
  }
}


buildIndexSNP <- function(snpFile,indexPath,genome,genomeFormat,alnModeID,cacheDir,checkMD5=FALSE){

  fastaOutFileR <- paste(snpFile,basename(genome),"R","fa",sep=".")
  fastaOutFileA <- paste(snpFile,basename(genome),"A","fa",sep=".")
 
  if(!file.exists(fastaOutFileR) | !file.exists(fastaOutFileA)){
    # read in the SNPs
    message(paste("Reading and processing the SNP file:",snpFile))
    snps <- read.table(snpFile,colClasses=c("factor","numeric","character","character"))
    colnames(snps) <- c("chrom","pos","ref","alt")
 
    snps[,3] <- toupper(snps[,3]) # convert ref nucleotides to upper case if not already the case
    snps[,4] <- toupper(snps[,4]) # convert alt nucleotides to upper case if not already the case

    # check if there are only regular nucleotides (ACGT) in the snp file
    if(!all(unique(c(snps[,3],snps[,4])) %in% c("A","C","G","T"))){stop("There are non-regular nucleoides in snpFile. Only ACGT are allowed.",call.=FALSE)}

    snpsL <- split(snps,snps[,1]) # split SNPs accoring to chromosome (first column)

    # check for duplicate snp entries (not allowed)
    if(sum(sapply(snpsL,function(x){sum(duplicated(x[,2]))})) > 0){stop("There are duplicate SNP positions in snpFile.",call.=FALSE)}

    if(genomeFormat=="file"){
      idx <- scanFaIndex(genome)
      allChrs <- as.character(seqnames(idx))
    }else{
      genomeObj <- get(genome) # access the BSgenome
      allChrs <- seqnames(genomeObj)
    }
    if(!all(names(snpsL) %in% allChrs)){stop("The snpFile contains chromosomes that are not present in the genome",call.=FALSE)}
  }

  for(j in 1:2){
    fastaOutFile <- c(fastaOutFileR,fastaOutFileA)[j]

    if(!file.exists(fastaOutFile)){
      append <- FALSE
      message(paste("Creating the genome fasta file containing the SNPs:",fastaOutFile))
      for(i in 1:length(allChrs)){
        if(genomeFormat=="file"){
          seq <- scanFa(genome, idx[i])[[1]]
        }else{
          if(is.null(masks(genomeObj[[allChrs[i]]]))){
            seq <- genomeObj[[allChrs[i]]]
          }else{
            seq <- injectHardMask(genomeObj[[allChrs[i]]], letter="N")
          }
        }
        snpsLC <- snpsL[[allChrs[i]]]
        # check if the ref nucleotide in snpFile matches the genome
        if(!is.null(snpsLC)){
          # inject the SNPs
          seqSNP <- replaceLetterAt(seq,snpsLC[,2],snpsLC[,2+j])
          seqSNP_SS <- DNAStringSet(seqSNP)
        }else{
          seqSNP_SS <- DNAStringSet(seq)
        }
        names(seqSNP_SS) <- allChrs[i]

        writeXStringSet(seqSNP_SS, filepath=fastaOutFile, format="fasta", append=append)
        append <- TRUE
      }
    }
  }

  # create .fai file for the snp genome
  for(fastaOutFile in c(fastaOutFileR,fastaOutFileA)){
    if(!file.exists(paste(fastaOutFile,"fai",sep="."))){
      message(paste("Creating a .fai file for the snp genome:",fastaOutFile))
      if(class(try(indexFa(fastaOutFile)))=="try-error"){
        stop("Cannot write into the directory where ",fastaOutFile," is located. Make sure you have the right permissions",call.=FALSE)
      }
    }
  }

  # create The two indices
  buildIndex(fastaOutFileR,paste(snpFile,basename(genome),"R","fa",alnModeID,sep="."),alnModeID,cacheDir)
  buildIndex(fastaOutFileA,paste(snpFile,basename(genome),"A","fa",alnModeID,sep="."),alnModeID,cacheDir)

}



# build an index based on a given fasta file with one or more sequences.
buildIndex <- function(seqFile,indexPath,alnModeID,cacheDir,checkMD5=FALSE){

  # check if the directory exists but contains no ref_md5Sum.txt file. This means that a last index builder call did not
  # finish properly. Delete the directory containing the partial index
  if(file.exists(indexPath)){
    if(!("ref_md5Sum.txt" %in% dir(indexPath))){
       message(paste("Deleting an incomplete index at:",indexPath))
       file.remove(dir(indexPath,full.names=TRUE))
       unlink(indexPath,recursive=TRUE)
    }else{
      # if checkMD5==TRUE, check consistency between the sequences and the index
      # delete index in the case of inconsistency
      if(checkMD5){
        MD5_fromIndexTab <- read.delim(file.path(indexPath,"ref_md5Sum.txt"),header=FALSE)
        if(!all(dim(MD5_fromIndexTab)==c(1,1))){stop("The md5sum file does not have the right format: ",file.path(indexPath,"ref_md5Sum.txt"),call.=FALSE)}
        MD5_fromIndex <- MD5_fromIndexTab[1,1]
        MD5_fromSeq <- tools::md5sum(seqFile)
        if(MD5_fromIndex != MD5_fromSeq){
          message(paste("The sequence file",seqFile,"was changed. Updating the index: ",indexPath)) 
          file.remove(dir(indexPath,full.names=TRUE))
          unlink(indexPath,recursive=TRUE)
        }
      }
    }
  }
  if(!file.exists(indexPath)){
    # create the directory for the index
    if(!dir.create(indexPath)){stop("Cannot create the directory: ",indexPath,call.=FALSE)}

    # Create all the various indices
    message(paste("Creating an",alnModeID,"index for",seqFile))

    if(alnModeID=="Rbowtie"){
      ret <- buildIndex_Rbowtie(seqFile,indexPath)
    }else if(alnModeID=="RbowtieCtoT"){
      ret <- buildIndex_RbowtieCtoT(seqFile,indexPath,cacheDir)
    }else if(alnModeID=="RbowtieCs"){
      ret <- buildIndex_RbowtieCs(seqFile,indexPath)
    }else{stop("Fatal error 2374027")}

    if(ret==0){
      indexMD5 <- tools::md5sum(seqFile)
      write.table(indexMD5,file=file.path(indexPath,"ref_md5Sum.txt"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
      message("Finished creating index")
    }else{
      # the execution of index-build failed, delete the directory.
      file.remove(dir(indexPath,full.names=TRUE))
      unlink(indexPath,recursive=TRUE)
      stop("The execution of the index builder failed. No index was created",call.=FALSE)
    }
  }
}


# build index for Rbowtie, base space, non-bisulfite converted
buildIndex_Rbowtie <- function(seqFile,indexPath){
  indexFullPath <- file.path(indexPath,"bowtieIndex")

  ret <- system2(file.path(system.file(package="Rbowtie"),"bowtie-build"),c(shQuote(seqFile),shQuote(indexFullPath)), stdout=TRUE, stderr=TRUE)
  if(!(grepl("^Total time for backward call to driver", ret[length(ret)]))){
    ret <- 1
  }else{
    ret <- 0
  }
 
  return(ret)
}

# build index for Rbowtie, base space, bisulfite converted
buildIndex_RbowtieCtoT <- function(seqFile,indexPath,cacheDir){
  indexFullPath <- file.path(indexPath,"bowtieIndex")

  # read the reference sequences
  idx <- scanFaIndex(seqFile)

  # create two temporary sequences files, C->T the plus strand and G->A for the minus strand
  outFilePlus <- tempfile(tmpdir=cacheDir, fileext=".fa")
  outFileMinus <- tempfile(tmpdir=cacheDir, fileext=".fa")

  on.exit(unlink(c(outFilePlus,outFileMinus)))

  # read one chromosome after the other, convert and write to disk
  append <- FALSE
  for(i in 1:length(idx)){
    seq <- scanFa(seqFile, idx[i])
    plus_strand <- chartr("C", "T", seq)
    minus_strand <- chartr("G", "A", seq)
    writeXStringSet(plus_strand, filepath=outFilePlus, format="fasta", append=append)
    writeXStringSet(minus_strand, filepath=outFileMinus, format="fasta", append=append)
    append <- TRUE
  }

  # execute bowtie twice to create the two indices
  ret1 <- system2(file.path(system.file(package="Rbowtie"),"bowtie-build"),c(shQuote(outFilePlus),shQuote(file.path(indexPath,"bowtieIndexCtoT"))), stdout=TRUE, stderr=TRUE)
  if(!(grepl("^Total time for backward call to driver", ret1[length(ret1)]))){
    ret1 <- 1
  }else{
    ret1 <- 0
  }

  ret2 <- system2(file.path(system.file(package="Rbowtie"),"bowtie-build"),c(shQuote(outFileMinus),shQuote(file.path(indexPath,"bowtieIndexGtoA"))), stdout=TRUE, stderr=TRUE)
  if(!(grepl("^Total time for backward call to driver", ret2[length(ret2)]))){
    ret2 <- 1
  }else{
    ret2 <- 0
  }
 
  if((ret1==0) & (ret2==0)){
    return(0)
  }else{
    return(1)
  }
}

# build index for Rbowtie, color space, non-bisulfite converted
buildIndex_RbowtieCs <- function(seqFile,indexPath){
  indexFullPath <- file.path(indexPath,"bowtieIndexCs")

  ret <- system2(file.path(system.file(package="Rbowtie"),"bowtie-build"),c(shQuote(seqFile),shQuote(indexFullPath),"-C"), stdout=TRUE, stderr=TRUE)
  if(!(grepl("^Total time for backward call to driver", ret[length(ret)]))){
    ret <- 1
  }else{
    ret <- 0
  }

  return(ret)
}


# flush all chromosomes of a BSgenome into a file
BSgenomeSeqToFasta <- function(bsgenome, outFile)
{
    if(!is(bsgenome, "BSgenome")){stop("The variable 'bsgenome' is not a BSgenome")}
    append <- FALSE
    for(chrT in seqnames(bsgenome)){
        if(is.null(masks(bsgenome[[chrT]])))
            chrSeq <- DNAStringSet(bsgenome[[chrT]])
        else
            chrSeq <- DNAStringSet(injectHardMask(bsgenome[[chrT]], letter="N"))
        names(chrSeq) <- chrT
        writeXStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
        append <- TRUE
    }
    return(outFile)
}



createSeedList <- function(genome, aligner, indexPackageName)
{
    seed <- list(##package seeds
                 PKGTITLE=paste(aligner, "Index of", bsgenomeName(genome)),
                 PKGDESCRIPTION=paste(aligner, "Index of", bsgenomeName(genome)),
                 PKGVERSION="1.0",
                 DATE=format(Sys.time(), "%Y-%M-%d"),
                 AUTHOR="QuasR",
                 MAINTAINER="This package was automatically created <yourfault@somewhere.net>",
                 LIC=paste("see",bsgenomeName(genome)),
                 PKGDETAILS="Storing genome index for BSgenome which is needed for alignments.",
                 PKGEXAMPLES="No examples",

                 ##genome seeds
                 GENOMENAME=bsgenomeName(genome),
                 PROVIDER=provider(genome),
                 PROVIDERVERSION=providerVersion(genome),
                 RELEASEDATE=releaseDate(genome),
                 RELEASENAME=releaseName(genome),
                 ORGANISM=organism(genome),
                 SPECIES=species(genome),
                 SRCDATAFILES=bsgenomeName(genome),
                 ORGANISMBIOCVIEW=gsub(" ", "_", organism(genome)),

                 #aligner seeds
                 ALIGNER=aligner,
                 ALIGNERVERSION=installed.packages()[aligner, 'Version']
                 )

    return(seed)
}

