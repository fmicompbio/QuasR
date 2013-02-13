# bowtie alignments: --best --strata -m 1 (return one alignment, but only if it is a unique mapper)
# plus -q/-f [--phred33/64-quals] [-C] -S --quiet


# paired:
# valid values are NULL, "no", "fr", "rf", "ff"
# in case of NULL, QuasR determines the paired status automatically from the given sampleFile (2 columns:unpaired, 3 columns:paired).
# The orientation for paired end is set to "fr" by default. In case of bam input files, the user needs to specify manually the 
# status ("no","fr","rf","ff")
# snpFile need four columns (chrom, pos, ref, alt) with no header

qAlign <- function(sampleFile, genome, auxiliaryFile=NULL, aligner="Rbowtie", maxHits=1, paired=NULL,
                   splicedAlignment=FALSE, snpFile=NULL, bisulfite="no", alignmentParameter=NULL, projectName="qProject",
                   alignmentsDir=NULL, lib.loc=NULL, cacheDir=NULL, clObj=NULL){

  # check if the user provided a samplefile and a genome
  if(missing(sampleFile)){stop("Missing 'sampleFile' parameter.")}
  if(missing(genome)){stop("Missing 'genome' parameter.")}

  # create a qProject, perform various tests, install required BSgenome and load aligner package
  # create fasta indices (.fai) files for all the reference sequences (genome & aux) as well as all md5 sum files
  proj <- createQProject(sampleFile, genome, auxiliaryFile, aligner, maxHits, paired, splicedAlignment, snpFile, bisulfite, alignmentParameter, projectName, alignmentsDir, lib.loc, cacheDir)

  # display the status of the project
  missingFilesMessage(proj)

  # check if there are any genomic alignment files missing. if so, create index and align
  if(any(is.na(proj@alignments$FileName))){
    # build genome index
    if(is.na(proj@snpFile)){
      if(proj@genomeFormat=="file"){
        buildIndex(proj@genome, paste(proj@genome,proj@alnModeID,sep="."), proj@alnModeID, resolveCacheDir(proj@cacheDir))
      }else if(proj@genomeFormat=="BSgenome"){
        buildIndexPackage(proj@genome,proj@aligner,proj@alnModeID,resolveCacheDir(proj@cacheDir),proj@lib.loc)
      }else{stop("Fatal error 45906")}
    }else{
      # create indices for the two secondary genomes specified by snps
      buildIndexSNP(proj@snpFile,paste(proj@snpFile,proj@alnModeID,sep="."),proj@genome,proj@genomeFormat,proj@alnModeID,resolveCacheDir(proj@cacheDir))
    }

    # align to the genome (qProject gets upadated with the alignment file names)
    proj <- createGenomicAlignments(proj, clObj)

  }

  # create indices for the auxiliary files if necessary (if there are any alignment files missing)
  if(any(is.na(proj@auxAlignments))){
    for(i in 1:nrow(proj@aux)){
      buildIndex(proj@aux$FileName[i], paste(proj@aux$FileName[i],proj@alnModeID,sep="."), proj@alnModeID, resolveCacheDir(proj@cacheDir),checkMD5=TRUE)
    }

    proj <- createAuxAlignments(proj, clObj)
  }

  return(proj)
}




# inspect the qProject and give an overview of all the files that need to be computed
missingFilesMessage <- function(proj){

  genomeIndexNeedsToBeCreated <- FALSE
  genomicAlignmentsNeedToBeCreated <- sum(is.na(proj@alignments$FileName))
  auxAlignmentsNeedToBeCreated <- sum(is.na(proj@auxAlignments))

  if(any(is.na(proj@alignments$FileName))){
    if(is.na(proj@snpFile)){
      if(proj@genomeFormat=="file"){
        indexPath <- paste(proj@genome,proj@alnModeID,sep=".")
        if(file.exists(indexPath)){
          if(!("ref_md5Sum.txt" %in% dir(indexPath))){
            genomeIndexNeedsToBeCreated <- TRUE
          } 
        }else{
          genomeIndexNeedsToBeCreated <- TRUE
        }
      }else if(proj@genomeFormat=="BSgenome"){
       indexPackageName <- paste(proj@genome,proj@alnModeID,sep=".")
       if(!(indexPackageName %in% installed.packages())){
         genomeIndexNeedsToBeCreated <- TRUE
       }
      }else{stop("Fatal error 45906")}
    }else{
      indexPathR <- paste(proj@snpFile,basename(proj@genome),"R","fa",proj@alnModeID,sep=".")
      indexPathA <- paste(proj@snpFile,basename(proj@genome),"A","fa",proj@alnModeID,sep=".")

      if(file.exists(indexPathR) & file.exists(indexPathA)){
        if( !("ref_md5Sum.txt" %in% dir(indexPathR)) | !("ref_md5Sum.txt" %in% dir(indexPathA)) ){
          genomeIndexNeedsToBeCreated <- TRUE
        }
      }else{
        genomeIndexNeedsToBeCreated <- TRUE
      }
    }
  }

  if(!genomeIndexNeedsToBeCreated & genomicAlignmentsNeedToBeCreated==0 & auxAlignmentsNeedToBeCreated==0){
    message("all necessary alignment files found")
  }else{
    message("alignment files missing - need to:")
    if(genomeIndexNeedsToBeCreated){message("    create alignment index for the genome")}
    if(genomicAlignmentsNeedToBeCreated>0){message("    create ",genomicAlignmentsNeedToBeCreated," genomic alignment(s)")}
    if(auxAlignmentsNeedToBeCreated>0){message("    create ",auxAlignmentsNeedToBeCreated," auxiliary alignment(s)")}

    if(interactive()){
      message("will start in ", appendLF=FALSE); flush.console()
      for(i in 9:1) { message("..",i,"s",appendLF=FALSE); flush.console();
      Sys.sleep(1) }
      message("")
    }
  }
}



# read all the input information, perform various checks and compile the data into a qProject object
createQProject <- function(sampleFile, genome, auxiliaryFile, aligner, maxHits, paired, splicedAlignment, 
                           snpFile, bisulfite, alignmentParameter, projectName, alignmentsDir, 
                           lib.loc, cacheDir){

  # instantiate a qProject
  proj <-new("qProject")

  # -----------DETERMINE THE FORMAT OF THE GENOME AND DOWNLOAD FROM BIOCONDUCTOR IF NEEDED -------------

  if(length(genome) != 1){stop("The genome parameter should contain only a single entry",call.=FALSE);}

  # if there exists a / or \\ at the end of the genome string, remove it. this is required for the
  # consistent behavior of file.exists on windows systems. file.exists("dir/") != file.exists("dir").
  genome <- sub("(\\\\|/)$","",genome)

  if(file.exists(genome)){
    # genome is a directory or a file
    if(!(file.info(genome)$isdir)){
       # genome is a file, check if it is a fasta file
      if(consolidateFileExtensions(genome)=="fasta"){
        proj@genome <- tools::file_path_as_absolute(genome)
        proj@genomeFormat <- "file"
      }else{stop("The specified genome ",genome," does not have the extension of a fasta file (fa,fasta,fna)",call.=FALSE)}
    }else{
      # genome is a directory
      stop("The specified genome has to be a file and not a directory: ", genome,call.=FALSE)
    }
  }else{
    # check if the genome is a BSgenome.
    if(genome %in% installed.packages(lib.loc=lib.loc)[,'Package']){
      # BSgenome is already installed, load it
      if(!require(genome, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)){stop("The BSgenome ",genome," is installed but cannot be loaded. The version of the BSgenome might be too old",call.=FALSE)}

      proj@genome <- genome
      proj@genomeFormat <- "BSgenome"
    }else{
      message("The specified genome is not a fasta file or an installed BSgenome.")
      message("Connecting to Bioconductor and searching for a matching genome (internet connection required)...", appendLF = FALSE)

      if(testBioCConnection()){
        # Connection to Bioconductor OK
        message("OK")
        if(genome %in% available.genomes()){
          # The genome is available in Bioconductor, install it
          message("Downloading genome... ",genome," from Bioconductor")
          biocLite(genome, suppressUpdates=TRUE, lib=lib.loc)

          # BSgenome has been installed, load it
          if(!require(genome, character.only=TRUE, quietly=TRUE, lib.loc=lib.loc)){stop("Fatal error 23445")}

          proj@genome <- genome
          proj@genomeFormat <- "BSgenome"
        }else{
          # The genome is not available in Bioconductor
          stop(genome," is not available in Bioconductor. Type available.genomes() for a complete list",call.=FALSE)
        }
      }else{
        # No connection to Bioconductor
        message("FAILED")
        stop("Could not find the specified genome: ",genome,call.=FALSE)
      }
    }
  }


  # ---------------------------------------- PARSE THE PAIRED PARAMETER ---------------------------------
  if(!is.null(paired)){
    if(!(paired %in% c("no","fr","rf","ff"))){
      stop("The parameter paired supports only NULL, 'no', 'fr', 'rf' and 'ff'",call.=FALSE)
    }
  }

  # ---------------------------------------- PARSE THE SAMPLE FILE -------------------------------------
  if(!file.exists(sampleFile)){stop("Cannot open ",sampleFile,call.=FALSE)}
  samples <- read.delim(sampleFile,header=TRUE,as.is=TRUE)

  if(nrow(samples)==0){stop(sampleFile," is either empty or there is a header missing",call.=FALSE)}

  # get absolute path of the sample.txt file
  sampleFileAbsolutePath <- dirname(tools::file_path_as_absolute(sampleFile))

  # Check if the reads are single end
  if(ncol(samples)==2){
    if(all(colnames(samples) == c("FileName","SampleName"))){
      # this is a valid single end samples file. check if listed files exist
      for(i in 1:nrow(samples)){
        pathRet <- pathAsAbsoluteRedirected(samples$FileName[i],sampleFileAbsolutePath)
        if(!is.null(pathRet)){
          if(samples$SampleName[i]==""){stop(samples$FileName[i]," listed in ",sampleFile," has no sample name",call.=FALSE)}
          samples$FileName[i] <- pathRet
        }else{stop(samples$FileName[i]," listed in ",sampleFile," does not exist",call.=FALSE)}
      }
      # deterime format of the files (fa,fasta,fna)
      proj@samplesFormat <- determineSamplesFormat(samples$FileName)

      if(proj@samplesFormat == "bam"){
        # samples are provided as bam files
        if(is.null(paired)){stop("Paired status of the provided bam files cannot be determined automatically, please set the paired parameter",call.=FALSE)}

        if(paired=="no"){
          proj@reads <- samples
          proj@reads$FileName=NA_character_
        }else{
          proj@reads <- data.frame(FileName1=NA_character_,FileName2=NA_character_,SampleName=samples$SampleName,stringsAsFactors=FALSE)
        }

        proj@paired <- paired
        proj@alignments <- samples
      }else{
        # samples are provided as reads
        if(!is.null(paired)){
          if(paired!="no"){stop(sampleFile," contains 2 columns which represents unpaired end data. Please reset the paired parameter",call.=FALSE)}
        }
        if(length(unique(samples$FileName)) != nrow(samples)){
          stop("There are duplicate files in sampleFile. This is not allowed as it would result in non-unique alignment file names",call.=FALSE)
        }

        proj@paired <- "no"
        proj@reads <- samples
        proj@alignments <- samples
        proj@alignments$FileName=NA_character_
      }
    }else{stop(sampleFile," should contain column names FileName,SampleName (single end data)",call.=FALSE)} # incorrect format
  # Check if the reads are paired data
  }else if(ncol(samples)==3){
    if(all(colnames(samples) == c("FileName1","FileName2","SampleName"))){
      # this is a valid paired end samples file. check if listed files exist
      for(i in 1:nrow(samples)){
        pathRet <- pathAsAbsoluteRedirected(samples$FileName1[i],sampleFileAbsolutePath)
        if(!is.null(pathRet)){
          if(samples$SampleName[i]==""){stop(samples$FileName1[i]," listed in ",sampleFile," has no sample name",call.=FALSE)}
          samples$FileName1[i] <- pathRet
        }else{stop(samples$FileName1[i]," listed in ",sampleFile," does not exist",call.=FALSE)}
      }
      for(i in 1:nrow(samples)){
        pathRet <- pathAsAbsoluteRedirected(samples$FileName2[i],sampleFileAbsolutePath)
        if(!is.null(pathRet)){
          if(samples$SampleName[i]==""){stop(samples$FileName2[i]," listed in ",sampleFile," has no sample name",call.=FALSE)}

          samples$FileName2[i] <- pathRet
        }else{stop(samples$FileName2[i]," listed in ",sampleFile," does not exist",call.=FALSE)}
      }
      # deterime format of the files (fa,fasta,fna)
      proj@samplesFormat <- determineSamplesFormat(c(samples$FileName1,samples$FileName2))

      if(proj@samplesFormat == "bam"){stop("Bam files need to be listed in a two column file: ",sampleFile,call.=FALSE)}
      if(!is.null(paired)){
        if(paired=="no"){stop(sampleFile," contains 3 columns which represents paired end data. Please reset the paired parameter",call.=FALSE)}
      }
      if(length(unique(samples$FileName1)) != nrow(samples)){
        stop("There are duplicate files in sampleFile. This is not allowed as it would result in non-unique alignment file names",call.=FALSE)
      }

      if(is.null(paired)){
        proj@paired <- "fr"
      }else{
        proj@paired <- paired
      }

      proj@reads <- samples
      proj@alignments <- data.frame(FileName=NA_character_,SampleName=samples$SampleName,stringsAsFactors=FALSE)

    }else{stop(sampleFile," should contain the column names FileName1,FileName2,SampleName (paired end data)",call.=FALSE)} # incorrect format
  # The format of the file is not supported
  }else{stop(sampleFile," should be a tab delimited file with either 2 or 3 columns",call.=FALSE)}


  # --------------- FOR FASTQ FILES, READ A SMALL CHUNK TO GUESS THE QUALITY FORMAT (phred33 or phred64) -------
  if(proj@samplesFormat == "fastq"){
    proj@reads <- data.frame(proj@reads,phred=NA_character_,stringsAsFactors=FALSE) # add an additional quality column to the reads table
    for(i in 1:nrow(proj@reads)){
      fastq_fs <- FastqStreamer(proj@reads[i,1], n=2000,readerBlockSize=5e5); # take only first column even for paired end data

      tryCatch({
        fastq_sq <- yield(fastq_fs)
      }, error = function(ex) {
        close(fastq_fs)
        stop("Incorrect format for the fastq file: ",proj@reads[i,1],call.=FALSE)
      })
      close(fastq_fs)
      # determine the format of the quality values
      if(inherits(quality(fastq_sq),"FastqQuality")){
        proj@reads$phred[i] <- "33"
      }else if(inherits(quality(fastq_sq),"SFastqQuality")){
        proj@reads$phred[i] <- "64"
      }else{stop("The quality values of the provided sequences files are not interpretable: ",proj@reads[i,1],call.=FALSE)}
    }
  }


  # -------------------- CALCULATE MD5 SUB SUMS FOR ALL READ FILES ----------------------------------------
  if(proj@samplesFormat != "bam"){
    if(proj@paired=="no"){
      proj@reads_md5subsum=data.frame(md5subsum=md5subsum(proj@reads$FileName),stringsAsFactors=FALSE)
    }else{
      proj@reads_md5subsum=data.frame(md5subsum1=md5subsum(proj@reads$FileName1),md5subsum2=md5subsum(proj@reads$FileName2),stringsAsFactors=FALSE)
    }
  }else{
    if(proj@paired=="no"){
      proj@reads_md5subsum=as.data.frame(matrix(NA_character_,nrow=nrow(proj@reads),ncol=1),stringsAsFactors=FALSE)
      colnames(proj@reads_md5subsum) <- "md5subsum"
    }else{
      proj@reads_md5subsum=as.data.frame(matrix(NA_character_,nrow=nrow(proj@reads),ncol=2),stringsAsFactors=FALSE)
      colnames(proj@reads_md5subsum) <- c("md5subsum1","md5subsum2")
    }
  }

  # ----------------------------------- PARSE THE SNP FILE ------------------------------------  

  if(!is.null(snpFile)){
    if(!file.exists(snpFile)){stop("Cannot open ",snpFile,call.=FALSE)}
    proj@snpFile=tools::file_path_as_absolute(snpFile)
  }else{proj@snpFile=NA_character_}


  # ---------------------------------------- PARSE THE AUXILIARY FILE -------------------------------------
  if(!is.null(auxiliaryFile)){
    if(proj@samplesFormat == "bam"){stop("The option 'auxiliaryFile' is not supported for bam input files",call.=FALSE)}
    if(!is.na(proj@snpFile)){stop("The option 'auxiliaryFile' is not supported in allelic mode",call.=FALSE)}
    if(!file.exists(auxiliaryFile)){stop("Cannot open ",auxiliaryFile,call.=FALSE)}
    auxiliaries <- read.delim(auxiliaryFile,header=TRUE,as.is=TRUE)

    if(nrow(auxiliaries)==0){stop(auxiliaryFile," is either empty or there is a header missing",call.=FALSE)}

    # get absolute path of the auxiliary.txt file
    auxFileAbsolutePath <- dirname(tools::file_path_as_absolute(auxiliaryFile))

    if(ncol(auxiliaries)==2){
      if(all(colnames(auxiliaries) == c("FileName","AuxName"))){
        # this is a valid auxiliaries file. check if listed files exist
        for(i in 1:nrow(auxiliaries)){
          pathRet <- pathAsAbsoluteRedirected(auxiliaries$FileName[i],auxFileAbsolutePath)
          if(!is.null(pathRet)){
            if(auxiliaries$AuxName[i]==""){stop(auxiliaries$FileName[i]," listed in ",auxiliaryFile," has no name",call.=FALSE)}
            auxiliaries$FileName[i] <- pathRet
          }else{stop(auxiliaries$FileName[i]," listed in ",auxiliaryFile," does not exist",call.=FALSE)}
        }
        # test the file extensions
        allAuxFileExts <- unique(consolidateFileExtensions(auxiliaries$FileName))
        if(!all(consolidateFileExtensions(auxiliaries$FileName)=="fasta")){stop("All auxiliary files need to have a fasta extension (fa,fasta,fna)",call.=FALSE)}

        if(length(unique(auxiliaries$AuxName))!=nrow(auxiliaries)){stop("All auxiliary files need to have a unique AuxName",call.=FALSE)}

        # store the auxiliary files
        proj@aux <- auxiliaries
      }else{stop(auxiliaryFile," should contain column names FileName,AuxName",call.=FALSE)} # incorrect format
    }else{stop(auxiliaryFile," should be a tab delimited file with 2 columns",call.=FALSE)}

    # Fill data into auxAlignments
    auxAlignments <- data.frame(matrix(NA_character_,nrow=nrow(auxiliaries),ncol=nrow(proj@reads)),stringsAsFactors=FALSE)
    rownames(auxAlignments) <- auxiliaries$AuxName
    colnames(auxAlignments) <- proj@reads$SampleName
    proj@auxAlignments <- auxAlignments

  }



  # ----------------------------------- PARSE BISULFITE PARAMETER -----------------------------------

  if(bisulfite %in% c("no","dir","undir")){
    proj@bisulfite <- bisulfite
    if(bisulfite != "no" && !(proj@paired %in% c("no","fr"))){
      stop("Bisulfite paired-end mode only supports pair orientation 'fr'",call.=FALSE)
    }
  }else{stop("Bisulfite mode only supports 'no', 'dir' and 'undir'",call.=FALSE)}
  
  # ----------------------------------- PARSE MAXHITS PARAMETER --------------------------------------
  proj@maxHits <- maxHits

  #------------------------------------ PARSE SPLICED ALIGNMENT PARAMETER ---------------------------
  proj@splicedAlignment <- splicedAlignment
  if(proj@splicedAlignment & (proj@bisulfite!="no")){stop("The spliced alignment mode is not supported for bisulfite samples")}
  if(proj@splicedAlignment & !(proj@paired %in% c("no","fr"))){stop("The spliced alignment mode only supports the pair orientation 'fr'")}

  #------------------------------------ PARSE THE ALIGNMENT PRAMETERS ----------------------------------
  if(is.null(alignmentParameter)){
    if(!proj@splicedAlignment){
      # Test for the case where no merge reorder is going to be executed later on. In that case maxhits needs to 
      # to be reinforced by bowtie.
      if((proj@bisulfite == "no") && (is.na(proj@snpFile))){
        proj@alignmentParameter <- paste("-m",proj@maxHits,"--best --strata")
      }else{
        proj@alignmentParameter <- paste("-k",proj@maxHits+1,"--best --strata")
      }
      # For the allelic case, ignore qualities. Anyways the assignment to the R or A genome is based on sequence.
      if((proj@samplesFormat == "fasta") || !is.null(snpFile)){ 
        proj@alignmentParameter <- paste(proj@alignmentParameter,"-v 2")
      }
      if(proj@paired != "no"){
        proj@alignmentParameter <- paste(proj@alignmentParameter,"--maxins 500")
      }

    }else{
      if(is.na(proj@snpFile)){
        proj@alignmentParameter <- "-max_intron 400000 -min_intron 20000 -max_multi_hit 10 -selectSingleHit TRUE -seed_mismatch 1 -read_mismatch 2 -try_hard yes"
      }else{
        proj@alignmentParameter <- "-max_intron 400000 -min_intron 20000 -max_multi_hit 10 -selectSingleHit FALSE -seed_mismatch 1 -read_mismatch 2 -try_hard yes"
      }
    }
  }else{
    if(length(alignmentParameter)==1){
      proj@alignmentParameter <- alignmentParameter
    }else{stop("The alignmentParameter should contain only a single character string",call.=FALSE)}
  }

  #------------------------------------ PARSE THE PROJECT NAME ----------------------------------
  if(is.null(projectName)){
    proj@projectName <- NA_character_
  }else{
    if(length(projectName)==1){
      proj@projectName <- projectName
    }else{stop("The projectName should contain only a single entry",call.=FALSE)}
  }

  # ---------------------------------- PARSE THE BAMFILE DIRECTORY --------------------------
  if(is.null(alignmentsDir)){
    proj@alignmentsDir <- NA_character_
  }else{
     alignmentsDir <- sub("(\\\\|/)$","",alignmentsDir)  # for windows systems
    if(file.exists(alignmentsDir)){
      # it's a directory or a file
      if(file.info(alignmentsDir)$isdir){
        proj@alignmentsDir <- tools::file_path_as_absolute(alignmentsDir)
      }else{stop("alignmentsDir ",alignmentsDir, " is not a directory",call.=FALSE)}
    }else{stop("alignmentsDir ",alignmentsDir, " does not exist",call.=FALSE)}
  }

  # ---------------------------------- PARSE THE LIB.LOC DIRECTORY --------------------------
  if(is.null(lib.loc)){
    proj@lib.loc <- NA_character_
  }else{
     lib.loc <- sub("(\\\\|/)$","",lib.loc)  # for windows systems
    if(file.exists(lib.loc)){
      # it's a directory or a file
      if(file.info(lib.loc)$isdir){
        proj@lib.loc <- tools::file_path_as_absolute(lib.loc)
      }else{stop("lib.loc ",lib.loc, " is not a directory",call.=FALSE)}
    }else{stop("lib.loc ",lib.loc, " does not exist",call.=FALSE)}
  }

  # ---------------------------------- PARSE THE CACHE DIRECTORY --------------------------

  if(is.null(cacheDir)){
    proj@cacheDir <- NA_character_
  }else{
    cacheDir <- sub("(\\\\|/)$","",cacheDir)  # for windows systems
    if(file.exists(cacheDir)){
      # it's a directory or a file
      if(file.info(cacheDir)$isdir){
        # if not already present, create a subdirectory in the cacheDir, this prevent collisions 
        # when the the same cacheDir is being used by multiple users
        proj@cacheDir <- file.path(tools::file_path_as_absolute(cacheDir),basename(tempdir()))

        if(!file.exists(proj@cacheDir)){
          if(!dir.create(path=proj@cacheDir, showWarnings=FALSE)){
            stop("No permissions to write in the cacheDir: ",proj@cacheDir,call.=FALSE)
          }
        }
      }else{stop("cacheDir ",cacheDir, " is not a directory",call.=FALSE)}
    }else{stop("cacheDir ",cacheDir, " does not exist",call.=FALSE)}
  }


  # --------------- SET THE ALIGNERMODE ID AND LOAD THE ALIGNER PACKAGE IF REQUIRED --------------------

  supportedAligners <- c("Rbowtie")
  if(aligner %in% supportedAligners){
    proj@aligner <- aligner
  }else{stop("The specified aligner is unknown, please select one of the following: ",paste(supportedAligners,collapse=","),call.=FALSE)}

  if((proj@samplesFormat %in% c("fasta","fastq")) & (proj@bisulfite=="no")){
    alnModeID <- proj@aligner
  }else if((proj@samplesFormat %in% c("fasta","fastq")) & (proj@bisulfite!="no")){
    alnModeID <- paste(proj@aligner,"CtoT",sep="")
  }else if((proj@samplesFormat %in% c("csfasta","csfastq")) & proj@bisulfite=="no"){
    alnModeID <- paste(proj@aligner,"Cs",sep="")
  }else if((proj@samplesFormat %in% c("csfasta","csfastq")) & proj@bisulfite!="no"){
    stop("Bisulfite alignment mode is not available for color space reads")
  }else if(proj@samplesFormat == "bam"){
    alnModeID <- NA_character_
    proj@aligner <- NA_character_
  }else{stop("Fatal error 2340923")}

  proj@alnModeID <- alnModeID

  # load the aligner package in the case where there are missing alignments (genomic or aux)
  if(any(is.na(proj@alignments$FileName)) | any(is.na(proj@auxAlignments))){

    pkgname <- aligner
    
    # these test are not needed while Rbowtie is in "Depends"
    #if(!(pkgname %in% installed.packages())){stop(pkgname, " package is required for the alignments but not installed on this system",call.=FALSE)}
    #if(!require(pkgname, character.only=TRUE, quietly=TRUE)){stop("Fatal error 340954")}
    #pkg_version <- installed.packages()[pkgname, 'Version']
    #QuasR_suggests <- unlist(strsplit(installed.packages()['QuasR','Suggests'], ","))
    #QuasR_suggests_pkg <- grep(pkgname, QuasR_suggests, value=T)
  }


  # ---------------------------------- CREATE .fai AND .md5 FILES --------------------------

  # create fasta indices (.fai) files for all the reference sequences (genome & aux)
  # create all the .md5 files necessary to identify preexisting bam files
  createReferenceSequenceIndices(proj)


  # --------------------- SEARCH FOR BAM FILES THAT HAVE BEEN CREATED PREVIOUSLY ----------------------

  # search for genomic alignments
  for(i in 1:nrow(proj@reads)){
    if(is.na(proj@alignments$FileName[i])){
       projBamInfo <- bamInfoOnlyBaseName(qProjectBamInfo(proj,i))
       if(is.na(proj@alignmentsDir)){bamDir <- dirname(proj@reads[i,1])}else{bamDir <- proj@alignmentsDir}
       samplePrefix <- basename(tools::file_path_sans_ext(proj@reads[i,1],compression=TRUE))
       filesInBamDir <- list.files(bamDir)
       bamFilesToInspect <- filesInBamDir[nchar(sub(paste(samplePrefix,"\\_[^\\_]+.bam$",sep=""),"",filesInBamDir))==0]
       bamTxtFilesToInspect <- paste(file.path(bamDir,bamFilesToInspect),"txt",sep=".")
       bamTxtFilesToInspectExist <- bamTxtFilesToInspect[file.exists(bamTxtFilesToInspect)]
       bamFilesToInspectWithTxt <- file.path(bamDir,bamFilesToInspect[file.exists(bamTxtFilesToInspect)])
       
       compatibleBamFileInd=NULL
       if(length(bamTxtFilesToInspectExist)>0){
         for(m in 1:length(bamTxtFilesToInspectExist)){
           bamInfoT_DF <- read.delim(bamTxtFilesToInspectExist[m],header=FALSE,row.names=1,stringsAsFactors=FALSE)
           bamInfoT <- bamInfoT_DF[,1]
           names(bamInfoT) <- rownames(bamInfoT_DF)
           bamInfoT <- bamInfoOnlyBaseName(bamInfoT)

           # compare the actual parameters to the one from the bam file on disk
           if(identical(projBamInfo,bamInfoT)){
             compatibleBamFileInd[length(compatibleBamFileInd)+1] <- m
           }
         }
         if(length(compatibleBamFileInd)>1){
           for(k in 1:length(compatibleBamFileInd)){message(bamFilesToInspectWithTxt[k])}
           stop("Multiple bam files exist with same alignment parameters (see above list). QuasR is unable to decide which one to use. Please delete manually the respective bam files",call.=FALSE)
         }
         if(length(compatibleBamFileInd)==1){
           proj@alignments$FileName[i] <- bamFilesToInspectWithTxt[compatibleBamFileInd]
         }
       }
    }
  }

  # search for aux alignments
  if(nrow(proj@auxAlignments)>0){
    for(i in 1:ncol(proj@auxAlignments)){
      for(j in 1:nrow(proj@auxAlignments)){
        projBamInfo <- bamInfoOnlyBaseName(qProjectBamInfo(proj,i,j))
          if(is.na(proj@auxAlignments[j,i])){
           if(is.na(proj@alignmentsDir)){bamDir <- dirname(proj@reads[i,1])}else{bamDir <- proj@alignmentsDir}
           samplePrefix <- basename(tools::file_path_sans_ext(proj@reads[i,1],compression=TRUE))
           filesInBamDir <- list.files(bamDir)
           bamFilesToInspect <- filesInBamDir[nchar(sub(paste(samplePrefix,"\\_[^\\_]+.bam$",sep=""),"",filesInBamDir))==0]
           bamTxtFilesToInspect <- paste(file.path(bamDir,bamFilesToInspect),"txt",sep=".")
           bamTxtFilesToInspectExist <- bamTxtFilesToInspect[file.exists(bamTxtFilesToInspect)]
           bamFilesToInspectWithTxt <- file.path(bamDir,bamFilesToInspect[file.exists(bamTxtFilesToInspect)])

           compatibleBamFileInd=NULL
           if(length(bamTxtFilesToInspectExist)>0){
             for(m in 1:length(bamTxtFilesToInspectExist)){
               bamInfoT_DF <- read.delim(bamTxtFilesToInspectExist[m],header=FALSE,row.names=1,stringsAsFactors=FALSE)
               bamInfoT <- bamInfoT_DF[,1]
               names(bamInfoT) <- rownames(bamInfoT_DF)
               bamInfoT <- bamInfoOnlyBaseName(bamInfoT)

               # compare the actual parameters to the one from the bam file on disk
               if(identical(projBamInfo,bamInfoT)){
                 compatibleBamFileInd[length(compatibleBamFileInd)+1] <- m
               }
             }
             if(length(compatibleBamFileInd)>1){
               for(k in 1:length(compatibleBamFileInd)){message(bamFilesToInspectWithTxt[k])}
               stop("Multiple bam files exist with same alignment parameters (see above list). QuasR is unable to decide which one to use. Please delete manually the respective bam files",call.=FALSE)
             }
             if(length(compatibleBamFileInd)==1){
               proj@auxAlignments[j,i] <- bamFilesToInspectWithTxt[compatibleBamFileInd]
             }
           }
         }
      }
    }
  }


  return(proj)
}



# create search index (.fai) files for all the reference sequences that are
# going to be used for the alignments (genome & aux)
# do nothing if the fai files exist already
# create an .md5 file that contains the md5sum of the referece sequences
# if snp file exists, calculate md5 and store it in a (.md5) file
createReferenceSequenceIndices <- function(proj){
  # create fasta index for the genome if it is not a BSgenome and if the .fai file is not present
  if(proj@genomeFormat=="file"){
    if(!file.exists(paste(proj@genome,"fai",sep="."))){
      # test if the sequence file is a valid fasta file using the fasta.info() function from biostrings
      # this is necessary because indexFa() from Rsamtools would crash R if it is passed an invalid fasta file
      message(paste("Creating .fai file for:",proj@genome))
      result <- try(fasta.info(proj@genome))
      if(class(result)=="try-error"){
        stop("The fasta file ",proj@genome," is not a valid fasta file",call.=FALSE)
      }
      faInfoNames <- sub("\\s+.+","",names(result),perl=TRUE) # remove everything after the white space after the id (not done by fasta.info)
      if(length(unique(faInfoNames)) != length(faInfoNames)){
        stop("Sequence names in the file: ",proj@genome," are not unique",call.=FALSE)
      }
      if(class(try(indexFa(proj@genome)))=="try-error"){
        stop("Cannot write into the directory where ",proj@genome," is located. Make sure you have the right permissions",call.=FALSE)
      }
    }else{
      faiSeqNames <- as.character(seqnames(scanFaIndex(proj@genome)))
      if(length(unique(faiSeqNames)) != length(faiSeqNames)){
        stop("Sequence names in the file: ",proj@genome," are not unique (information extracted from the .fai file",call.=FALSE)
      }
    }
    # create md5 sum file
    if(!file.exists(paste(proj@genome,"md5",sep="."))){
      write.table(tools::md5sum(proj@genome),paste(proj@genome,"md5",sep="."),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
  }
  # create fasta index for the auxilliaries if not present already
  if(nrow(proj@aux)>0){
    for(i in 1:nrow(proj@aux)){
      if(!file.exists(paste(proj@aux$FileName[i],"fai",sep="."))){
        # test if the sequence file is a valid fasta file using the fasta.info() function from biostrings
        # this is necessary because indexFa() from Rsamtools would crash R if it is passed an invalid fasta file (R version 2.14.1)
        result <- try(fasta.info(proj@aux$FileName[i]))
        if(class(result)=="try-error"){
          stop("The fasta file ",proj@aux$FileName[i]," is not a valid fasta file",call.=FALSE)
        }
        faInfoNames <- sub("\\s+.+","",names(result),perl=TRUE) # remove everything after the white space after the id (not done by fasta.info)
        if(length(unique(faInfoNames)) != length(faInfoNames)){
          stop("Sequence names in the file: ",proj@aux$FileName[i]," are not unique",call.=FALSE)
        }
        if(class(try(indexFa(proj@aux$FileName[i])))=="try-error"){
          stop("Cannot write into the directory where ",proj@aux$FileName[i]," is located. Make sure you have the right permissions",call.=FALSE)
        }
      }

      if(!file.exists(paste(proj@aux$FileName[i],"md5",sep="."))){
        # create md5 sum file
        write.table(tools::md5sum(proj@aux$FileName[i]),paste(proj@aux$FileName[i],"md5",sep="."),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
      }
    }
  }

  #create md5 sum file for the snp file
  if(!is.na(proj@snpFile)){
    if(!file.exists(paste(proj@snpFile,"md5",sep="."))){
      write.table(tools::md5sum(proj@snpFile),paste(proj@snpFile,"md5",sep="."),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
  }
}


# returns information for one bam file in qProject
# sampleNr specifies the sample for which the information is compiled
# if auxNr is not NULL, it returns information about aux bam file
qProjectBamInfo <- function(proj,sampleNr,auxNr=NULL){

  if(sampleNr > nrow(proj@reads)){stop("project has no sample ",sampleNr)}
  if(!is.null(auxNr)){
    if(auxNr > nrow(proj@aux)){stop("project has no aux ",auxNr)}
  }

  alnInfo <- NULL
  if(proj@paired=="no"){
    alnInfo["reads.FileName1"]=proj@reads$FileName[sampleNr]
    alnInfo["reads.FileName2"]=NA
    alnInfo["reads.md5subsum1"]=proj@reads_md5subsum$md5subsum[sampleNr]
    alnInfo["reads.md5subsum2"]=NA
  }else{
    alnInfo["reads.FileName1"]=proj@reads$FileName1[sampleNr]
    alnInfo["reads.FileName2"]=proj@reads$FileName2[sampleNr]
    alnInfo["reads.md5subsum1"]=proj@reads_md5subsum$md5subsum1[sampleNr]
    alnInfo["reads.md5subsum2"]=proj@reads_md5subsum$md5subsum2[sampleNr]
  }

  alnInfo["samplesFormat"]=proj@samplesFormat
  alnInfo["genome"]=proj@genome
  if(proj@genomeFormat=="file"){
    alnInfo["genome.md5"]=as.character(read.delim(paste(proj@genome,"md5",sep="."),header=FALSE)[1,1])
  }else{
    alnInfo["genome.md5"]=NA
  }
  alnInfo["genomeFormat"]=proj@genomeFormat
  alnInfo["aux"]=NA
  alnInfo["aux.md5"]=NA
  alnInfo["aligner"]=proj@aligner
  if(!is.na(proj@aligner)){
    alnInfo["aligner.version"]=installed.packages()[proj@aligner, 'Version']
  }else{
    alnInfo["aligner.version"]=NA
  }
  alnInfo["maxHits"]=proj@maxHits
  alnInfo["paired"]=proj@paired
  alnInfo["splicedAlignment"]=proj@splicedAlignment
  alnInfo["snpFile"]=proj@snpFile
  alnInfo["snpFile.md5"]=NA
  alnInfo["bisulfite"]=proj@bisulfite
  alnInfo["alignmentParameter"]=proj@alignmentParameter
  alnInfo["QuasR.version"]=installed.packages()['QuasR', 'Version']

  if(!is.null(auxNr)){
    alnInfo["aux"]=proj@aux$FileName[auxNr]
    alnInfo["aux.md5"]=as.character(read.delim(paste(proj@aux$FileName[auxNr],"md5",sep="."),header=FALSE)[1,1])
  }

  if(!is.na(proj@snpFile)){
    alnInfo["snpFile.md5"]=as.character(read.delim(paste(proj@snpFile,"md5",sep="."),header=FALSE)[1,1])
  }

  return(alnInfo)
}


# replace full file paths in bamInfo by base name
# remove the aligner and QuasR version. this allows using bam files generated with an older QuasR version
bamInfoOnlyBaseName <- function(bamInfo){
  bamInfo["reads.FileName1"] <- basename(bamInfo["reads.FileName1"])
  bamInfo["reads.FileName2"] <- basename(bamInfo["reads.FileName2"])
  bamInfo["genome"] <- basename(bamInfo["genome"])
  bamInfo["aux"] <- basename(bamInfo["aux"])
  bamInfo["snpFile"] <- basename(bamInfo["snpFile"])

  if("aligner.version" %in% names(bamInfo)){
    bamInfo <- bamInfo[!(names(bamInfo) %in% "aligner.version")]
  }

  if("QuasR.version" %in% names(bamInfo)){
    bamInfo <- bamInfo[!(names(bamInfo) %in% "QuasR.version")]
  }

  return(bamInfo)
}

# helper function that converts filepaths to absolute filepaths in the case where the working
# directory needs to be temporarily redirected. This is needed if e.g. the samples.txt file
# is not located in the working directory.
# returns the absolute path if the file exists
# returns NULL if the file does not exist
pathAsAbsoluteRedirected <- function(fileName,redirectPath){
  curwd <- getwd() # store the original working directory required to jump back
  on.exit(setwd(curwd)) # make sure that it will jump back in a case of error or not

  setwd(redirectPath) # jump to the redirected position
  if(file.exists(fileName)){
    return(tools::file_path_as_absolute(fileName))
  }else{return(NULL);}
}

# helper function that returns the temporary directory given the cacheDir (from qProject)
# if the user specified cacheDir, then it returns it. but if the user did not specify a cacheDir
# it returns the R temporary directory (which can be different for different machines)
resolveCacheDir <- function(cacheDir){
  if(is.na(cacheDir))
    return(tempdir())
  else
    return(cacheDir)
}


# given a vector of filenames, the function determines the final format 
# taking into account all the different types of extensions for 
# sequence files. Throughs an error if one ore more samples do not contain 
# a valid file extension
determineSamplesFormat <- function(filename){
  fileExtension <- consolidateFileExtensions(filename,compressed=TRUE) 

  validExtsSel <- fileExtension %in% c("fasta","fastq","bam","csfasta","csfastq")

  if(all(validExtsSel)){
    if(length(unique(fileExtension))==1){
      return(unique(fileExtension))
    }else{stop("Mixed sample file formats are not allowed: ",paste(unique(fileExtension),collapse=","),call.=FALSE)}
  }else{
    stop(filename[!validExtsSel][1]," does not have a valid file extension (fa,fna,fasta,fq,fastq,bam,csfasta,csfastq)",call.=FALSE)
  } 
}


