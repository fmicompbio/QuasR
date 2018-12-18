### helper functions used by QuasR's unit test functions



createProjectAllelic <- function(){
    require(Rsamtools)

    # Create samfile with XV tag
    samfile <- tempfile()
    cat("@HD\tVN:1.0\tSO:unsorted\n", file=samfile, append=FALSE)
    cat("@SQ\tSN:chrV\tLN:99\n", file=samfile, append=T)
    # fr XV=R
    cat("seq1\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:R\n", file=samfile, append=T)
    cat("seq1\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:R\n", file=samfile, append=T)
    # fr XV=U
    cat("seq2\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:U\n", file=samfile, append=T)
    cat("seq2\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:U\n", file=samfile, append=T)
    # fr XV=A
    cat("seq3\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\tXV:A:A\n", file=samfile, append=T)
    cat("seq3\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\tXV:A:A", file=samfile, append=T)

    # Create Bamfile
    bamfile <- asBam(samfile, samfile)
    
    # Create Sample File
    samplefile <- tempfile()
    write.table(data.frame(FileName=bamfile,
                           SampleName="Allelic", stringsAsFactors=F),
                sep="\t", quote=F, row.names=F, file=samplefile)
    
    # Create Genome
    genome <- tempfile(fileext=".fa")
    cat(">chrV\n", file=genome, append=FALSE)
    cat(paste(rep("G",99), collapse=""), file=genome, append=T)
    
    # Create SNP file
    snp <- tempfile(fileext=".txt")
    cat(">chrV\t10\tC\tG", file=snp, append=FALSE)     
    
    # 
    td <- tempdir()
    project <- qAlign(samplefile, genome, snpFile=snp, paired="fr", alignmentsDir=td)
    
    return(project)
}



createGtfRegion <- function()
{
    require("rtracklayer")
    annotationFile <- system.file("extdata", "hg19sub_annotation.gtf",
                                  package="QuasR")
    gtfRegion <- import.gff(annotationFile, format="gtf",
                            feature.type="exon")
    names(gtfRegion) <- mcols(gtfRegion)$gene_name
    return(gtfRegion)
}

createGenomeRegion <- function()
{
    genomeFile <- system.file("extdata", "hg19sub.fa", package="QuasR")
    genomeRegion <- scanFaIndex(genomeFile)
    return(genomeRegion)
}

createAuxRegion <- function()
{
    auxGenomeFile <- system.file("extdata", "NC_001422.1.fa", package="QuasR")
    auxRegion <- scanFaIndex(auxGenomeFile)
    return(auxRegion)
}

createReads <- function(genomeFile, destDir, paired=TRUE, snpFile=NULL, bisulfite=NULL, transcript=NULL, format="fastq"){
    library(Rsamtools)
    genome <- scanFa(genomeFile)
    destination <- tempfile(tmpdir=destDir)
    filename1 <- ifelse(format=="fastq", 
                        paste(destination,"_1.fq", sep=""),
                        paste(destination,"_1.fa", sep=""))
    filename2 <- ifelse(format=="fastq", 
                        paste(destination,"_2.fq", sep=""),
                        paste(destination,"_2.fa", sep=""))
    if(paired){
        suppressWarnings(file.remove(c(filename1, filename2)))
    }else{
        suppressWarnings(file.remove(filename1))
    }
    
    if(!is.null(snpFile)){
        snps <- read.table(snpFile, colClasses=c("factor","numeric","character","character"))
        snps[,3] <- toupper(snps[,3]) # convert ref nucleotides to upper case if not already the case
        snps[,4] <- toupper(snps[,4]) # convert alt nucleotides to upper case if not already the case
        genomeSNP <- lapply(split(snps, snps[,1]), function(s){
            chr <- as.character(s[1,1])
            as.character(replaceLetterAt(genome[[chr]], s[,2], s[,4]))
        })
        genomeSNP <- DNAStringSet(unlist(genomeSNP, use.names=T))
        names(genomeSNP) <- paste(names(genomeSNP), "A", sep="-")
        names(genome) <- paste(names(genome), "R", sep="-")
        genome <- append(genome, genomeSNP)
    } else {
        names(genome) <- paste(names(genome), "R", sep="-")
    }
    
    if(!is.null(bisulfite)){
        CpG <- lapply(seq_along(genome), function(i){
            seq <- genome[[i]]
            seqnames <- unlist(strsplit(names(genome[i]), "-"))[1]
            ranges <- ranges(bisulfite[seqnames(bisulfite) == seqnames])
            mask <- Mask(mask.width=length(seq), start=start(ranges), end=end(ranges))
            masks(seq) <- mask
            start(matchPattern("CG", seq))
        })
        names(CpG) <- names(genome)
        
        genomeCtoT <- chartr("C","T", genome)
        genomeGtoA <- chartr("G","A", genome)
        genome <- lapply(seq_along(CpG), function(i){
            seqCtoT <- genomeCtoT[[names(CpG)[i]]]
            seqCtoT[CpG[[i]]] <- "C"
            seqCtoT <- as(seqCtoT, "XStringSet")
            names(seqCtoT) <- paste(names(CpG)[i], "CtoT", sep="-")
            
            seqGtoA <- genomeGtoA[[names(CpG)[i]]]
            seqGtoA[CpG[[i]]+1] <- "G"
            seqGtoA <- as(seqGtoA, "XStringSet")
            names(seqGtoA) <- paste(names(CpG)[i], "GtoA", sep="-")
            
            append(seqCtoT, seqGtoA)                        
        })
        genome <- do.call(c, genome)
    }
    
    set.seed(111)   
    sampleRange <- 5L
    width=50L + sampleRange
    gap=10L #-width+coverage
    insertSize <- 90L
    append <- FALSE
    
    lapply(seq_along(genome), function(i){
        nreads=floor(length(genome[[i]])/(width+gap))
        
        # create read ranges
        r1 <- successiveIRanges(rep(width, nreads), gapwidth=gap)
        start(r1) <- start(r1) + sample(length(r1))%%sampleRange
        r2 <- successiveIRanges(rep(width, nreads), gapwidth=gap, from=width+insertSize)
        inGenome <- end(r2) < length(genome[[i]])
        # extract and write sequence
        read1 <- Views(genome[[i]], r1[inGenome])
        read2 <- Views(genome[[i]], r2[inGenome])
        ids <- paste("seq", names(genome)[i], start(read1), end(read1), start(read2), end(read2), sep="_")
        read1 <- as(read1,"XStringSet")
        read2 <- reverseComplement(as(read2,"XStringSet"))
        # swap half of the pairs
        if(!is.null(bisulfite)){
            names(read1) <- paste(ids, "b", "l", sep="_")
            names(read2) <- paste(ids, "b", "r", sep="_")
        }        
        swap <- sample(length(read1))%%2
        sread1 <- append(read1[swap!=1], read2[swap==1])
        sread2 <- append(read2[swap!=1], read1[swap==1])
        sids <- c(ids[swap!=1], ids[swap==1])
        # add qname
        if(is.null(bisulfite)){
            names(sread1) <- paste(sids, swap, "/1", sep="_")
            names(sread2) <- paste(sids, swap, "/2", sep="_")
        }
        if(paired){
            writeXStringSet(sread1, filename1, append=TRUE, format=format)
            writeXStringSet(sread2, filename2, append=TRUE, format=format)
        }else{
            writeXStringSet(sread1, filename1, append=TRUE, format=format)
            writeXStringSet(sread2, filename1, append=TRUE, format=format)            
        }
    })
    
    samplefile <- paste(destination, "sample.txt", sep="")
    
    if(!paired){
        write.table(data.frame(FileName=filename1,
                               SampleName="Normal", stringsAsFactors=F),
                    sep="\t", quote=F, row.names=F, file=samplefile)        
    } else {
        write.table(data.frame(FileName1=filename1,
                               FileName2=filename2,
                               SampleName="Normal", stringsAsFactors=F),
                    sep="\t", quote=F, row.names=F, file=samplefile)        
    }
    
    return(samplefile)
}
