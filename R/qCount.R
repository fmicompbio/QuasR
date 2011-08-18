qCount <- function(qProject, gRange, stranded=FALSE, ...)
{

#     ## count of regions
    if(missing(gRange)){
        isGTFFormat <- .fileExtension(qProject@annotations$filepath) %in% c("gtf")
        gtfFiles <- qProject@annotations[isGTFFormat,]$filepath
        if(length(gtfFiles) != 1)
            stop("There is more/less than one 'gtf' file in the annotation")
        require(rtracklayer)
        #targets <- import.gff2("data/Drosophila_melanogaster.BDGP5.25.62.gtf", asRangedData=FALSE)
        targets <- import.gff2(gtfFiles, asRangedData=FALSE)
        namesGTF <- paste("chr", seqlevels(targets), sep="")
        seqlevels(targets) <- namesGTF
        targets <- targets[seqnames(targets)!="chrdmel_mitochondrion_genome"]
        #gRange <- targets
    
        levels <- levels(elementMetadata(targets)[,"source"])
        snoRNA <- targets[elementMetadata(targets)[,"source"]=="snoRNA"]
        snRNA <- targets[elementMetadata(targets)[,"source"]=="snRNA"]
        tRNA <- targets[elementMetadata(targets)[,"source"]=="tRNA"]
        miRNA <- targets[elementMetadata(targets)[,"source"]=="miRNA"]
        rRNA <- targets[elementMetadata(targets)[,"source"]=="rRNA"]
        ##protein <- targets[elementMetadata(targets)[,"source"]=="protein_coding"]
        gRange <- c(snoRNA, snRNA, tRNA, miRNA, rRNA)
        rm(targets)
    }
    
    bamFiles <- unlist(qProject@alignments$genome)
    counts <- lapply(bamFiles, .countAlignments, gRange)
    counts <- as(counts,"DataFrame") 
    #counts <- as.data.frame(counts)
    colnames(counts) <- as.character(qProject@samples$name)
    rownames(counts) <- NULL
    values(gRange) <- IRanges::cbind(values(gRange), counts)
    #values(gRange) <- cbind.data.frame(as.data.frame(val), counts)
    return(gRange)

}
