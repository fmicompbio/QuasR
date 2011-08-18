qCount <- function(qProject, gRange, stranded=FALSE, ...)
{
    if(missing(gRange)){
        isGTFFormat <- .fileExtension(qProject@annotations$filepath) %in% c("gtf")
        gtfFiles <- qProject@annotations[isGTFFormat,]$filepath
        if(length(gtfFiles) != 1)
            stop("There is more/less than one 'gtf' file in the annotation")
        require(rtracklayer)
        #targets <- import.gff2("data/Drosophila_melanogaster.BDGP5.25.62.gtf", asRangedData=FALSE)
        #targets <- import.gff2(file(as.character(gtfFiles), "r"), asRangedData=FALSE)
        gRange <- import.gff2(as.character(gtfFiles), asRangedData=FALSE)
        ## HACK fix of wrong sequence name
        namesGTF <- paste("chr", seqlevels(gRange), sep="")
        seqlevels(gRange) <- namesGTF
        gRange <- gRange[seqnames(gRange)!="chrdmel_mitochondrion_genome"]
        ## subset GRAnges object with features from the annotation file
        levels <- levels(elementMetadata(gRange)[,"source"])
        queryTarget <- unlist(strsplit(as.character(project@annotations[isGTFFormat,]$feature), ","))
        #if(!queryTarget %in% levels)
        #    stop("The source column of the 'gtf' files contains '", levels, "' but you query for '", queryTarget, "'.")
        gRange <- targets[elementMetadata(gRange)[,"source"] == queryTarget]
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