setGeneric("qCount", function(qproject, query, stranded=FALSE, collapseSamples=TRUE) standardGeneric("qCount"))

setMethod("qCount",
          signature(qproject="qProject", query="missing"),
          function(qproject, query, stranded, collapseSamples)
      {
          qCount(qproject, "summary", stranded, collapseSamples)
      })

setMethod("qCount",
          signature(qproject="qProject", query="character"),
          function(qproject, query, stranded, collapseSamples)
      {

          if(query != "summary" && !query %in% qproject@annotations$feature)
              error("Wrong query statment")

          .progressReport(sprintf("Load annotation for %s", query), phase=-1)

          if(query == "summary")
              query <- unique(qproject@annotations$feature)

          annotations <- qproject@annotations$feature %in% query

          #isGTFFormat <- .fileExtension(qproject@annotations$filepath) %in% c("gtf")
          #gtfFiles <- qproject@annotations[isGTFFormat,]$filepath
          gtfFiles <- annotations[ annotations$filetype %in% "gtf", "filepath" ]
          gtfFiles <- unique(gtfFiles)
          if(length(gtfFiles) != 1)
              stop("There is more or less than one 'gtf' file in the annotation")
          require(rtracklayer)
          gRange <- import.gff2(as.character(gtfFiles), asRangedData=FALSE,
                                colnames=c("strand", "type", "source", "gene_id", "transcript_id", "exon_number"))

          gRange <- gRange[values(gRange)[, "type"] == "exon"]
          seqlevels(gRange) <- .mapSeqnames(names(getGenomeInformation(qproject)), seqlevels(gRange))
          ## subset GRAnges object with features from the annotation file
          levels <- levels(elementMetadata(gRange)[,"source"])
          #queryTarget <- unlist(strsplit(as.character(qproject@annotations[isGTFFormat,]$feature), ","))
          #if(!queryTarget %in% levels)
          #    stop("The source column of the 'gtf' files contains '", levels, "' but you query for '", queryTarget, "'.")
          gRange <- gRange[ elementMetadata(gRange)[,"source"] %in% query ]
          .progressReport("Successfully loaded the annotation.", phase=1)
          qCount(qproject, gRange, stranded, collapse)
      })

setMethod("qCount",
          signature(qproject="qProject", query="GRanges"),
          function(qproject, query, stranded, collapseSamples)
      {
          .progressReport("Starting count alignments", phase=-1)
          bamFiles <- unlist(qproject@alignments$genome)

          if(collapseSamples == TRUE){
              counts <- lapply(split(bamFiles, qproject@samples$name), .countAlignments, query)
              #names(counts) <- as.character(qproject@samples$name)
          }else{
              counts <- lapply(bamFiles, .countAlignments, query)
              names(counts) <- basename(qproject@samples$filepath)
          }
    
          #counts <- as(counts,"DataFrame")
          counts <- as.data.frame(counts)
          rownames(counts) <- names(query)
          #values(query) <- IRanges::cbind(values(query), counts)
          #values(gRange) <- cbind.data.frame(as.data.frame(val), counts)
          .progressReport("Successfully terminated the quasr counting.", phase=1)
          return(counts)
      })
