createProjectSingle <- function(allelic=FALSE) {
  requireNamespace("Rsamtools", quietly = TRUE)
  
  # create bam files
  samfile_plus <- tempfile(fileext = ".sam", tmpdir = "extdata")
  samfile_minus <- tempfile(fileext = ".sam", tmpdir = "extdata")
  cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chrV\tLN:800\n", file = samfile_plus)
  cat("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chrV\tLN:800\n", file = samfile_minus)
  pos <- seq.int(512)
  cat(paste("seq1\t0\tchrV", pos, pos %% 256, "10M\t*\t0\t0\t*",
            if (allelic) rep(c("*\tXV:A:R","*\tXV:A:A","*\tXV:A:U"), each = length(pos)) else "*", 
            sep = "\t", collapse = "\n"), file = samfile_plus, append = TRUE)
  cat(paste("seq1\t16\tchrV", pos, pos %% 256, "10M\t*\t0\t0\t*",
            if (allelic) rep(c("*\tXV:A:R","*\tXV:A:A","*\tXV:A:U"), each = length(pos)) else "*",
            sep = "\t", collapse = "\n"), file = samfile_minus, append = TRUE)
  bamfile_plus <- Rsamtools::asBam(samfile_plus, indexDestination = TRUE)
  bamfile_minus <- Rsamtools::asBam(samfile_minus, indexDestination = TRUE)
  
  # create genome
  genome <- tempfile(fileext = ".fa", tmpdir = "extdata")
  cat(">chrV\n", paste(rep("G", 600), collapse = ""), "\n", file = genome)

  # create SNP file
  if (allelic) {
    snp <- tempfile(fileext = ".txt", tmpdir = "extdata")
    cat(">chrV\t10\tC\tG\n", file = snp)
  }
  
  # create sample file
  samplefile <- tempfile(fileext = ".txt", tmpdir = "extdata")
  write.table(data.frame(FileName = basename(c(bamfile_plus, bamfile_minus)),
                         SampleName = c("Sample", "Sample"), stringsAsFactors = FALSE),
              sep = "\t", quote = FALSE, row.names = FALSE, file = samplefile)
  
  # qAlign
  td <- tempdir()
  if (!allelic) {
    proj <- qAlign(samplefile, genome, paired = "no", alignmentsDir = td)
  } else {
    proj <- qAlign(samplefile, genome, snpFile = snp, paired = "no", alignmentsDir = td)        
  }
  return(proj)
}

createProjectPaired <- function() {
  requireNamespace("Rsamtools", quietly = TRUE)
  
  # create samfile
  sfile <- tempfile(fileext = ".sam", tmpdir = "extdata")
  cat("@HD\tVN:1.0\tSO:unsorted\n", file = sfile, append = FALSE)
  cat("@SQ\tSN:chrV\tLN:99\n", file = sfile, append = TRUE)
  # fr R1->left R2->right
  cat("seq1\t99\tchrV\t4\t255\t5M\tchrV\t12\t13\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq1\t147\tchrV\t12\t255\t5M\tchrV\t4\t-13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq2\t99\tchrV\t4\t255\t5M\tchrV\t13\t14\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq2\t147\tchrV\t13\t255\t5M\tchrV\t4\t-14\t*\t*\n", file = sfile, append = TRUE)
  # fr R2->left R1->right
  cat("seq3\t163\tchrV\t24\t255\t5M\tchrV\t32\t13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq3\t83\tchrV\t32\t255\t5M\tchrV\t24\t-13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq4\t163\tchrV\t24\t255\t5M\tchrV\t33\t14\t*\t*\n", file = sfile, append = TRUE)
  cat("seq4\t83\tchrV\t33\t255\t5M\tchrV\t24\t-14\t*\t*\n", file = sfile, append = TRUE)
  # ff R1->left R2->right
  cat("seq5\t67\tchrV\t44\t255\t5M\tchrV\t52\t13\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq5\t131\tchrV\t52\t255\t5M\tchrV\t44\t-13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq6\t67\tchrV\t44\t255\t5M\tchrV\t53\t14\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq6\t131\tchrV\t53\t255\t5M\tchrV\t44\t-14\t*\t*\n", file = sfile, append = TRUE)
  # rr R2->left R1->right
  cat("seq7\t115\tchrV\t64\t255\t5M\tchrV\t72\t13\t*\t*\n",  file = sfile, append = TRUE)
  cat("seq7\t179\tchrV\t72\t255\t5M\tchrV\t64\t-13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq8\t115\tchrV\t64\t255\t5M\tchrV\t73\t14\t*\t*\n",  file = sfile, append = TRUE)
  cat("seq8\t179\tchrV\t73\t255\t5M\tchrV\t64\t-14\t*\t*\n", file = sfile, append = TRUE)
  # rf R1->left R2->right
  cat("seq9\t83\tchrV\t84\t255\t5M\tchrV\t92\t13\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq9\t163\tchrV\t92\t255\t5M\tchrV\t84\t-13\t*\t*\n", file = sfile, append = TRUE)
  cat("seq0\t83\tchrV\t84\t255\t5M\tchrV\t93\t14\t*\t*\n",   file = sfile, append = TRUE)
  cat("seq0\t163\tchrV\t93\t255\t5M\tchrV\t84\t-14\t*\t*\n", file = sfile, append = TRUE)
  
  # create bamfile
  bfile <- Rsamtools::asBam(sfile)
  
  # create sample file
  samplefile <- tempfile(fileext = ".txt", tmpdir = "extdata")
  write.table(data.frame(FileName = basename(bfile),
                         SampleName = "Normal", stringsAsFactors = FALSE),
              sep = "\t", quote = FALSE, row.names = FALSE, file = samplefile)
  
  # create genome
  genome <- tempfile(fileext = ".fa", tmpdir = "extdata")
  cat(">chrV\n", file = genome, append = FALSE)
  cat(paste(rep("G",99), collapse = ""), file = genome, append = TRUE)
  
  # qAlign
  td <- tempdir()
  project <- qAlign(samplefile, genome, paired = "fr", alignmentsDir = td)
  
  return(project)
}

createTiles <- function() {
  #        101  151  201  251  301  351  401  451  501  551  601  651  701     X = 10x 3reads(=R,U,A) 
  #   H1 :   XXXXXXXXXX          XXXXXXXXXX          XXXXXXXXXX            --> 30X => 900
  #   H2 :        XXXXXXXXXX          XXXXXXXXX           XXXXXXXXXX       --> 30X => 900
  #   H3 :             XXXXXXXXXX          XXXXXXXXXX          XXXXXXXXXX  --> 30X => 900
  #   H4 :                  XXXXXXXXXX          XXXXXXXXXX          XXXXX  --> 24X => 750
  
  requireNamespace("GenomicRanges", quietly = TRUE)
  
  width = 100L
  gap = 100L
  times = 600L / (width + gap)
  h1 <- GenomicRanges::GRanges("chrV", IRanges::successiveIRanges(rep(width,times),gap, from = 101))
  names(h1) <- rep("H1", length(h1))
  h2 <- GenomicRanges::GRanges("chrV", IRanges::successiveIRanges(rep(width,times),gap, from = 151))
  names(h2) <- rep("H2", length(h2))
  h3 <- GenomicRanges::GRanges("chrV", IRanges::successiveIRanges(rep(width,times),gap, from = 201))
  names(h3) <- rep("H3", length(h3))
  h4 <- GenomicRanges::GRanges("chrV", IRanges::successiveIRanges(rep(width,times),gap, from = 251))
  names(h4) <- rep("H4", length(h4))
  end(h4[3]) <- 700
  tiles <- sort(c(h1, h2, h3, h4))    
  return(tiles)
}

createGtfGr <- function() {
  requireNamespace("rtracklayer")
  requireNamespace("GenomicRanges")
  gtfFile <- system.file("extdata", "hg19sub_annotation.gtf", package = "QuasR")
  gr <- rtracklayer::import(gtfFile, format = "gtf", feature.type = "exon")
  names(gr) <- mcols(gr)$gene_id
  gr <- gr[order(names(gr))]
  return(gr)
}

createTxDb <- function() {
  requireNamespace("GenomicFeatures")
  gtfGr <- createGtfGr()
  txdb <- GenomicFeatures::makeTxDbFromGRanges(gtfGr)
  return(txdb)
}

# createBSgenome <- function() {
#   requireNamespace("BSgenome")
#   
#   # prepare sequences and seed files
#   seqdir      <- tempfile(pattern = "BSgenome.seqs",    tmpdir = "extdata")
#   bsgenomedir <- tempfile(pattern = "BSgenome.hg19sub", tmpdir = "extdata")
#   rlibdir     <- tempfile(pattern = "Rlib",             tmpdir = "extdata")
# 
#   pkgseed <- new(Class = "BSgenomeDataPkgSeed",
#                  Package = "BSgenome.HSapiens.QuasR.hg19sub",
#                  Title = "Nucleotide sequences of subregions from hg19",
#                  Description = "Nucleotide sequences of three short chromosomal subregions from hg19",
#                  Version = "0.1.0",
#                  Author = "Michael Stadler",
#                  Maintainer = "Michael Stadler <michael.stadler@fmi.ch>",
#                  License = "GPL-2",
#                  organism = "Homo Sapiens",
#                  common_name = "Human",
#                  provider = "QuasR",
#                  provider_version = "hg19sub",
#                  release_date = "Feb. 2009",
#                  release_name = "Genome Reference Consortium GRCh37",
#                  organism_biocview = "Homo Sapiens",
#                  BSgenomeObjname = "hg19sub",
#                  seqnames = "paste0('chr',1:3)")
#   
#   chrs <- Biostrings::readDNAStringSet(file.path("extdata", "hg19sub.fa"))
#   dir.create(seqdir)
#   for (i in seq_along(chrs))
#     Biostrings::writeXStringSet(chrs[i], format = "fasta",
#                                 filepath = file.path(seqdir, paste0(names(chrs[i]), ".fa")))
#   
#   # build and install the package
#   dir.create(bsgenomedir)
#   uniqueLetters <- Biostrings::uniqueLetters # without this, BSgenome::forgeBSgenomeDataPkg fails
#   BSgenome::forgeBSgenomeDataPkg(x = pkgseed, seqs_srcdir = seqdir, destdir = bsgenomedir)
#   dir.create(rlibdir)
#   utils::install.packages(pkgs = file.path(bsgenomedir, "BSgenome.hg19sub"),
#                           lib = rlibdir, repos = NULL, type = "source")
#   
#   return(rlibdir)
# }
