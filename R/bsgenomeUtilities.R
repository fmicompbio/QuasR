downloadBSgenome <- function (graphics=getOption("menu.graphics"), overwrite=FALSE, ...) 
{
    if(!interactive()) 
        stop("cannot choose a BSgenome non-interactively")
    require(BSgenome, quietly=TRUE)
    genome <- character()
    m <- available.genomes()
    res <- menu(m, graphics, "Available BSgnomes")
    if(res > 0L) {
        genome <- m[res]
        if(genome %in% installed.genomes() && !overwrite){
          return(genome)
        }
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(genome, ...)
    }
    return(genome)
}

choseBSgenome <- function (graphics=getOption("menu.graphics")) 
{
    if(!interactive()) 
        stop("cannot choose a BSgenome non-interactively")
    require(BSgenome, quietly=TRUE)
    genome <- character()
    m <- installed.genomes()
    m <- c(m, "Download an other BSgenome")
    res <- menu(m, graphics, "Installed BSgnomes")
    if(res == length(m)){
      genome <- downloadBSgenome()
    } else {
      if(res > 0L) {
        genome <- m[res]
      }
    }
    return(genome)
}
