<img src="vignettes/QuasR-logo.png" align="left" alt="" />
<br>

# Quantify and annotate short reads in R 

<br>

## Overview

`QuasR` aims to cover the whole analysis workflow of typical ultra-high throughput
sequencing experiments, starting from the raw sequence reads, over pre-processing and
alignment, up to quantification. A single **R** script can contain all steps of a complete
analysis, making it simple to document, reproduce or share the workflow containing all
relevant details.

The `QuasR` package integrates the functionality of several **R** packages
(such as `IRanges` and `Rsamtools`) and external software (e.g. `bowtie`,
through the `Rbowtie` package).

The package was originally created for in house use at the Friedrich Miescher Institute.
Current contributors include:

- [Michael Stadler](https://github.com/mbstadler)
- Dimos Gaidatzis
- [Charlotte Soneson](https://github.com/csoneson)
- Anita Lerch
- Florian Hahne

## Functionality

The current `QuasR` release supports the analysis of single read and
paired-end ChIP-seq (chromatin immuno-precipitation combined with sequencing), RNA-seq
(gene expression profiling by sequencing of RNA) and Bis-seq (measurement of DNA
methylation by sequencing of bisulfite-converted genomic DNA) experiments. Allele
specific alignment and quantification is also supported in all of these experiment
types. `QuasR` has been successfully used with data from Illumina, 454 Life Technologies
and SOLiD sequencers, the latter by using bam files created externally `QuasR`.

## Reference
`QuasR` has been described in:  

"QuasR: Quantify and Annotate Short Reads in R"  
Gaidatzis, D. and Lerch, A. and Hahne, F. and Stadler, M.B.  
*Bioinformatics* **2015**, 31(7):1130-1132.  
[PubMed: 25417205](https://www.ncbi.nlm.nih.gov/pubmed/25417205), [doi: 10.1093/bioinformatics/btu781](https://doi.org/10.1093/bioinformatics/btu781)

## Download from Bioconductor
[QuasR download page](https://bioconductor.org/packages/QuasR/)

## Software status

| Platforms          |  OS              | R CMD check      | Coverage         | 
|:------------------:|:----------------:|:----------------:|:----------------:|
| Travis CI          | Linux            | [![Travis CI build status](https://travis-ci.com/fmicompbio/QuasR.svg?branch=master)](https://travis-ci.com/fmicompbio/QuasR) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
| Github-Actions     | Multiple         |  [![R build status](https://github.com/fmicompbio/QuasR/workflows/R-CMD-check/badge.svg)](https://github.com/fmicompbio/QuasR/actions) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
| Bioc ([_devel_](http://bioconductor.org/packages/devel/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/QuasR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/QuasR) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
| Bioc ([_release_](http://bioconductor.org/packages/release/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/QuasR.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/QuasR) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
