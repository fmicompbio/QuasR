// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "utilities.h"
#include "samtools/sam.h"

SEXP split_sam_chr(SEXP samFile, SEXP outDir);
