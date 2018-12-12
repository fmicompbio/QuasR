// prevent remapping (e.g. Rf_error to error), which causes conflicts under windows
#define STRICT_R_HEADERS 1

// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "samtools/sam.h"

SEXP cat_bam(SEXP inbam, SEXP outbam);
