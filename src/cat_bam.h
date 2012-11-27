// prevent remapping (e.g. Rf_error to error), which causes conflicts under windows
#define STRICT_R_HEADERS 1

#include <Rdefines.h>
#include "samtools/sam.h"

SEXP cat_bam(SEXP inbam, SEXP outbam);
