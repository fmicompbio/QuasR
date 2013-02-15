// prevent remapping of "Rf_length" to "length" in Rdefines.h, which clashes with fstream::length
#define R_NO_REMAP
// prevent remapping e.g. Ralloc, which causes conflicts under windows
#define STRICT_R_HEADERS 1

#include <Rdefines.h>
#include "samtools/sam.h"

SEXP count_junctions(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP allelic);
