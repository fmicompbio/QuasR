// prevent remapping of "Rf_length" to "length" in Rdefines.h, which clashes with fstream::length
#define R_NO_REMAP
// prevent remapping e.g. Ralloc, which causes conflicts under windows
#define STRICT_R_HEADERS 1

#include <Rdefines.h>
#include <R.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <queue>
#include <map>
#include <algorithm>
#include "samtools/sam.h"


int _merge_reorder_sam(const char** fnin, int nin, const char* fnout, int mode, int maxhits);

#ifdef __cplusplus
extern "C" {
#endif
    SEXP merge_reorder_sam(SEXP infiles, SEXP outfile, SEXP mode, SEXP maxhits);
#ifdef __cplusplus
}
#endif
