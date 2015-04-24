// prevent remapping (e.g. Rf_error to error), which causes conflicts under windows
#define STRICT_R_HEADERS 1

#include <stdlib.h>
#include <math.h>
#include "zlib.h"
#include <Rdefines.h>
#include "utilities.h"


SEXP bamfile_to_wig(SEXP _bam_in, SEXP _wig_out, SEXP _paired, SEXP _binsize, SEXP _shift,
                    SEXP _strand, SEXP _norm_factor, SEXP _tracknames, SEXP _log2p1,
                    SEXP _colors, SEXP _compress, SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);
