// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "htslib/sam.h"
#include "utilities.h"

SEXP profile_alignments_non_allelic(SEXP bamfile, SEXP targetprofile, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                    SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                    SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                    SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

SEXP profile_alignments_allelic(SEXP bamfile, SEXP targetprofile, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

