#include <Rdefines.h>
#include "samtools/sam.h"
#include "utilities.h"

SEXP profile_alignments_non_allelic(SEXP bamfile, SEXP targetprofile, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                    SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                    SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                    SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

SEXP profile_alignments_allelic(SEXP bamfile, SEXP targetprofile, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

