// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "samtools/sam.h"
#include "utilities.h"

SEXP count_alignments_non_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                      SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

SEXP count_alignments_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                      SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax);

