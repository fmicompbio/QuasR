#include <Rdefines.h>
#include "samtools/sam.h"
#include "utilities.h"

SEXP count_alignments_non_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                      SEXP mapqMin, SEXP mapqMax);

SEXP count_alignments_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                      SEXP mapqMin, SEXP mapqMax);

