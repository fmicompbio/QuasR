#include <Rdefines.h>
#include "samtools/sam.h"

SEXP count_alignments_non_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden);

SEXP count_alignments_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand, 
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden);

