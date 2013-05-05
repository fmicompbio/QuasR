#include <Rdefines.h>
#include "samtools/sam.h"
#include "utilities.h"


SEXP count_alignments_subregions(SEXP bam_in, SEXP idx_in, SEXP tid, SEXP min_regions,
                                 SEXP max_regions, SEXP regions, SEXP shift, SEXP broaden,
                                 SEXP overlap_type, SEXP includeSpliced);
