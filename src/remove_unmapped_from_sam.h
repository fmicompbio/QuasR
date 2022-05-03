// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include <samtools-1.7-compat.h>
#include <stdio.h>

SEXP remove_unmapped_from_sam_and_convert_to_bam(SEXP inSam, SEXP outBam);

