#include <ctype.h>
#include <assert.h>
// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include <samtools-1.7-compat.h>
#include "htslib/khash.h"
#include "htslib/ksort.h"
#ifdef _USE_KNETFILE
#include "htslib/knetfile.h"
#endif

SEXP idxstats_bam(SEXP inBam);

