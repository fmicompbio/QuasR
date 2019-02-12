#include <ctype.h>
#include <assert.h>
// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "sam.h"
#include "htslib/khash.h"
#include "htslib/ksort.h"
#include "bam_endian.h"
#ifdef _USE_KNETFILE
#include "htslib/knetfile.h"
#endif

SEXP idxstats_bam(SEXP inBam);

