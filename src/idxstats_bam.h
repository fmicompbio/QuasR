#include <ctype.h>
#include <assert.h>
// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include "samtools/bam.h"
#include "samtools/khash.h"
#include "samtools/ksort.h"
#include "samtools/bam_endian.h"
#ifdef _USE_KNETFILE
#include "samtools/knetfile.h"
#endif

SEXP idxstats_bam(SEXP inBam);

