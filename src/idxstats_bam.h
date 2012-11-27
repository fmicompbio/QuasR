#include <ctype.h>
#include <assert.h>
#include <Rdefines.h>
#include "samtools/bam.h"
#include "samtools/khash.h"
#include "samtools/ksort.h"
#include "samtools/bam_endian.h"
#ifdef _USE_KNETFILE
#include "samtools/knetfile.h"
#endif

SEXP idxstats_bam(SEXP inBam);

