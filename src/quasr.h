#ifndef _QUASR_H_
#define _QUASR_H_

//#include "samtools/bam.h"
#include "samtools/sam.h"
//#include <R.h>
//#include <Rinternals.h>
#include <Rdefines.h>

//void _write_buffered_alignment(samfile_t *fout, bam1_t **buf, uint16_t n);
//void _count_alignments(const char *fn_in, const char *fn_out, int max_hits, size_t max_mem);

//utilities.c
samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux);
SEXP seqname(SEXP filename);
SEXP getListElement(SEXP list, const char *str);

//weight_alignments.c
SEXP get_allowed_max_hits();
SEXP weight_alignments(SEXP bam_in, SEXP bam_out, SEXP max_hits);

//count_alignments.c
SEXP count_alignments(SEXP bam_in, SEXP idx_in, SEXP regions, SEXP stranded, SEXP overlap);
#endif /* _QUASR_H_ */
