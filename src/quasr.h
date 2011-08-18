#ifndef _QUASR_H_
#define _QUASR_H_

//#include "samtools/bam.h"
#include "samtools/sam.h"
//#include <R.h>
//#include <Rinternals.h>
#include <Rdefines.h>


//utilities.c
samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux);
SEXP seqname(SEXP filename);

//weight_alignments.c
SEXP get_allowed_max_hits();
SEXP weight_alignments(SEXP bam_in, SEXP bam_out, SEXP max_hits);
//void _write_buffered_alignment(samfile_t *fout, bam1_t **buf, uint16_t n);
//void _count_alignments(const char *fn_in, const char *fn_out, int max_hits, size_t max_mem);

#endif /* _QUASR_H_ */
