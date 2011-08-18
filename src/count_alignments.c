#include "quasr.h"


// callback for bam_fetch()
static int fetch_func(const bam1_t *hit, void *data){
    //TODO check if whole alignment is in region hit->core.pos
    double *sum = (double*)data;
    int32_t c = bam_aux2i(bam_aux_get(hit,"IH"));
    *sum += 1.0/c; //TODO check if c is really 0 because of division
    return 0;  
} 

double _count_alignments(samfile_t *fin, bam_index_t *idx, int tid, int start, int end)
{
    double sum = 0.0;
    int res = bam_fetch(fin->x.bam, idx, tid, start, end, &sum, fetch_func);
    return sum;
}

SEXP count_alignments(SEXP bam_in, SEXP idx_in, SEXP tid, SEXP start, SEXP end)
{
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
        Rf_error("'bam_in' must be character(1)");

    if (!IS_CHARACTER(idx_in) || 1 != LENGTH(idx_in))
        Rf_error("'idx_in' must be character(1)");
    
    if (!IS_INTEGER(tid))
        Rf_error("'tid' must be of type integer");

    if (!IS_INTEGER(start))
        Rf_error("'start' must be of type integer");

    if (!IS_INTEGER(end))
        Rf_error("'end' must be of type integer");

    if (!LENGTH(tid) == LENGTH(start) && LENGTH(tid) == LENGTH(end))
        Rf_error("'tid', 'start', 'end' must be of equal length");

    if (1 != LENGTH(tid))
        Rf_error("Length of 'tid', 'start' or 'end' must be 1");

    samfile_t *fin = _bam_tryopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);    
    //TODO error handling. f_in leaks if this fails
    bam_index_t *idx = bam_index_load(translateChar(STRING_ELT(idx_in, 0))); // f_in leaks if this fails 

    //int tid = 1;
    //int start = 0;
    //int end = fin->header->target_len[tid]-1;
    double count = _count_alignments(fin, idx, asInteger(tid), asInteger(start), asInteger(end));
    
    samclose(fin);
    bam_index_destroy(idx);
    
    return ScalarReal(count);
}