/*!
  @header

  TODO.
 
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "quasr.h"

/*! @function
  @abstract  callback for bam_fetch()
  @param  todo  todo
  @return    todo
 */
static int fetch_func(const bam1_t *hit, void *data){
    //TODO check if whole alignment is in region hit->core.pos
    double *sum = (double*)data;
    int32_t c = bam_aux2i(bam_aux_get(hit,"IH"));
    *sum += 1.0/c; //TODO check if c is really 0 because of division
    return 0;  
} 

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
double _count_alignments(samfile_t *fin, bam_index_t *idx, int tid, int start, int end)
{
    double sum = 0.0;
    int res = bam_fetch(fin->x.bam, idx, tid, start, end, &sum, fetch_func);
    return sum;
}

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
SEXP count_alignments(SEXP bam_in, SEXP idx_in, SEXP regions, SEXP type)
{
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
        Rf_error("'bam_in' must be character(1)");

    if (!IS_CHARACTER(idx_in) || 1 != LENGTH(idx_in))
        Rf_error("'idx_in' must be character(1)");
    samfile_t *fin = _bam_tryopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);    
    //TODO error handling. f_in leaks if this fails
    bam_index_t *idx = bam_index_load(translateChar(STRING_ELT(idx_in, 0))); // f_in leaks if this fails 

    SEXP tid = getListElement(regions, "tid");
    SEXP start = getListElement(regions, "start");
    SEXP end = getListElement(regions, "end");
    
    if (!IS_INTEGER(tid))
        Rf_error("Column 'tid' must be of type integer");

    if (!IS_INTEGER(start))
        Rf_error("Column 'start' must be of type integer");

    if (!IS_INTEGER(end))
        Rf_error("Column 'end' must be of type integer");

    if (!LENGTH(tid) == LENGTH(start) && LENGTH(tid) == LENGTH(end))
        Rf_error("Columns 'tid', 'start', 'end' must be of equal length");

    int num_regions = LENGTH(tid);

    SEXP count;
    PROTECT(count = NEW_NUMERIC(num_regions));
    for(int i = 0; i < num_regions; i++){
        REAL(count)[i] = _count_alignments(fin, idx, INTEGER(tid)[i], INTEGER(start)[i], INTEGER(end)[i]);
    }
    
    samclose(fin);
    bam_index_destroy(idx);
    UNPROTECT(1);

    return count;
}
