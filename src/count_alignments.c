/*!
  @header

  TODO.
 
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "quasr.h"

/*! @typedef
  @abstract Structure for the todo.
  @field todo   todo
  @field todo todo
 */
typedef struct {
    double sum;
    int start;
    int end;
    const char * strand;
    int stranded;
    const char * overlap;
} fetch_param;

/*! @function
  @abstract  callback for bam_fetch()
  @param  todo  todo
  @return    todo
 */
static int fetch_func(const bam1_t *hit, void *data){
    fetch_param *fparam = (fetch_param*)data;
    int check_strand = 1;
    int check_within = 1;
    // check if whole alignment is in region
    if(strcmp(fparam->overlap, "within") == 0){
	if(!(fparam->start <= hit->core.pos) 
	   || !(bam_calend(&hit->core, bam1_cigar(hit)) <= fparam->end)){
	    check_within = 0;
//	    Rprintf("%swithin-", bam1_qname(hit));
	}
    }
    // check strand
    if(fparam->stranded && strcmp(fparam->strand, "*") != 0){
	if((strcmp(fparam->strand,"-") == 0) != ((hit->core.flag & BAM_FREVERSE) == 16)){
	    check_strand = 0;
//	    Rprintf("%sstrand%i-", fparam->strand, (hit->core.flag & BAM_FREVERSE));
	}
    }
    // sum up weight of alignment
    if(check_strand && check_within){
	int32_t w = bam_aux2i(bam_aux_get(hit,"IH"));
	//TODO eventual check if w is really not 0 because of division
	fparam->sum += 1.0/w;
    }
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
SEXP count_alignments(SEXP bam_in, SEXP idx_in, SEXP regions, SEXP stranded, SEXP overlap)
{
    // check bam and index file and open it
    if(!IS_CHARACTER(bam_in))
        Rf_error("'bam_in' must be character");

    if(!IS_CHARACTER(idx_in))
        Rf_error("'idx_in' must be character");

    if(LENGTH(bam_in) != LENGTH(idx_in))
	Rf_error("The 'bam' file list and 'bam index' file list must be of equal length.");

    //if(LENGTH(bam_in) != 1)
//	Rf_error("The 'bam' file list must be of length 1.");

    int num_files = LENGTH(bam_in);
    samfile_t **fin;
    bam_index_t **idx;

    fin = (samfile_t**)calloc(num_files, sizeof(samfile_t*));
    idx = (bam_index_t**)calloc(num_files, sizeof(bam_index_t*));

    for(int j = 0; j < num_files; j++){
	fin[j] = _bam_tryopen(translateChar(STRING_ELT(bam_in, j)), "rb", NULL);
	idx[j] = bam_index_load(translateChar(STRING_ELT(idx_in, j)));
	//TODO error handling. f_in leaks if this fails
    }

    // check parameter stranded and overlap type
    if(!IS_LOGICAL(stranded) || LENGTH(stranded) != 1)
	Rf_error("'stranded must be of type locical(1)");

    if(!IS_CHARACTER(overlap) || LENGTH(overlap) != 1)
	Rf_error("'type must be of type character(1)");
    
    // check parameter region and get direct pointer to the elements
    SEXP tid = getListElement(regions, "tid");
    SEXP start = getListElement(regions, "start");
    SEXP end = getListElement(regions, "end");
    SEXP strand = getListElement(regions, "strand");
    
    if(!IS_INTEGER(tid))
        Rf_error("Column 'tid' must be of type integer");

    if(!IS_INTEGER(start))
        Rf_error("Column 'start' must be of type integer");

    if(!IS_INTEGER(end))
        Rf_error("Column 'end' must be of type integer");

    if(!IS_CHARACTER(strand))
        Rf_error("Column 'strand' must be of type character");

    int num_regions = LENGTH(tid);
    
    if(num_regions != LENGTH(start) || num_regions != LENGTH(end) || num_regions != LENGTH(strand))
        Rf_error("The columns 'tid', 'start', 'end', 'stand' must have equal length.");

    // initialise and run bam_fetch
    SEXP count;
    PROTECT(count = NEW_NUMERIC(num_regions));
    fetch_param fparam;
    fparam.stranded = LOGICAL(stranded)[0];
    fparam.overlap = translateChar(STRING_ELT(overlap, 0));

    for(int i = 0; i < num_regions; i++){
	fparam.sum = 0.0;
	fparam.start = INTEGER(start)[i];
	fparam.end = INTEGER(end)[i];
	fparam.strand =	translateChar(STRING_ELT(strand, i));
	
	for(int j = 0; j < num_files; j++){
	    bam_fetch(fin[j]->x.bam, idx[j], INTEGER(tid)[i], INTEGER(start)[i], INTEGER(end)[i], &fparam, fetch_func);
	    //REAL(count)[i] = _count_alignments(fin, idx, INTEGER(tid)[i], INTEGER(start)[i], INTEGER(end)[i]);
	}
	REAL(count)[i] = fparam.sum;
    }

    for(int j = 0; j < num_files; j++){
	samclose(fin[j]);
	bam_index_destroy(idx[j]);;
    }    
    //samclose(fin);
    //bam_index_destroy(idx);
    UNPROTECT(1);

    return count;
}
