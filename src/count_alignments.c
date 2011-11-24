/*!
  @header

  TODO.
 
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "quasr.h"
#include <stdlib.h>


/*! @typedef
  @abstract Structure for the todo.
  @field todo   todo
  @field todo todo
 */
typedef struct {
    double sum;
    int start;
    int end;
    int shift;
    const char * strand;
    int minoverlap;
    int maxhit;
} fetch_param;

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
int32_t get_inverse_weight(const bam1_t *b){
    uint8_t *w = bam_aux_get(b, "IH");
    //Rprintf("Ptr-%p ", w);
    if(w != 0)
	return bam_aux2i(w);
    else
	return 1;
    return 0;
}

/*! @function
  @abstract  callback for bam_fetch() to sum up weight of alignment
  @param  todo  todo
  @return    todo
 */
static int _fetch_any(const bam1_t *hit, void *data){
    fetch_param *fparam = (fetch_param*)data;
    //int32_t w = bam_aux2i(bam_aux_get(hit,"IH"));
    //fparam->sum += 1.0/w;
    fparam->sum += 1.0/get_inverse_weight(hit);
    return 0;
}

static int _fetch_any_shift(const bam1_t *hit, void *data){
    fetch_param *fparam = (fetch_param*)data;
    if((hit->core.flag & BAM_FREVERSE) != 16){
	//plus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "+") == 0))
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) - (fparam->start - fparam->shift) >= fparam->minoverlap
	   && (fparam->end - fparam->shift) - hit->core.pos >= fparam->minoverlap)
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit) //TODO length(hit) if minoverlap filterout also to short reads
		fparam->sum += 1.0/w;
	}
    }else{
	//minus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "-") == 0))
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) - (fparam->start + fparam->shift) >= fparam->minoverlap
	   && (fparam->end + fparam->shift) - hit->core.pos >= fparam->minoverlap)
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }
    return 0;
}

// no overlap check needed because it must be completly within
static int _fetch_within_shift(const bam1_t *hit, void *data){
//Rprintf("within");
    fetch_param *fparam = (fetch_param*)data;
    if((hit->core.flag & BAM_FREVERSE) != 16)
    {
	//plus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "+") == 0))
	   && (fparam->start - fparam->shift) <= hit->core.pos
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) <= (fparam->end - fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }else{
	//minus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "-") == 0))
	   && (fparam->start + fparam->shift) <= hit->core.pos 
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) <= (fparam->end + fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }
    return 0;
}

static int _fetch_startwithin_shift(const bam1_t *hit, void *data){
//Rprintf("startwithin");
    fetch_param *fparam = (fetch_param*)data;
    if((hit->core.flag & BAM_FREVERSE) != 16)
    {
	//plus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "+") == 0))
	   && (fparam->start - fparam->shift) <= hit->core.pos
	   && hit->core.pos <= (fparam->end - fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }else{
	//minus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "-") == 0))
	   && (fparam->start + fparam->shift) <= (int)bam_calend(&hit->core, bam1_cigar(hit)) 
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) <= (fparam->end + fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }
    return 0;
}

static int _fetch_endwithin_shift(const bam1_t *hit, void *data){
//Rprintf("endwithin");
    fetch_param *fparam = (fetch_param*)data;
    if((hit->core.flag & BAM_FREVERSE) != 16)
    {
	//plus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "+") == 0))
	   && (fparam->start - fparam->shift) <= (int)bam_calend(&hit->core, bam1_cigar(hit))
	   && (int)bam_calend(&hit->core, bam1_cigar(hit)) <= (fparam->end - fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }else{
	//minus strand
	if((strcmp(fparam->strand, "*") == 0 || (strcmp(fparam->strand, "-") == 0))
	   && (fparam->start + fparam->shift) <= hit->core.pos 
	   && hit->core.pos <= (fparam->end + fparam->shift))
	{
	    int32_t w = get_inverse_weight(hit);
	    if(fparam->maxhit == 0 || w <= fparam->maxhit)
		fparam->sum += 1.0/w;
	}
    }
    return 0;
}

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
SEXP count_alignments(SEXP bam_in, SEXP idx_in, SEXP regions, SEXP stranded, 
		      SEXP overlap_type, SEXP min_overlap, SEXP shift, SEXP maxhit)
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

    // check parameter stranded, overlap type, shift, maxhit and minoverlap
    if(!IS_LOGICAL(stranded) || LENGTH(stranded) != 1)
	Rf_error("'stranded must be of type locical(1)");
    if(!IS_CHARACTER(overlap_type) || LENGTH(overlap_type) != 1)
	Rf_error("'type must be of type character(1)");
    if(!IS_INTEGER(shift) && LENGTH(shift) != num_files)
	Rf_error("'type must be of type integer and equal length as 'bam_in'");
    if(!IS_INTEGER(maxhit) || LENGTH(maxhit) != 1)
	Rf_error("'type must be of type integer(1)");
    if(!IS_INTEGER(min_overlap) || LENGTH(min_overlap) != 1)
	Rf_error("'type must be of type integer(1)");
    
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

    // initialise fetch_param with constant values
    fetch_param fparam;
    fparam.minoverlap = INTEGER(min_overlap)[0];
    fparam.shift = 0;
    fparam.maxhit = INTEGER(maxhit)[0];
    // ev. TODO if stranded false then strand always *, do once not for each region

 
    // chose fetch function
    bam_fetch_f fetch_func = 0;
    switch(translateChar(STRING_ELT(overlap_type, 0))[0]){
    case 'a': // any
	if(INTEGER(shift)[0] == 0 && !LOGICAL(stranded)[0] 
	   && INTEGER(min_overlap)[0] == 1 && INTEGER(maxhit)[0] < 1)
	    fetch_func = _fetch_any;
	else
	    fetch_func = _fetch_any_shift;
	break;
    case 'w': // within
	fetch_func = _fetch_within_shift;
	break;
    case 's': // startWithin
	fetch_func = _fetch_startwithin_shift;
	break;
    case 'e': // endWithin
	fetch_func =  _fetch_endwithin_shift;
      	break;
    default:
	Rf_error("The value of 'overlap_type' not supportet.");
	break;
    }
    
    // run bam_fetch
    SEXP count;
    PROTECT(count = NEW_NUMERIC(num_regions));
    for(int i = 0; i < num_regions; i++){
	fparam.sum = 0.0;
	fparam.start = INTEGER(start)[i];
	fparam.end = INTEGER(end)[i];
	if(LOGICAL(stranded)[0])
	    fparam.strand = translateChar(STRING_ELT(strand, i));
	else
	    fparam.strand = "*";
	
	for(int j = 0; j < num_files; j++){
	    fparam.shift = INTEGER(shift)[j]; // TODO ev. move outside of loop
	    bam_fetch(fin[j]->x.bam, idx[j], INTEGER(tid)[i], 
		      (INTEGER(start)[i]-abs(INTEGER(shift)[j])), 
		      INTEGER(end)[i]+abs(INTEGER(shift)[j]), 
		      &fparam, fetch_func); // start-1 because bam_fetch is exclusiv end of reads
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
    // TODO ev. free mem of fparam

    return count;
}
