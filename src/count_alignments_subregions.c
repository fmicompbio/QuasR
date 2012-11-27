/*!
  @header

  Counts alignments in a given set of regions which are located in a subspace of the genome
 
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "count_alignments_subregions.h"
#include "utilities.h"
#include <stdlib.h>

/*! @typedef
  @abstract Structure to provid the data to the bam_fetch() functions.
  @field cov_plus   Coverage vector for plus strand
  @field cov_minus  Coverage vector for minus strand
  @field start      Start of the fetch region
  @field end        End of the fetch region
 */
typedef struct {
    int* cov_plus;
    int* cov_minus;
    int start; // offset
    int end;   // offset+width
    int shift; // shift of the reads
} fetch_param;


/*! @function
  @abstract  callback for bam_fetch() to ....
  @param  b     the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _add_start_to_coverage_vector(const bam1_t *hit, void *data)
{
    fetch_param *fparam = (fetch_param*)data;
    static int start_pos = 0;

    if(((hit->core.flag & BAM_FREVERSE) != 16)){
	// plus strand, start position of the read on left side
	start_pos = (int)hit->core.pos + fparam->shift;
	if(fparam->start <= start_pos && start_pos < fparam->end)
	    fparam->cov_plus[1 + start_pos - fparam->start] += 1;
    }else{
	// minus strand, start position of the read on right side
	start_pos = (int)bam_calend(&hit->core, bam1_cigar(hit)) - fparam->shift;
	if(fparam->start < start_pos && start_pos <= fparam->end)
	    fparam->cov_minus[start_pos - fparam->start] += 1;
    }

    return 0;
}

/*! @function
  @abstract  callback for bam_fetch() to ....
  @param  b     the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _add_end_to_coverage_vector(const bam1_t *hit, void *data)
{
    fetch_param *fparam = (fetch_param*)data;
    static int end_pos = 0;

    if(((hit->core.flag & BAM_FREVERSE) != 16)){
	// plus strand, end position of the read on right side
	end_pos = (int)bam_calend(&hit->core, bam1_cigar(hit)) + fparam->shift;
	if(fparam->start < end_pos && end_pos <= fparam->end)
	    fparam->cov_plus[end_pos - fparam->start] += 1;
    }else{
	// minus strand, end position of the read on left side
	end_pos = (int)hit->core.pos - fparam->shift;
	if(fparam->start <= end_pos && end_pos < fparam->end)
	    fparam->cov_minus[1 + end_pos - fparam->start] += 1;
    }

    return 0;
}

/*! @function
  @abstract  callback for bam_fetch() to...
  @param  b     the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _add_mid_to_coverage_vector(const bam1_t *hit, void *data)
{
    fetch_param *fparam = (fetch_param*)data;
    static int mid_pos = 0;

    if(hit->core.isize > 0){
	// leftmost fragment
        mid_pos = (int)floor((double)hit->core.pos + ((double)hit->core.isize-1)/2);
	if(fparam->start <= mid_pos && mid_pos < fparam->end)
	    fparam->cov_plus[1 + mid_pos - fparam->start] += 1;

    } else if(hit->core.isize < 0){
	// rightmost fragment
        mid_pos = (int)floor(((double)bam_calend(&hit->core, bam1_cigar(hit))) + ((double)hit->core.isize-1)/2);
	if(fparam->start <= mid_pos && mid_pos < fparam->end)
	    fparam->cov_minus[1 + mid_pos - fparam->start] += 1;
    } else {
	// insert size is zero. don't count
    }

    return 0;
}

/*! @function
  @abstract  Counts alignments in a given set of regions which are located in a subspace of the genome
  @param  bam_in        Name of the bamfile
  @param  idx_in        Name of the bam index file without the '.bai' extension
  @param  tid           Reference sequence name identifier of the bamfile
  @param  min_regions   Minimum coordinate of regions
  @param  max_regions   Maximum coordinate of regions
  @param  regions       Coordinates of fetch regions
  @param  shift         Shift size of the reads 
  @param  broaden       Broaden size of the regions
  @param  overlap_type  Type of the overlap criterion
  @return               Vector of the alignment counts
 */
SEXP count_alignments_subregions(SEXP bam_in, SEXP idx_in, SEXP tid,  SEXP min_regions, SEXP max_regions, SEXP regions, SEXP shift, SEXP broaden, SEXP overlap_type)
{
    // check bam_in and idx_in parameters
    if(!IS_CHARACTER(bam_in) || LENGTH(bam_in) != 1)
        Rf_error("'bam_in' must be character(1)");
    if(!IS_CHARACTER(idx_in)|| LENGTH(idx_in) != 1)
        Rf_error("'idx_in' must be character(1)");

    // open bam file
    samfile_t *fin = 0;
    fin = samopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);
    if (fin == 0)
	Rf_error("failed to open BAM file: '%s'", translateChar(STRING_ELT(bam_in, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
	samclose(fin);
	Rf_error("BAM header missing or empty of file: '%s'", translateChar(STRING_ELT(bam_in, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(translateChar(STRING_ELT(idx_in, 0)));
    if (idx == 0){
	samclose(fin);
	Rf_error("failed to open BAM index file: '%s'", translateChar(STRING_ELT(bam_in, 0)));
    }

    // check parameter stranded, overlap type, shift and minoverlap
    if(!IS_CHARACTER(overlap_type) || LENGTH(overlap_type) != 1)
        Rf_error("'overlap_type' must be of type character(1)");
    if(!IS_INTEGER(tid) && LENGTH(tid) != 1)
        Rf_error("'tid' must be of type integer(1)");
    if(!IS_INTEGER(min_regions) && LENGTH(min_regions) != 1)
        Rf_error("'min_regions' must be of type integer(1)");
    if(!IS_INTEGER(max_regions) && LENGTH(max_regions) != 1)
        Rf_error("'max_regions' must be of type integer(1)");
    if(!IS_INTEGER(shift) && LENGTH(shift) != 1)
        Rf_error("'shift' must be integer(1)");
    if(!IS_INTEGER(broaden) && LENGTH(broaden) != 1)
        Rf_error("'broaden' must be integer(1)");
    if(INTEGER(broaden)[0] < 0)
        Rf_error("'broaden' must be a positive value.");
        
    // check parameter region and get direct pointer to the elements
    SEXP start = _getListElement(regions, "start");
    SEXP end = _getListElement(regions, "end");
    SEXP strand = _getListElement(regions, "strand");
    if(!IS_INTEGER(start))
        Rf_error("Column 'start' must be of type integer");
    if(!IS_INTEGER(end))
        Rf_error("Column 'end' must be of type integer");
    if(!IS_CHARACTER(strand))
        Rf_error("Column 'strand' must be of type character");
    int num_regions = LENGTH(start);
    if( num_regions != LENGTH(end) || num_regions != LENGTH(strand))
        Rf_error("The columns 'start', 'end', 'stand' must have equal length.");
 
    // chose fetch function
    bam_fetch_f fetch_func = 0;
    switch(translateChar(STRING_ELT(overlap_type, 0))[0]){
    case 's': // startWithin
	fetch_func = _add_start_to_coverage_vector;
	break;
    case 'e': // endWithin
	fetch_func =  _add_end_to_coverage_vector;
      	break;
    case 'm': // midWithin
	fetch_func = _add_mid_to_coverage_vector;
	break;
    default:
	Rf_error("The value of 'overlap_type' not supportet.");
	break;
    }

    // set up return parameter
    SEXP cnt;
    PROTECT(cnt = allocVector(INTSXP, num_regions));

    // set up fetch region coordinats
    int cov_start = INTEGER(min_regions)[0]; // offset
    int cov_end = INTEGER(max_regions)[0];
    int width = cov_end - cov_start + 1; // first vector position is for the initial zero in the cumsum vector

    // set up coverage vectors
    int *cov_plus = (int*) R_Calloc(width, int);
    int *cov_minus = (int*) R_Calloc(width, int);

    // initialise fetch_param with constant values
    fetch_param fparam;
    fparam.cov_plus = cov_plus;
    fparam.cov_minus = cov_minus;
    fparam.start = cov_start;
    fparam.end = cov_end;
    fparam.shift = INTEGER(shift)[0];

    // run fetch
    bam_fetch(fin->x.bam, idx, 
	      INTEGER(tid)[0], 
	      cov_start - INTEGER(broaden)[0] - abs(INTEGER(shift)[0]), 
	      cov_end + INTEGER(broaden)[0] + abs(INTEGER(shift)[0]), 
	      &fparam, fetch_func);
    
    // run cumsum
    unsigned int i = 0;
    for(i=1; i<width; i++){
	cov_plus[i] = cov_plus[i-1] + cov_plus[i];
	cov_minus[i] = cov_minus[i-1] + cov_minus[i];
    }

    // get counts per region
    static int cnt_plus, cnt_minus;
    for(i=0; i < num_regions; i++){
        cnt_plus = cov_plus[ INTEGER(end)[i] - cov_start ] - cov_plus[ INTEGER(start)[i] - cov_start ];
        cnt_minus = cov_minus[ INTEGER(end)[i] - cov_start ] - cov_minus[ INTEGER(start)[i] - cov_start ];
	switch(translateChar(STRING_ELT(strand, i))[0]){
	case '+': 
	    INTEGER(cnt)[i] = cnt_plus;
	    break;
	case '-': 
	    INTEGER(cnt)[i] = cnt_minus;
	    break;
	default:    
	    INTEGER(cnt)[i] = cnt_plus + cnt_minus;
	    break;
	}
    }

    // clean up
    samclose(fin);
    bam_index_destroy(idx);
    R_Free(cov_plus);
    R_Free(cov_minus);

    UNPROTECT(1);
    return cnt;
}
