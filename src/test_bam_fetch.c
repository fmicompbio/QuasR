/*!
This funciton is used to test the coordinate system in bam files and also the coordinates 
used when performing fetch. It corresponds to the same convention as bed files:
zero based end non-inclusive.

When fetching from bam files, the coordinates that are returned are:
startReadBam = hit-core.pos
endReadBam = (int)bam_calend(&hit->core, bam1_cigar(hit))


One based genome (UCSC browser)
1 2 3 4 5 6 7 8 9
    ---      read (3-4, 2 nt length)

      -----  region of interest (4-6, 3 nt length)

Bam file (startReadBam,endReadBam):
0 1 2 3 4 5 6 7 8
    ---      read (2-4, 2 nt length)

      -----  region of interest (3-6, 3 nt length)

Fetch expamples (startReadBam,endReadBam)
0 1 2 3 4 5 6 7 8
  ---             (1-3)
      -----       (3-6)  -> no hit

    ---           (2-4)
      -----       (3-6)  -> hit

            ---   (6-8)
      -----       (3-6)  -> no hit

          ---     (5-7)
      -----       (3-6)  -> hit

 */

#include "count_alignments.h"
#include "utilities.h"
#include <stdlib.h>
#include <time.h>

typedef struct {
    int start; // offset
    int end;// width
} fetch_param;


static int _bam_fetch_function(const bam1_t *hit, void *data){

    fetch_param *fparam = (fetch_param*)data;
	Rprintf("%i %i %i\n", hit->core.pos, (int)bam_calend(&hit->core, bam1_cigar(hit)), hit->core.flag & BAM_FREVERSE );
	
    return 0;
}


SEXP test_bam_fetch(SEXP bam_in, SEXP tid, SEXP start, SEXP end)
{

    // check bam_in 
    if(!IS_CHARACTER(bam_in) || LENGTH(bam_in) != 1)
        Rf_error("'bam_in' must be character(1)");

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
    idx = bam_index_load(translateChar(STRING_ELT(bam_in, 0)));
    if (idx == 0){
	samclose(fin);
	Rf_error("failed to open BAM index file: '%s'", translateChar(STRING_ELT(bam_in, 0)));
    }

    // check parameter stranded, overlap type, shift and minoverlap
    if(!IS_INTEGER(tid) && LENGTH(tid) != 1)
	Rf_error("'tid' must be of type integer(1)");
    if(!IS_INTEGER(start) && LENGTH(start) != 1)
	Rf_error("'start' must be of type integer(1)");
    if(!IS_INTEGER(end) && LENGTH(end) != 1)
	Rf_error("'end' must be of type integer(1)");
 
    // initialise fetch_param with constant values
    fetch_param fparam;
    fparam.start = INTEGER(start)[0];
    fparam.end = INTEGER(end)[0];
     
    // set fetch function
    bam_fetch_f fetch_func = _bam_fetch_function;

    // run fetch
    bam_fetch(fin->x.bam, idx, INTEGER(tid)[0], INTEGER(start)[0], INTEGER(end)[0], &fparam, fetch_func);

    samclose(fin);
    bam_index_destroy(idx);

    return Rf_ScalarInteger(0);
}
