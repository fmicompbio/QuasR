/*!
  @header

  TODO.
 
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "quasr.h"

/*! @abstract maximal allowed alignments per read resp. the maximal inverse weight
 */
static const int MAXHITS = 32766; // Allowed upper range for signed short int minus 1

/*! @function
  @abstract  Get the maximal allowed alignments per read resp. the maximal inverse weight.
  @return    the maximal inverse weight as integer
 */
SEXP get_allowed_max_hits(){
    return ScalarInteger(MAXHITS);
}

/*! @function
  @abstract  todo
  @param  fout  todo
  @param  todo  todo
  @return    todo
 */
static int _write_buffered_alignment(samfile_t *fout, bam1_t **fifo, int r_idx, int w_idx, int32_t fifo_size, int32_t c)
{
    uint8_t *ih_ptr;
    while(w_idx != r_idx){
        ih_ptr = bam_aux_get(fifo[w_idx],"IH");
        if(ih_ptr != 0)
            bam_aux_del(fifo[w_idx], ih_ptr);
    	bam_aux_append(fifo[w_idx], "IH", 'i', 4, (uint8_t*)&c);         
    	samwrite(fout, fifo[w_idx]);
    	w_idx = (w_idx + 1) % fifo_size;
    }
    return w_idx;
}

static int _weight_alignments(samfile_t *fin, samfile_t *fout, int max_hits)
{
    if (max_hits > MAXHITS){
	max_hits = MAXHITS;
    }
    int fifo_size = max_hits+1; //fifo size is plus 1 that read and write index are never equal
    int32_t k= 0; //counts of a query
    int r = 0; //return value of samread
    int rd = 0; //read index
    int w = 0; //write index
    int count = 0; //counter of samread
    const char *qname = ""; //current query name
    
    bam1_t **fifo;
    fifo = (bam1_t**)calloc(fifo_size, sizeof(bam1_t*));
    for(int i = 0; i < fifo_size; i++)	
	fifo[i] = bam_init1();
    
    while (0 <= (r = samread(fin, fifo[rd]))) {
	//check if same query
	if(strcmp(qname, bam1_qname(fifo[rd]))!= 0){
	    //check if unmapped 
	    if((fifo[w]->core.flag & BAM_FUNMAP) != 0){
		w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, 0);
	    }else{    
		w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, k);
	    }
	    k = 0;
	    qname = bam1_qname(fifo[rd]);
	}
	//check if there are more alignments per query then max_hits
	if (k >= max_hits){
	    //TODO
	    Rf_warning("Max Hits %i exceeded", max_hits);
	    w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, 0); //TODO count is 0 or discart
	} 
	rd = (rd + 1) % fifo_size;
	k++;  
	count++;	  
    }
    //to empty fifo queue
    if((fifo[w]->core.flag & BAM_FUNMAP) != 0){
	w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, 0);
    }else{    
	w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, k);
    }  
    
    for(int i = 0; i < fifo_size; i++)	
	bam_destroy1(fifo[i]);
    return r >= -1 ? count : -1 * count;
}

SEXP weight_alignments(SEXP bam_in, SEXP bam_out, SEXP max_hits)
{
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
	Rf_error("'bam_in' must be character(1)");
    if (!IS_CHARACTER(bam_out) || 1 != LENGTH(bam_out))
	Rf_error("'bam_out' must be character(1)");
    
    samfile_t *fin = _bam_tryopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);
    if (fin->header == 0) {
	samclose(fin);
	Rf_error("invalid header");
    }
    
    // TODO add manipulation to header
    samfile_t *fout = _bam_tryopen(translateChar(STRING_ELT(bam_out, 0)), "wb", fin->header); // f_in leaks if this fails 
    
    int count = _weight_alignments(fin, fout, asInteger(max_hits));
    
    samclose(fin);
    samclose(fout);
    if (count < 0)
	Rf_error("truncated input file at record %d", -1 * count);
    
    return bam_out;
}
