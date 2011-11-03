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
  @abstract  Write out buffered alignment to bam file
  @param  *fout      Pointer to the output SAM/BAM file handler 
  @param  **fifo     Pointer to the FIFO queue
  @param  r_idx      Index to the current tail of the queues
  @param  w_idx      Index to the current head of the queues
  @param  fifo_size  Size of the FIFO queue
  @param  c          Inverse weight of the alignments
  @return    Index to the current head of the queues
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

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
static int* _weight_alignments(samfile_t *fin, samfile_t *fout, int max_hits, int* hitcount)
{
    if (max_hits > MAXHITS){
        max_hits = MAXHITS;
    }
    int fifo_size = max_hits+2; //plus 2, we read one ahead and  that read and write index are never equal
    int32_t k = 0; //counts of a query
    int r = 0; //return value of samread
    int rd = 0; //read index, head of fifo queue
    int w = 0; //write index, tail of fifo queue
    int count = 0; //counter of samread
    const char *qname = ""; //current query name
    
    bam1_t **fifo;
    fifo = (bam1_t**)calloc(fifo_size, sizeof(bam1_t*));
    for(int i = 0; i < fifo_size; i++)	
        fifo[i] = bam_init1();

    // initialize the array
    for(int i = 0; i < max_hits+2; i++)
        hitcount[i] = 0;    
    
    while (0 <= (r = samread(fin, fifo[rd]))) {
	//check if same query
	if(strcmp(qname, bam1_qname(fifo[rd]))!= 0){
	    //check if unmapped 
	    if((fifo[w]->core.flag & BAM_FUNMAP) != 0){
		if(bam_aux2i(bam_aux_get(fifo[w], "XM")) == 0)
		    hitcount[0] += 1;
		else
		    hitcount[max_hits+1] += 1;
		w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, 0);
	    }else{
		hitcount[k] += 1;
		w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, k);
	    }
	    k = 0;
	    qname = bam1_qname(fifo[rd]);
	}	
	k++;  
	count++;
	rd = (rd + 1) % fifo_size;
	//check if there are more alignments per query then max_hits
	if (k > max_hits){
	    if(k = max_hits+1){
		Rf_warning("Max Hits %i exceeded", max_hits);
		hitcount[max_hits+1] += 1;
	    }
	    w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, MAXHITS); // MAXHITS instead of zero because it is an inverse weight
	} 

	  
    }
    //to empty fifo queue
    if((fifo[w]->core.flag & BAM_FUNMAP) != 0){
	hitcount[0] += 1;
        w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, 0); 
    } else {
	hitcount[k] += 1;
        w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, k);
    }  
    
    for(int i = 0; i < fifo_size; i++)	
        bam_destroy1(fifo[i]);

    if (r < -1)
        Rf_error("truncated input file at record %d", count);
    return hitcount;
}

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
SEXP weight_alignments(SEXP bam_in, SEXP bam_out, SEXP header, SEXP max_hits)
{
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
        Rf_error("'bam_in' must be character(1)");
    if (!IS_CHARACTER(bam_out) || 1 != LENGTH(bam_out))
        Rf_error("'bam_out' must be character(1)");
   if (!IS_CHARACTER(header) || 1 != LENGTH(header))
        Rf_error("'header' must be character(1)");
    
    samfile_t *fin = _bam_tryopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);
    if (fin->header == 0) {
        samclose(fin);
        Rf_error("invalid header");
    }
    samfile_t *hout = _bam_tryopen(translateChar(STRING_ELT(header, 0)), "r", NULL);
    if (hout->header == 0) {
        samclose(hout);
        Rf_error("invalid header");
    }  
    samfile_t *fout = _bam_tryopen(translateChar(STRING_ELT(bam_out, 0)), "wb", hout->header); // f_in leaks if this fails 

    //int hitcout[asInteger(max_hits)+2];
    SEXP hitcount;
    PROTECT(hitcount = NEW_INTEGER(asInteger(max_hits)+2));
    int* hit = _weight_alignments(fin, fout, asInteger(max_hits), INTEGER(hitcount));
    
    samclose(fin);
    samclose(fout);
    samclose(hout);
    UNPROTECT(1);

    return hitcount;
}
