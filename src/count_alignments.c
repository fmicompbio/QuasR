
#include "quasr.h"

const int MAXHITS = 32766; // Allowed upper range for signed short int minus 1

SEXP
get_allowed_max_hits(){
  return ScalarInteger(MAXHITS);
}

samfile_t *
_bam_tryopen(const char *filename, const char *filemode, void *aux)
{
    samfile_t *sfile = samopen(filename, filemode, aux);
    if (sfile == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", 
		 filename);
    if (sfile->header == 0 || sfile->header->n_targets == 0) {
        samclose(sfile);
        Rf_error("SAM/BAM header missing or empty\n  file: '%s'", 
                 filename);
    }
    return sfile;
}

int
_write_buffered_alignment(samfile_t *fout, bam1_t **fifo, int r_idx, int w_idx, int32_t fifo_size, int32_t c)
{
Rprintf("\nWrite buffer r%i w%i c%i size%i",r_idx, w_idx, c, fifo_size); 
  while(w_idx != r_idx){
    Rprintf("-qname %s w%i c%i-",bam1_qname(fifo[w_idx]), w_idx, c);
    bam_aux_append(fifo[w_idx], "IH", 'i', 4, (uint8_t*)&c);         
    samwrite(fout, fifo[w_idx]);
    w_idx = (w_idx + 1) % fifo_size;
    
  }
  return w_idx;
}

int
_count_alignments(samfile_t *fin, samfile_t *fout, int max_hits)
{
  if (max_hits > MAXHITS){
    max_hits = MAXHITS;
  }
  int fifo_size = max_hits+1;
  int32_t k= 0;
	int r, rd, w, count = 0;
	const char *qname = "";
	
	bam1_t **fifo;
  fifo = (bam1_t**)calloc(fifo_size, sizeof(bam1_t*));
  for(int i = 0; i < fifo_size; i++)	
    fifo[i] = bam_init1();

Rprintf("\nAfter alloc r%i w%i c%i size%i",rd, w, k, max_hits);  
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
      Rprintf("\nQname %s",qname); 
    }
    //check if there are more alignments per query then max_hits
    if (k >= max_hits){
      //TODO
//Rprintf("\nMax hits %s r%i w%i c%i size%i", qname, rd, w, k, max_hits); 
Rprintf("\nError: Max Hits exceeded");
	    w = _write_buffered_alignment(fout, fifo, rd, w, fifo_size, max_hits);
	    //return (-1 * count);
	  } 
    rd = (rd + 1) % fifo_size;
    k++;  
	  count++;	  
  }
  w = _write_buffered_alignment(fout, fifo, rd, w, max_hits, k);
  
  for(int i = 0; i < fifo_size; i++)	
    bam_destroy1(fifo[i]);
  return r >= -1 ? count : -1 * count;
}

SEXP
count_alignments(SEXP bam_in, SEXP bam_out, SEXP max_hits)
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

  int count = _count_alignments(fin, fout, asInteger(max_hits));

  samclose(fin);
  samclose(fout);
  if (count < 0)
      Rf_error("truncated input file at record %d", -1 * count);

  return bam_out;
}

