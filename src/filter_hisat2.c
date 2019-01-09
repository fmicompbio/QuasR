#include "filter_hisat2.h"

/*
 * parse a SAM/BAM file (record by record)
 * and remove any secondary alignments, as well as
 * alignments for reads with more than maxHits alignments
 * (as obtained from the NH tag)
 * 
 * arguments:
 *   samFile : input SAM/BAM file to filter (must contain header)
 *   outFile : output SAM/BAM file to write to
 *   maxHits : integer scalar with the maximum allowed hits per read
 *   
 * return value:
 *   integer(2) SEXP with
 *   - the number of suppressed secondary alignments
 *   - the number of suppressed alignments with more than maxHits hits
 *
 * for documentation on samtools C API, see:
 *   http://samtools.sourceforge.net/samtools/bam/
 *   http://samtools.sourceforge.net/samtools/sam/ 
 *   https://github.com/lh3/samtools/blob/master/bam.h
 */ 
// 
//  
SEXP filter_hisat2(SEXP samFile, SEXP outFile, SEXP maxHits) {
  if (!Rf_isString(samFile) || 1 != Rf_length(samFile))
    Rf_error("'samFile' must be character(1)");
  if (!Rf_isString(outFile) || 1 != Rf_length(outFile))
    Rf_error("'outFile' must be character(1)");
  if (!Rf_isInteger(maxHits) || 1 != Rf_length(maxHits))
    Rf_error("'maxhits' must be integer(1)");

  int max_hits = Rf_asInteger(maxHits), n_secondary = 0, n_overmapped = 0;    
  const char * sam_file =  Rf_translateChar(STRING_ELT(samFile, 0));
  const char * out_file =  Rf_translateChar(STRING_ELT(outFile, 0));

  // open the input file
  samfile_t *fin = _bam_tryopen(sam_file, "r", NULL);

  // remove \r from header if exists (for windows)
  int j, k = 0;
  for(j = 0; j<fin->header->l_text; j++){
    if(fin->header->text[j] != '\r'){
      fin->header->text[k++] = fin->header->text[j];
    }
  }
  if(j != k){
    fin->header->text[k] = '\0';
    fin->header->l_text = (uint32_t)strlen(fin->header->text);
  }

  // open the output file handle
  samfile_t *fout = _bam_tryopen(out_file, "wh", fin->header);

  // walk through alignments and filter
  int32_t nh;
  bam1_t *b = bam_init1();
  bam1_t *bu = bam_init1(); // for unaligned reads
  bu->core.tid = -1;
  bu->core.pos = -1;
  bu->core.bin = 0;
  bu->core.qual = 0;
  bu->core.n_cigar = 0;
  bu->core.mtid = -1;
  bu->core.mpos = -1;
  bu->core.isize = 0;
  bu->l_aux = 0;
  
  // for each alignment ...
  while (0 <= samread(fin, b)) {
    
    // if the alignment is “primary”
    if (!(b->core.flag & BAM_FSECONDARY)) {
      
      // read the 'NH' tag
      nh = bam_aux2i(bam_aux_get(b, "NH"));
    
      if (nh >= max_hits + 1) {
        // primary alignment, but more than maxHits alignments
        // ... modify b to represent an unaligned read
        // ... ... flag
        bu->core.flag = b->core.flag;
        if (bu->core.flag & BAM_FPROPER_PAIR)
          bu->core.flag -= BAM_FPROPER_PAIR;
        if (!(bu->core.flag & BAM_FUNMAP))
          bu->core.flag += BAM_FUNMAP;
        if ((bu->core.flag & BAM_FPAIRED) & !(bu->core.flag & BAM_FMUNMAP))
          bu->core.flag += BAM_FMUNMAP;
        if (bu->core.flag & BAM_FREVERSE)
          bu->core.flag -= BAM_FREVERSE;
        if (bu->core.flag & BAM_FMREVERSE)
          bu->core.flag -= BAM_FMREVERSE;
        
        // ... ... length of qname and qseq
        bu->core.l_qname = b->core.l_qname;
        bu->core.l_qseq = b->core.l_qseq;
        
        // ... ... data
        // b->data structure: qname-cigar-seq-qual-aux, here: only need qname-seq-qual
        bu->data_len = b->core.l_qname + (b->core.l_qseq + 1)/2 + b->core.l_qseq; // new size of variable data (qname-seq-qual), see bam1_aux in bam.h
        if (bu->m_data < bu->data_len) { // allocate space if necessary
          bu->m_data = bu->data_len;
          kroundup32(bu->m_data);
          bu->data = (uint8_t*)realloc(bu->data, bu->m_data);
        }
        memcpy(bu->data, bam1_qname(b), b->core.l_qname); // copy qname
        memcpy(bam1_seq(bu), bam1_seq(b), (b->core.l_qseq + 1)/2 + b->core.l_qseq); // copy seq-qual

        // ... and output
        n_overmapped++;
        samwrite(fout, bu);

      } else {
        // output the alignment as is
        samwrite(fout, b);
      }
    } else {
      // secondary alignment; do nothing
      n_secondary++;
    }
  }
  
  // clean up and close file handles
  bam_destroy1(b);
  bam_destroy1(bu);
  samclose(fout);
  samclose(fin);
  
  // return counts
  SEXP counts, attrib;
  PROTECT(counts = Rf_allocVector(INTSXP, 2));
  PROTECT(attrib = Rf_allocVector(STRSXP, 2));
  INTEGER(counts)[0] = n_secondary;
  INTEGER(counts)[1] = n_overmapped;
  SET_STRING_ELT(attrib, 0, Rf_mkChar("n_secondary"));
  SET_STRING_ELT(attrib, 1, Rf_mkChar("n_overmapped"));
  Rf_setAttrib(counts, R_NamesSymbol, attrib);
  UNPROTECT(2);
  
  return counts;
}
