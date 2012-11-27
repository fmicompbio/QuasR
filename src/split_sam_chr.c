#include "split_sam_chr.h"

/* // moved to utilities.c
samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux)
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
 */

int _walk_through_sam_and_split(samfile_t * fin, samfile_t **foutList)
{
    bam1_t *b = bam_init1();
    int r, count = 0;

    while (0 <= (r = samread(fin, b))) {
      if(b->core.tid > -1){
        samwrite(foutList[b->core.tid], b);
      }else{
        samwrite(foutList[fin->header->n_targets], b);
      }
      count++;
    }
    bam_destroy1(b);

    return r >= -1 ? count : -1 * count;
}

char * _assemble_file_name(const char *out_dir, const char *chr){
  char* full_file_name = malloc(strlen(out_dir)+strlen(chr)+4+1+1);
  strcpy(full_file_name, out_dir);
  strcat(full_file_name, "/");
  strcat(full_file_name, chr);
  strcat(full_file_name, ".sam");
  
  return full_file_name;
}


// Splits a sam file into individual files, one per chromosome. The files are created in the specified directory.
// Each splitted file contains the same header as the original samFile. The generated files may only contain
// a header if no alignments to that chromosome exist. The names of the split files come from the sam header
// with an additional .sam extension. The unmapped alignments are collected in the file splitChrSam_unaligned.sam
// Returns the chromosome names in the order in which they occur in the sam file header
SEXP split_sam_chr(SEXP samFile, SEXP outDir)
{
  if (!Rf_isString(samFile) || 1 != Rf_length(samFile)){
    Rf_error("'samFile' must be character(1)");
  }

  if (!Rf_isString(outDir) || 1 != Rf_length(outDir)){
    Rf_error("'outDir' must be character(1)");
  }

  const char * sam_file =  Rf_translateChar(STRING_ELT(samFile, 0));
  const char * out_dir =  Rf_translateChar(STRING_ELT(outDir, 0));

  // open the input sam file
  samfile_t *fin = _bam_tryopen(sam_file, "r", NULL);
  if (fin->header == 0) {
    samclose(fin);
    Rf_error("invalid header");
  }

  // remove \r from header if exists (for windows)
  int j, k = 0;
  for(j = 0; j<fin->header->l_text; j++){
    if(fin->header->text[j] != '\r'){
      fin->header->text[k++] = fin->header->text[j];
    }
  }
  if(j != k){
    fin->header->text[k] = '\0';
    fin->header->l_text = strlen(fin->header->text);
  }

  // allocate memory for a list of filehandles (n+1 because of the unaligned reads)
  samfile_t **foutList = (samfile_t**)calloc((fin->header->n_targets+1), sizeof(samfile_t*));

  // open the output file handles (n+1 due to the unaligned reads)
  int i;
  SEXP chrNames;
  PROTECT(chrNames = allocVector(STRSXP, (fin->header->n_targets+1))); // protect from garbage collector

  for (i = 0; i < (fin->header->n_targets); i++) {
    foutList[i] = _bam_tryopen(_assemble_file_name(out_dir,fin->header->target_name[i]), "wh", fin->header);
    SET_STRING_ELT(chrNames, i, mkChar(fin->header->target_name[i]));
  }
  foutList[fin->header->n_targets] = _bam_tryopen(_assemble_file_name(out_dir,"splitChrSam_unaligned"), "wh", fin->header);
  SET_STRING_ELT(chrNames, fin->header->n_targets, mkChar("splitChrSam_unaligned"));

  // split the sam file based on chromosome
  _walk_through_sam_and_split(fin,foutList);

  // close all the file handles
  for (i = 0; i < (fin->header->n_targets+1); i++){samclose(foutList[i]);}
  samclose(fin);

  UNPROTECT(1); // release
  return chrNames;
}

