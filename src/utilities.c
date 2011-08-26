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
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
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

/*! @function
  @abstract  todo
  @param  todo  todo
  @return    todo
 */
SEXP seqname(SEXP bam_in)
{
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
    Rf_error("'bam_in' must be character(1)");
    samfile_t *fin = _bam_tryopen(translateChar(STRING_ELT(bam_in, 0)), "rb", NULL);

    SEXP ans, tid, name, attrib;
    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(tid = allocVector(INTSXP, fin->header->n_targets));
    PROTECT(name = allocVector(STRSXP, fin->header->n_targets));
    PROTECT(attrib = allocVector(STRSXP, 2));
    
    for(int i=0; i < fin->header->n_targets; i++){
        INTEGER(tid)[i] = i;
        SET_STRING_ELT(name, i, mkChar(fin->header->target_name[i]));
    }
    
    SET_STRING_ELT(attrib, 0, mkChar("seqnames"));
    SET_STRING_ELT(attrib, 1, mkChar("tid"));
    
    SET_VECTOR_ELT(ans, 0, name);
    SET_VECTOR_ELT(ans, 1, tid);
    setAttrib(ans, R_NamesSymbol, attrib);
    
    UNPROTECT(4);
    samclose(fin);

    return ans;
}

/* get the list element named str, or return NULL */
     
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue;
    SEXP names = getAttrib(list, R_NamesSymbol);
    
    for (R_len_t i = 0; i < length(list); i++)
	if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    return elmt;
}
