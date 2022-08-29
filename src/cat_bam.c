#include "cat_bam.h"

int bam_cat(int, char * const *, sam_hdr_t *, const char*, char *, int);

SEXP cat_bam(SEXP inbam, SEXP outbam)
{
    bam_hdr_t *h = NULL;
	int ret = 0;

	// check parameter inbam
    if(!Rf_isString(inbam)){
        Rf_error("'inbam' must be character()");
    }
    // number of input bamfiles
    int nfin = Rf_length(inbam);
    // convert/copy STRSXP to a const char **
    const char ** fin = (const char **)R_Calloc(nfin, const char *);
    for(int i = 0; i < nfin; i++){
    	fin[i] = Rf_translateChar(STRING_ELT(inbam, i));
    }

    // check parameter outbam
    if(!Rf_isString(outbam) || 1 != Rf_length(outbam)){
        Rf_error("'outbam' must be character(1)");
    }
    const char * fout = Rf_translateChar(STRING_ELT(outbam, 0));

    ret = bam_cat(nfin, (char * const *)fin, h, fout, NULL, 1);
    if (ret != 0)
        Rf_error("call to bam_cat() returned a non-zero value");

    R_Free(fin);
    return Rf_ScalarInteger(ret);
}
