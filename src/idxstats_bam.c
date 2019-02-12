/*!

Retrieve and print stats in the index file.
This code is derived from the bam_idxstats function in samtools.

 */

#include "idxstats_bam.h"

SEXP idxstats_bam(SEXP inBam)
{
    if (!Rf_isString(inBam) || 1 != Rf_length(inBam)){
	Rf_error("'inBam' must be character(1)");
    }
    
    const char * in_bam =  Rf_translateChar(STRING_ELT(inBam, 0));
    
    bam_index_t *idx;
    bam_header_t *header;
    bamFile fp;
    int i;
    
    fp = bam_open(in_bam, "r");
    if (fp == 0) { 
	Rf_error("[%s] fail to open BAM.\n", __func__); 
    }
    header = bam_header_read(fp);
    bam_close(fp);
    
    idx = bam_index_load(in_bam);
    if (idx == 0) { 
	Rf_error("[%s] fail to load the index.\n", __func__); 
    }
    
    SEXP stats, attrib, refnames, reflengths, mapped, unmapped;
    PROTECT(stats = allocVector(VECSXP, 4));
    PROTECT(refnames = allocVector(STRSXP, hts_idx_get_n(idx)+1));
    PROTECT(reflengths = allocVector(INTSXP, hts_idx_get_n(idx)+1));
    PROTECT(mapped = allocVector(INTSXP, hts_idx_get_n(idx)+1));
    PROTECT(unmapped = allocVector(INTSXP, hts_idx_get_n(idx)+1));
    PROTECT(attrib = allocVector(STRSXP, 4));
    
    for (i = 0; i < hts_idx_get_n(idx); ++i) {
        uint64_t mapped_i, unmapped_i;
	SET_STRING_ELT(refnames, i, mkChar(header->target_name[i]));
	INTEGER(reflengths)[i] = (int)(header->target_len[i]);
        hts_idx_get_stat(idx, i, &mapped_i, &unmapped_i);
	INTEGER(mapped)[i] = (int) mapped_i;
        INTEGER(unmapped)[i] = (int) unmapped_i;
    }
    SET_STRING_ELT(refnames, i, mkChar("*"));
    INTEGER(reflengths)[i] = 0;
    INTEGER(mapped)[i] = 0;
    INTEGER(unmapped)[i] = (int) hts_idx_get_n_no_coor(idx);
    
    SET_STRING_ELT(attrib, 0, mkChar("seqname"));
    SET_STRING_ELT(attrib, 1, mkChar("seqlength"));
    SET_STRING_ELT(attrib, 2, mkChar("mapped"));
    SET_STRING_ELT(attrib, 3, mkChar("unmapped"));
    
    SET_VECTOR_ELT(stats, 0, refnames);
    SET_VECTOR_ELT(stats, 1, reflengths);
    SET_VECTOR_ELT(stats, 2, mapped);
    SET_VECTOR_ELT(stats, 3, unmapped);
    setAttrib(stats, R_NamesSymbol, attrib);
    
    bam_header_destroy(header);
    bam_index_destroy(idx);
    
    UNPROTECT(6);
    
    return stats;
}
