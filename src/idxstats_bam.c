/*!

Retrieve and print stats in the index file.
This code is derived from the bam_idxstats function in samtools.

 */

#include "idxstats_bam.h"

#define BAM_MAX_BIN 37450 // =(8^6-1)/7+1

typedef struct {
    uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
//KSORT_INIT(off, pair64_t, pair64_lt)

typedef struct {
    uint32_t m, n;
    pair64_t *list;
} bam_binlist_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} bam_lidx_t;


KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct __bam_index_t {
    int32_t n;
    uint64_t n_no_coor; // unmapped reads without coordinate
    khash_t(i) **index;
    bam_lidx_t *index2;
};


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
    PROTECT(refnames = allocVector(STRSXP, idx->n+1));
    PROTECT(reflengths = allocVector(INTSXP, idx->n+1));
    PROTECT(mapped = allocVector(INTSXP, idx->n+1));
    PROTECT(unmapped = allocVector(INTSXP, idx->n+1));
    PROTECT(attrib = allocVector(STRSXP, 4));
    
    for (i = 0; i < idx->n; ++i) {
	khint_t k;
	khash_t(i) *h = idx->index[i];
	SET_STRING_ELT(refnames, i, mkChar(header->target_name[i]));
	INTEGER(reflengths)[i] = (int)(header->target_len[i]);
	k = kh_get(i, h, BAM_MAX_BIN);
	if (k != kh_end(h)){
	    INTEGER(mapped)[i] = (int)kh_val(h, k).list[1].u;
	    INTEGER(unmapped)[i] = (int)kh_val(h, k).list[1].v;
	} else  {
	    INTEGER(mapped)[i] = 0;
	    INTEGER(unmapped)[i] = 0;		
	}
    }
    SET_STRING_ELT(refnames, i, mkChar("*"));
    INTEGER(reflengths)[i] = 0;
    INTEGER(mapped)[i] = 0;
    INTEGER(unmapped)[i] = (int)idx->n_no_coor;	
    
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
