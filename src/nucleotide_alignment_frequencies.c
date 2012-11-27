#include "nucleotide_alignment_frequencies.h"
#include "extract_unmapped_reads.h"
#include <stdlib.h>


/*! @typedef
  @abstract Structure to provid the data to the bam_fetch() functions.
  @field mm_dist    Vector containing the frequencies
  @field start      Start of the reference sequence in the genome
  @field end        End of the reference sequence in the genome
  @field ref        The reference sequence
  @field len        Maximal length of the reads
 */
typedef struct {
    int *mm_dist;
    int *mm_dist2;
    unsigned int mm_dist_len;
    int *frag_dist;
    unsigned int frag_dist_len;
    int start;
    int end;
    const char * ref;
    int len;
} fetch_param;

/*! @function
  @abstract  callback for bam_fetch() to count the nucleotide alignment frequencies.
  @param  b     the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _nucleotide_alignment_frequencies(const bam1_t *b, void *data)
{
    fetch_param *fparam = (fetch_param*)data;
    int *mm_dist;
    // separate mm_dist vector for read 1 and 2
    if ((b->core.flag & BAM_FREAD2) == 0)
	mm_dist = fparam->mm_dist;
    else
	mm_dist = fparam->mm_dist2;
    int start = fparam->start;
    int end = fparam->end;
    const char *ref = fparam->ref;

    if(start <= (int)b->core.pos && ((int)bam_calend(&b->core, bam1_cigar(b))) < end){
	uint8_t *seq = bam1_seq(b);
	uint32_t *cigar = bam1_cigar(b);
	int i = 0; // current position of cigar operations
	int j; // j count of current cigar operation
	int l; // l number of current cigar operation
	int x; // leftmost coordinate of the current cigar operation in reference string
	int y; // leftmost coordinate of the current cigar operation in read
	int z; // current position of the read sequence
	int u = 0; // count of matched position since last mismatch
	int op; // current cigar operation type
	int c1, c2; // current base int value of reference and read
	int32_t nm = 0; // number of mismatch
	int qlen = 0; // length of the query sequence
	int isize = 0; // fragment size
    
	// 4-bit encoding to index value
	// NACNGNNNTNNNNNNN -> A=0, C=1, G=2, T=3, N=4
	static int bit2idx[16] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

	// get length of the query sequence
	qlen = bam_cigar2qlen(&b->core, cigar);
	if(fparam->len < qlen)
	    fparam->len = qlen;
	qlen = qlen - 1; // vector index starts with 0

	// get length of the fragment
	isize = abs((int)b->core.isize);
	if((b->core.flag & BAM_FREAD1) && isize != 0){
	    if(isize < (fparam->frag_dist_len)) // check vector overflow
	    	fparam->frag_dist[isize - 1] += 1;
	    else
	    	// fragment size is bigger than vector length 
		// use last vector element for counting 
	    	fparam->frag_dist[fparam->frag_dist_len - 1] += 1; 
	}

	for (i = y = 0, x = b->core.pos-start; i < b->core.n_cigar; ++i) {
	    l = cigar[i]>>4, op = cigar[i]&0xf;
	    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
		for (j = 0; j < l; ++j) {
		    z = y + j;
		    c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
		    if (ref[x+j] == 0) break; // out of boundary
		    // calculate index position in the mm_dist vector
		    // c1 und c2 in the 4-bit encoding space
		    // A->1, C->2, G->4, T->8
		    if (bam1_strand(b) == 0)
			mm_dist[ bit2idx[c2] + 5*bit2idx[c1] + 25*z ] += 1;
		    else
			mm_dist[ bit2idx[c2] + 5*bit2idx[c1] + 25*(qlen-z) ] += 1;
		    // check if match or mismatch
		    if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
			++u;
		    } else {
			u = 0; ++nm;
		    }
		}
		if (j < l) break;
		x += l; y += l;
	    } else if (op == BAM_CDEL) {
		for (j = 0; j < l; ++j) {
		    if (ref[x+j] == 0) break;
		}
		u = 0;
		if (j < l) break;
		x += l; nm += l;
	    } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
		y += l;
		if (op == BAM_CINS) nm += l;
	    } else if (op == BAM_CREF_SKIP) {
		x += l;
	    }
	}
    }

    return 0;
}

/*! @function
  @abstract  Counts the nucleotide alignments frequencies by cycle in regions.
  @param  bamfile        Name of the bamfile
  @param  refSequence    Sequence string
  @param  refChr         tid of the reference sequence
  @param  refStart       Start coordination of the sequence string
  @param  mmDist         Vector containing the frequencies
  @return                Vector of the frequencies
 */
SEXP nucleotide_alignment_frequencies(SEXP bamfile, SEXP refSequence, SEXP refChr, SEXP refStart, SEXP mmDist)
{
    // check bamfile parameter
    if(!Rf_isString(bamfile) || Rf_length(bamfile) != 1)
        Rf_error("'bamfile' must be character(1)");
    const char * bam_in = translateChar(STRING_ELT(bamfile, 0));

    // check reference parameters
    if(!Rf_isString(refSequence) || Rf_length(refSequence) != 1)
        Rf_error("'refSequence' must be character(1)");
    const char *ref = translateChar(STRING_ELT(refSequence, 0));
    if(!Rf_isInteger(refChr) || Rf_length(refChr) != 1)
        Rf_error("'refChr' must be integer(1)");
    int tid = INTEGER(refChr)[0];
    if(!Rf_isInteger(refStart) && Rf_length(refStart) != 1)
	Rf_error("'refStart' must be of type integer(1)");
    int start = INTEGER(refStart)[0] - 1; // 0-based start
    int end = start + strlen(ref);

    // check mmDist parameter
    if(!Rf_isInteger(VECTOR_ELT(mmDist,0))){
        Rf_error("'mmDist[[0]]' must be integer");
    }
    if(!Rf_isInteger(VECTOR_ELT(mmDist,1))){
        Rf_error("'mmDist[[1]]' must be integer");
    }
    if(!Rf_isInteger(VECTOR_ELT(mmDist,2))){
        Rf_error("'mmDist[[2]]' must be integer");
    }
    //int *mm_dist = INTEGER(mmDist);

    // open bam file
    samfile_t *fp = 0;
    fp = samopen(bam_in, "rb", NULL);
    if (fp == 0)
	Rf_error("failed to open BAM file: '%s'", bam_in);
    if (fp->header == 0 || fp->header->n_targets == 0)
	Rf_error("BAM header missing or empty of file: '%s'", bam_in);

    // open index file
    bam_index_t *idx = 0;
    idx = bam_index_load(bam_in);
    if (idx == 0){
	samclose(fp);
	Rf_error("failed to open BAM index file: '%s.bai'", bam_in);
    }

    // create fetch data structure
    fetch_param fparam;
    fparam.mm_dist = INTEGER(VECTOR_ELT(mmDist,0));
    fparam.mm_dist2 = INTEGER(VECTOR_ELT(mmDist,1));
    fparam.mm_dist_len = Rf_length(VECTOR_ELT(mmDist,0));
    fparam.frag_dist = INTEGER(VECTOR_ELT(mmDist,2));
    fparam.frag_dist_len = Rf_length(VECTOR_ELT(mmDist,2));
    fparam.start = start;
    fparam.end = end;
    fparam.ref = ref;
    fparam.len = 0;

    // set fetch function
    bam_fetch_f fetch_func = _nucleotide_alignment_frequencies;

    // execute fetch
    bam_fetch(fp->x.bam, idx, tid, start, end, &fparam, fetch_func);

    // clean up
    samclose(fp);
    bam_index_destroy(idx);

    if(fparam.len == 0)
	return  Rf_ScalarInteger(0);
    else
	return Rf_ScalarInteger(fparam.len);
}
