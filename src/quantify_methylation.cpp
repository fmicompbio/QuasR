/*
  quantify methylation from bisulfite-seq bam files: count C-to-C and C-to-T events
*/

#include "quantify_methylation.h"

using namespace std;

/*
inline char comp(char element) {
    static const char charMapNoIUPAC[] = {
        'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', // A to M
        'N', 'N', 'N', 'N', 'N', 'N', 'A', 'A', 'N', 'N', 'N', 'N', 'N'  // N to Z
    };
    return charMapNoIUPAC[toupper(element) - 'A'];
}

inline char baseIntToChar(const int x) {
    static const int intMapNoIUPAC[] = {
	'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
    };
    return intMapNoIUPAC[ x % 15 ];
}
*/

typedef struct { // for use with addHitToCounts(), bam_fetch callback function of quantify_methylation()
    int *Tp;  // total      (plus)
    int *Mp;  // methylated (plus)
    int *Tm;  // total      (minus)
    int *Mm;  // methylated (minus)
    bool *op; // output     (plus)
    bool *om; // output     (minus)
    int offset; // region offset
} methCounters;

typedef struct { // for use with addHitToSNP(), bam_fetch callback function of detect_SNVs()
    // count (mis-)matches on opposite strand (on the strand that was not altered in the bisulfite conversion)
    unsigned int *match;
    unsigned int *total;
    bool *targetC;
    bool *targetG;
    int offset; // region offset
} snpCounters;

typedef struct { // for use with addHitToCountsAllele(), bam_fetch callback function of quantify_methylation_allele()
    int *Tp[3];  // total      (plus, R/U/A)
    int *Mp[3];  // methylated (plus, R/U/A)
    int *Tm[3];  // total      (minus, R/U/A)
    int *Mm[3];  // methylated (minus, R/U/A)
    bool *op; // output     (plus)
    bool *om; // output     (minus)
    int offset; // region offset
} methCountersAllele;

typedef struct { // for use with addHitToCountsSingleAlignments(), bam_fetch callback function of quantify_methylation_singleAlignments()
    //vector<int> aid;  // vector of alignment identifiers
    vector<string> aid;  // vector of alignment identifiers
    vector<int> Cid;  // vector of C (DNA base) identifiers
    vector<char> strand; // vector of strands
    vector<int> meth; // vector of methylation status (0 or 1)
    bool *op; // output     (plus)
    bool *om; // output     (minus)
    int offset; // region offset
} methCountersSingleAlignments;

const inline int alleleFlagToInt(char xv) {
    static int xvi;
    switch(xv) {
    case 'R':
	xvi = 0;
	break;
    case 'A':
	xvi = 2;
	break;
    default:
	xvi = 1;
    };
    return xvi;
}

static int addHitToCounts(const bam1_t *hit, void *data) { // bam_fetch callback function of quantify_methylation()
  // REMARKS:
  //   assume ungapped read-global alignment
  //   count events separately for +/- strands (collapse strands during output if necessary)
  //   hit->core.pos and count vectors are zero-based (add one during output)
  //   bam_calend() return zero-based, exclusive end
  static uint8_t *hitseq=NULL;
  static unsigned long int i=0;
  static unsigned long int iend=0;
  static int j=0;
  static methCounters *cnt = NULL;

  hitseq = bam1_seq(hit);
  cnt = (methCounters*) data;
  iend = bam_calend(&(hit->core), bam1_cigar(hit)) - cnt->offset;

  if ((hit->core.flag & BAM_FPROPER_PAIR) && (hit->core.isize > 0) && (iend > hit->core.mpos - cnt->offset))
      // left fragment of a paired alignment --> make sure iend does not overlap alignment of right fragment
      iend = hit->core.mpos - cnt->offset;

  if (hit->core.flag & BAM_FREVERSE) {       // alignment on minus strand (reads are reverse complemented, look for G-A mismatches)
      //Rprintf("\nminus strand alignment %d-%d (offset %d), id=%s\n", hit->core.pos+1, bam_calend(&(hit->core), bam1_cigar(hit)), cnt->offset, bam1_qname(hit));
    for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
      if(cnt->om[i]) {                   //  target base is 'G'
	  //char Twobit2base[] = {'X', 'A', 'C', 'X', 'G', 'X', 'X', 'X', 'T', 'X', 'X', 'X', 'X', 'X', 'X', 'N'};
	  //Rprintf("  adding to genomic position %d (read pos %d has %c)\n", i+cnt->offset+1, j+1, Twobit2base[bam1_seqi(hitseq, j)]);
	if(bam1_seqi(hitseq, j)==4) {        //  query base is 'G'
	  cnt->Tm[i]++;
	  cnt->Mm[i]++;
	} else if(bam1_seqi(hitseq, j)==1) { //  query base is 'A'
	  cnt->Tm[i]++;
	}
      }

  } else {                                   // alignment on plus strand (look for C-T mismatches)
      //Rprintf("\nplus strand alignment %d-%d (offset %d), id=%s\n", hit->core.pos+1, bam_calend(&(hit->core), bam1_cigar(hit)), cnt->offset, bam1_qname(hit));
    for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
      if(cnt->op[i]) {                    //  target base is 'C'
	  //char Twobit2base[] = {'X', 'A', 'C', 'X', 'G', 'X', 'X', 'X', 'T', 'X', 'X', 'X', 'X', 'X', 'X', 'N'};
	  //Rprintf("  adding to genomic position %d (read pos %d has %c)\n", i+cnt->offset+1, j+1, Twobit2base[bam1_seqi(hitseq, j)]);
	if(bam1_seqi(hitseq, j)==2) {        //  query base is 'C'
	  cnt->Tp[i]++;
	  cnt->Mp[i]++;
	} else if(bam1_seqi(hitseq, j)==8) { //  query base is 'T'
	  cnt->Tp[i]++;
	}
      }
  }

  return 0;
}


static int addHitToSNP(const bam1_t *hit, void *data) { // bam_fetch callback function of detect_SNVs()
    // REMARKS:
    //   assume ungapped read-global alignment
    //   hit->core.pos and count vectors are zero-based (add one during output)
    static uint8_t *hitseq=NULL;
    static unsigned long int i=0;
    static unsigned long int iend=0;
    static int j=0;
    static snpCounters *cnt = NULL;

    hitseq = bam1_seq(hit);
    cnt = (snpCounters*) data;
    iend = bam_calend(&(hit->core), bam1_cigar(hit)) - cnt->offset;

    if ((hit->core.flag & BAM_FPROPER_PAIR) && (hit->core.isize > 0) && (iend > hit->core.mpos - cnt->offset))
	// left fragment of a paired alignment --> make sure iend does not overlap alignment of right fragment
	iend = hit->core.mpos - cnt->offset;

    if (hit->core.flag & BAM_FREVERSE) {       // alignment on minus strand (reads are reverse complemented, look for C-C matches on opposite strand)
	for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
	    if(cnt->targetC[i]) {                 //  target base is 'C'
		if(bam1_seqi(hitseq, j)==2) {        //  query base is 'C'
		    cnt->total[i]++;
		    cnt->match[i]++;
		} else {                             //  query base is not 'C'
		    cnt->total[i]++;
		}
	    }

    } else {                                   // alignment on plus strand (look for G-G matches on opposite strand)
	for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
	    if(cnt->targetG[i]) {                 // target base is 'G'
		if(bam1_seqi(hitseq, j)==4) {        //  query base is 'G'
		    cnt->total[i]++;
		    cnt->match[i]++;
		} else {                             //  query base is not 'G'
		    cnt->total[i]++;
		}
	    }
    }

    return 0;
}


static int addHitToCountsAllele(const bam1_t *hit, void *data) { // bam_fetch callback function of quantify_methylation_allele()
  // REMARKS:
  //   assume ungapped read-global alignment
  //   count events separately for +/- strands (collapse strands during output if necessary) and for allele flag (R/U/A)
  //   hit->core.pos and count vectors are zero-based (add one during output)
  static uint8_t *hitseq=NULL;
  static unsigned long int i=0;
  static unsigned long int iend=0;
  static int j=0, a=0;
  static methCountersAllele *cnt = NULL;

  hitseq = bam1_seq(hit);
  cnt = (methCountersAllele*) data;
  iend = bam_calend(&(hit->core), bam1_cigar(hit)) - cnt->offset;
  a = alleleFlagToInt((char)*(uint8_t*)(bam_aux_get(hit,"XV") + 1));

  if ((hit->core.flag & BAM_FPROPER_PAIR) && (hit->core.isize > 0) && (iend > hit->core.mpos - cnt->offset))
      // left fragment of a paired alignment --> make sure iend does not overlap alignment of right fragment
      iend = hit->core.mpos - cnt->offset;

  if (hit->core.flag & BAM_FREVERSE) {       // alignment on minus strand (reads are reverse complemented, look for G-A mismatches)
    for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
      if(cnt->om[i]) {                   //  target base is 'G'
	if(bam1_seqi(hitseq, j)==4) {        //  query base is 'G'
	  cnt->Tm[a][i]++;
	  cnt->Mm[a][i]++;
	} else if(bam1_seqi(hitseq, j)==1) { //  query base is 'A'
	  cnt->Tm[a][i]++;
	}
      }

  } else {                                   // alignment on plus strand (look for C-T mismatches)
    for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
      if(cnt->op[i]) {                    //  target base is 'C'
	if(bam1_seqi(hitseq, j)==2) {        //  query base is 'C'
	  cnt->Tp[a][i]++;
	  cnt->Mp[a][i]++;
	} else if(bam1_seqi(hitseq, j)==8) { //  query base is 'T'
	  cnt->Tp[a][i]++;
	}
      }
  }

  return 0;
}


static int addHitToCountsSingleAlignments(const bam1_t *hit, void *data) { // bam_fetch callback function of quantify_methylation_singleAlignments()
  // REMARKS:
  //   assume ungapped read-global alignment
  //   count events separately for +/- strands (collapse strands during output if necessary)
  //   hit->core.pos and count vectors are zero-based (add one during output)
  //   bam_calend() return zero-based, exclusive end
  static uint8_t *hitseq=NULL;
  static unsigned long int i=0;
  static unsigned long int iend=0;
  static int j=0;
  static methCountersSingleAlignments *cnt = NULL;

  hitseq = bam1_seq(hit);
  cnt = (methCountersSingleAlignments*) data;
  iend = bam_calend(&(hit->core), bam1_cigar(hit)) - cnt->offset;

  if ((hit->core.flag & BAM_FPROPER_PAIR) && (hit->core.isize > 0) && (iend > hit->core.mpos - cnt->offset))
      // left fragment of a paired alignment --> make sure iend does not overlap alignment of right fragment
      iend = hit->core.mpos - cnt->offset;

  // scan alignment
  if (hit->core.flag & BAM_FREVERSE) {       // alignment on minus strand (reads are reverse complemented, look for G-A mismatches)
      for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
	  if(cnt->om[i]) {                   //  target base is 'G'
	      if(bam1_seqi(hitseq, j)==4) {        //  query base is 'G' (methylated)
		  cnt->aid.push_back( string(bam1_qname(hit)) );
		  cnt->Cid.push_back( i );
		  cnt->strand.push_back( '-' );
		  cnt->meth.push_back( 1 );
	      } else if(bam1_seqi(hitseq, j)==1) { //  query base is 'A' (unmethylated)
		  cnt->aid.push_back( string(bam1_qname(hit)) );
		  cnt->Cid.push_back( i );
		  cnt->strand.push_back( '-' );
		  cnt->meth.push_back( 0 );
	      }
	  }

  } else {                                   // alignment on plus strand (look for C-T mismatches)
      for(i=hit->core.pos-cnt->offset, j=0; i<iend; i++, j++)
	  if(cnt->op[i]) {                    //  target base is 'C'
	      if(bam1_seqi(hitseq, j)==2) {        //  query base is 'C' (methylated)
		  cnt->aid.push_back( string(bam1_qname(hit)) );
		  cnt->Cid.push_back( i );
		  cnt->strand.push_back( '+' );
		  cnt->meth.push_back( 1 );
	      } else if(bam1_seqi(hitseq, j)==8) { //  query base is 'T' (unmethylated)
		  cnt->aid.push_back( string(bam1_qname(hit)) );
		  cnt->Cid.push_back( i );
		  cnt->strand.push_back( '+' );
		  cnt->meth.push_back( 0 );
	      }
	  }
  }

  return 0;
}


/*! @function
  @abstract  verify the parameters of the quantify_methylation, detect_SNPs, quantify_methylation_allele functions
  @param  bamfile        Name of the bamfile

  @return       0 if successful
 */
int _verify_parameters(SEXP infiles, SEXP regionChr, SEXP regionChrLen, SEXP regionStart,
		       SEXP seqstring, SEXP mode, SEXP returnZero){

    if (!Rf_isString(infiles))
	Rf_error("'infiles' must be a character vector");
    if (!Rf_isString(regionChr) || 1 != Rf_length(regionChr))
	Rf_error("'regionChr' must be a single character value");
    if (!Rf_isInteger(regionChrLen) || 1 != Rf_length(regionChrLen))
	Rf_error("'regionChrLen' must be integer(1)");
    if (!Rf_isInteger(regionStart) || 1 != Rf_length(regionStart))
        Rf_error("'regionStart' must be integer(1)");
    if (!Rf_isString(seqstring) || 1 != Rf_length(seqstring))
	Rf_error("'seqstring' must be a single character value");
    if (mode!=NULL && (!Rf_isInteger(mode) || 1 != Rf_length(mode)))
        Rf_error("'mode' must be integer(1)");
    if (returnZero!=NULL && (!Rf_isLogical(returnZero) || 1 != Rf_length(returnZero)))
        Rf_error("'returnZero' must be logical(1)");

    return 0;
}


/*!
  @function  _scanSeqForCG
  @abstract  scan dna sequence ('seq') for presence of C's and G's and initialize output arrays (outputPlus, ouputMinus)
             according to the analysis mode (mode_int)
  @param  mode          analysis mode:
                          mode == 0 : only C's in CpG context (+/- strands collapsed)
                                  1 : only C's in CpG context (+/- strands separate)
                                  2 : all C's (+/- strands separate)
                                  3 : SNP detection, only C's in CpG context (+/- strands separate) --> should never come here (see detect_SNPs)
  @param  seq           character array with target DNA sequence
  @param  seqlen        length of seq
  @param  leftextension offset of output arrays to account for alignments that are partially outside of region of interest
  @param  outputPlus    pointer to bool array with true at positions to be quantified on the plus strand
  @param  outputMinus   pointer to bool array with true at positions to be quantified on the minus strand
  @param  pnOutputPlus  pointer to int (number of true values in outputPlus)
  @param  pnOutputMinus pointer to int (number of true values in outputMinus)
  @param  pnOutput      pointer to int (number of true values in total)

  @return       0 if successful
 */
int _scanSeqForCG(int mode, const char* seq, int seqlen, int leftextension, bool* outputPlus, bool* outputMinus,
		  int* pnOutputPlus, int* pnOutputMinus, int* pnOutput)
{
    int i=0, nOutputPlus=0, nOutputMinus=0, nOutput=0;

    if(mode == 2) {
	// all C's (+/- strands separate)
	for(i=0; i<seqlen; i++) {
	    if(seq[i]=='C' || seq[i]=='c') {
		outputPlus[i+leftextension] = true;
		nOutputPlus++;
	    } else if(seq[i]=='G' || seq[i]=='g') {
		outputMinus[i+leftextension] = true;
		nOutputMinus++;
	    }
	}
	nOutput = nOutputPlus + nOutputMinus;
    
    } else if((mode == 1) || (mode == 0)) {
	// only C's in CpG context (+/- strands separate) OR
	// only C's in CpG context (+/- strands collapsed)
	for(i=0; i<seqlen-1; i++)
	    if((seq[i]=='C' || seq[i]=='c') && (seq[i+1]=='G' || seq[i+1]=='g')) {
		outputPlus[i+leftextension] = true;
		outputMinus[i+1+leftextension] = true;
		nOutputPlus++;
	    }
	nOutput = (mode == 1) ? 2*nOutputPlus : nOutputPlus;

    } else {
	Rf_error("unknown mode '%d', should be one of 0, 1, or 2.\n", mode);
	return 1;
    }

    (*pnOutputPlus) = nOutputPlus;
    (*pnOutputMinus) = nOutputMinus;
    (*pnOutput) = nOutput;

    return 0;
}


/*!
  @function  quantify_methylation
  @abstract  parse bis-seq alignments and quantify methylation states
  @param  infiles        character vector with one or several bam file names (counts will be summed)
  @param  regionChr      character(1) with target sequence (chromosome) name
  @param  regionChrLen   integer(1) with the length of the target sequence (chromosome)
  @param  regionStart    integer(1) with position on target sequence (chromosome) to start quantification of methylation states
                         ('regionEnd' is defined by the regionStart + strlen(seqstring))
  @param  seqstring      character(1) with the reference sequence [regionStart, ...] on regionChr
  @param  mode           analysis mode:
                             mode == 0 : only C's in CpG context (+/- strands collapsed)
                                     1 : only C's in CpG context (+/- strands separate)
                                     2 : all C's (+/- strands separate)
                                     3 : SNP detection, only C's in CpG context (+/- strands separate) --> should never come here (see detect_SNPs)
  @param  returnZero     if true, keep C's with zero counts in the return value

  @return list containing five vectors (one element for each C or CpG) with chr, position, strand, total and methylated counts
 */
SEXP quantify_methylation(SEXP infiles, SEXP regionChr, SEXP regionChrLen, SEXP regionStart,
			  SEXP seqstring, SEXP mode, SEXP returnZero) {
    // validate arguments
    _verify_parameters(infiles, regionChr, regionChrLen, regionStart, seqstring, mode, returnZero);

    // declare parameters
    SEXP regionChrFirst = STRING_ELT(regionChr, 0), strandPlus = Rf_mkChar("+"), strandMinus = Rf_mkChar("-"), strandAny = Rf_mkChar("*");
    const char *target_name = Rf_translateChar(regionChrFirst);
    const char *seq = Rf_translateChar(STRING_ELT(seqstring, 0));
    int i = 0, j = 0, tid = 0, mode_int = Rf_asInteger(mode), nbIn = Rf_length(infiles),
	start = Rf_asInteger(regionStart) - 1, seqlen = strlen(seq), end = 0;
    end = start + seqlen; // end: 0-based, exclusive
    bool keepZero = Rf_asLogical(returnZero);

    const char **inf = (const char**) R_Calloc(nbIn, char*);
    for(i=0; i<nbIn; i++)
	inf[i] = Rf_translateChar(STRING_ELT(infiles, i));

    int *cntPlusT = NULL, *cntPlusM = NULL, *cntMinusT = NULL, *cntMinusM = NULL;
    cntPlusT = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusM = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusT = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusM = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    bool *outputPlus = NULL, *outputMinus = NULL;
    outputPlus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    outputMinus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    int nOutputPlus = 0, nOutputMinus = 0, nOutput = 0;
    int leftextension = (start < MAX_READ_LENGTH ? start : MAX_READ_LENGTH);
    methCounters data;
    data.Tp = cntPlusT;
    data.Mp = cntPlusM;
    data.Tm = cntMinusT;
    data.Mm = cntMinusM;
    data.op = outputPlus;
    data.om = outputMinus;
    data.offset = start - leftextension;

    // scan seq and initialize counters
    _scanSeqForCG(mode_int, seq, seqlen, leftextension, outputPlus, outputMinus, &nOutputPlus, &nOutputMinus, &nOutput);

    // loop over infiles and add to counters
    bam1_t *hit = bam_init1();
    samfile_t *fin;
    bam_index_t *idx;

    for(i=0; i<nbIn; i++) {
	fin = _bam_tryopen(inf[i], "rb", NULL);
	idx = bam_index_load(inf[i]); // load BAM index
	if (idx == 0)
	    Rf_error("BAM index for '%s' unavailable\n", inf[i]);


	// get target id
	tid = 0;
	while(strcmp(fin->header->target_name[tid], target_name) && tid+1<fin->header->n_targets)
	    tid++;

	if(strcmp(fin->header->target_name[tid], target_name))
	    Rf_error("could not find target '%s' in bam header of '%s'.\n", target_name, inf[i]);


	// call addHitToCounts on all alignments in region
	bam_fetch(fin->x.bam, idx, tid, start, end, &data, &addHitToCounts);


	// clean bam file objects
	bam_index_destroy(idx);
	samclose(fin);
    }

    bam_destroy1(hit);


    // allocate result objects
    if(!keepZero) {
	// re-count the number of output elements (<= nOutput)
	j = 0;
	
	if((mode_int == 2) || (mode_int == 1)) {
	    // all C's (+/- strands separate) OR
	    // only C's in CpG context (+/- strands separate)
	    for(i=leftextension; i<seqlen+leftextension; i++) {
		if((outputPlus[i] && cntPlusT[i]>0) || (outputMinus[i] && cntMinusT[i]>0))
		    j++;
	    }
	} else if(mode_int == 0) {
	    // only C's in CpG context (+/- strands collapsed)
	    for(i=leftextension; i<seqlen-1+leftextension; i++) {
		if(outputPlus[i] && (cntPlusT[i]>0 || cntMinusT[i+1]>0))
		    j++;
	    }
	}
	nOutput = j;
    }

    SEXP resChr, resPos, resStrand, resT, resM, res, resNames;
    PROTECT(resChr = Rf_allocVector(STRSXP, nOutput));
    PROTECT(resPos = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resStrand = Rf_allocVector(STRSXP, nOutput));
    PROTECT(resT = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resM = Rf_allocVector(INTSXP, nOutput));
    PROTECT(res = Rf_allocVector(VECSXP, 5));
    PROTECT(resNames = Rf_allocVector(STRSXP, 5));

    int *resPos_int = INTEGER(resPos), *resT_int = INTEGER(resT), *resM_int = INTEGER(resM);


    // fill in results objects
    j = 0;

    if((mode_int == 2) || (mode_int == 1)) {
	// all C's (+/- strands separate) OR
	// only C's in CpG context (+/- strands separate)

	for(i=leftextension; i<seqlen+leftextension; i++) {
	    if(outputPlus[i] && (keepZero || cntPlusT[i]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandPlus);
		resT_int[j] = cntPlusT[i];
		resM_int[j] = cntPlusM[i];
		j++;
	    }
	    if(outputMinus[i] && (keepZero || cntMinusT[i]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandMinus);
		resT_int[j] = cntMinusT[i];
		resM_int[j] = cntMinusM[i];
		j++;
	    }
	}

    } else if(mode_int == 0) {
	// only C's in CpG context (+/- strands collapsed)
	for(i=leftextension; i<seqlen-1+leftextension; i++) {
	    if(outputPlus[i] && (keepZero || cntPlusT[i]>0 || cntMinusT[i+1]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandAny);
		resT_int[j] = cntPlusT[i] + cntMinusT[i+1];
		resM_int[j] = cntPlusM[i] + cntMinusM[i+1];
		j++;
	    }
	}
    }

    SET_VECTOR_ELT(res, 0, resChr);
    SET_VECTOR_ELT(res, 1, resPos);
    SET_VECTOR_ELT(res, 2, resStrand);
    SET_VECTOR_ELT(res, 3, resT);
    SET_VECTOR_ELT(res, 4, resM);

    SET_STRING_ELT(resNames, 0, Rf_mkChar("chr"));
    SET_STRING_ELT(resNames, 1, Rf_mkChar("position"));
    SET_STRING_ELT(resNames, 2, Rf_mkChar("strand"));
    SET_STRING_ELT(resNames, 3, Rf_mkChar("T"));
    SET_STRING_ELT(resNames, 4, Rf_mkChar("M"));
    Rf_setAttrib(res, R_NamesSymbol, resNames);

    // clean up
    R_Free(inf);
    R_Free(cntPlusT);
    R_Free(cntPlusM);
    R_Free(cntMinusT);
    R_Free(cntMinusM);
    R_Free(outputPlus);
    R_Free(outputMinus);

    // return
    UNPROTECT(7);
    return(res);
}

SEXP detect_SNVs(SEXP infiles, SEXP regionChr, SEXP regionChrLen, SEXP regionStart,
		 SEXP seqstring, SEXP returnZero) {

    // validate arguments
    _verify_parameters(infiles, regionChr, regionChrLen, regionStart, seqstring, NULL, returnZero);

    // declare parameters
    SEXP regionChrFirst = STRING_ELT(regionChr, 0);
    const char *target_name = Rf_translateChar(regionChrFirst);
    const char *seq = Rf_translateChar(STRING_ELT(seqstring, 0));
    int i = 0, j = 0, tid = 0, nbIn = Rf_length(infiles),
	start = Rf_asInteger(regionStart) - 1, seqlen = strlen(seq), end = 0;
    end = start + seqlen; // end: 0-based, exclusive
    bool keepZero = Rf_asLogical(returnZero);

    const char **inf = (const char**) R_Calloc(nbIn, char*);
    for(i=0; i<nbIn; i++)
	inf[i] = Rf_translateChar(STRING_ELT(infiles, i));

    unsigned int *cntMatch = NULL, *cntTotal = NULL;
    cntMatch = (unsigned int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, unsigned int);
    cntTotal = (unsigned int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, unsigned int);
    bool *targetC = NULL, *targetG = NULL;
    targetC = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    targetG = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    int nTarget = 0;
    int leftextension = (start < MAX_READ_LENGTH ? start : MAX_READ_LENGTH);
    snpCounters data;
    data.match   = cntMatch;
    data.total   = cntTotal;
    data.targetC = targetC;
    data.targetG = targetG;
    data.offset  = start - leftextension;

    // only C's or G's in CpG context
    for(i=0; i<seqlen-1; i++)
	if((seq[i]=='C' || seq[i]=='c') && (seq[i+1]=='G' || seq[i+1]=='g')) {
	    targetC[i] = true;
	    targetG[i+1] = true;
	    nTarget += 2;
	}

    // loop over infiles and add to counters
    bam1_t *hit = bam_init1();
    samfile_t *fin;
    bam_index_t *idx;

    for(i=0; i<nbIn; i++) {
	fin = _bam_tryopen(inf[i], "rb", NULL);
	idx = bam_index_load(inf[i]); // load BAM index
	if (idx == 0)
	    Rf_error("BAM index for '%s' unavailable\n", inf[i]);


	// get target id
	tid = 0;
	while(strcmp(fin->header->target_name[tid], target_name) && tid+1<fin->header->n_targets)
	    tid++;

	if(strcmp(fin->header->target_name[tid], target_name))
	    Rf_error("could not find target '%s' in bam header of '%s'.\n", target_name, inf[i]);

	// call addHitToCounts on all alignments in region
	bam_fetch(fin->x.bam, idx, tid, start, end, &data, &addHitToSNP);


	// clean bam file objects
	bam_index_destroy(idx);
	samclose(fin);
    }

    bam_destroy1(hit);


    // allocate result objects
    if(!keepZero) {
	// re-count the number of output elements (<= nTarget)
	j = 0;
	
	for(i=leftextension; i<seqlen+leftextension; i++) {
	    if((targetC[i] || targetG[i]) && cntTotal[i]>0)
	    j++;
	}
	nTarget = j;
    }

    SEXP resChr, resPos, resMatch, resTotal, res, resNames;
    PROTECT(resChr = Rf_allocVector(STRSXP, nTarget));
    PROTECT(resPos = Rf_allocVector(INTSXP, nTarget));
    PROTECT(resMatch = Rf_allocVector(INTSXP, nTarget));
    PROTECT(resTotal = Rf_allocVector(INTSXP, nTarget));
    PROTECT(res = Rf_allocVector(VECSXP, 4));
    PROTECT(resNames = Rf_allocVector(STRSXP, 4));

    int *resPos_int = INTEGER(resPos), *resMatch_int = INTEGER(resMatch), *resTotal_int = INTEGER(resTotal);


    // fill in results objects
    j = 0;
    for(i=leftextension; i<seqlen+leftextension; i++) {
	if((targetC[i] || targetG[i]) && (keepZero || cntTotal[i]>0)) {
	    SET_STRING_ELT(resChr, j, regionChrFirst);
	    resPos_int[j] = i + data.offset + 1;
	    resMatch_int[j] = cntMatch[i];
	    resTotal_int[j] = cntTotal[i];
	    j++;
	}
    }

    SET_VECTOR_ELT(res, 0, resChr);
    SET_VECTOR_ELT(res, 1, resPos);
    SET_VECTOR_ELT(res, 2, resTotal);
    SET_VECTOR_ELT(res, 3, resMatch);

    SET_STRING_ELT(resNames, 0, Rf_mkChar("chr"));
    SET_STRING_ELT(resNames, 1, Rf_mkChar("position"));
    SET_STRING_ELT(resNames, 2, Rf_mkChar("nTotal"));
    SET_STRING_ELT(resNames, 3, Rf_mkChar("nMatch"));
    Rf_setAttrib(res, R_NamesSymbol, resNames);

    // clean up
    R_Free(inf);
    R_Free(cntMatch);
    R_Free(cntTotal);
    R_Free(targetC);
    R_Free(targetG);

    // return
    UNPROTECT(6);
    return(res);
}

SEXP quantify_methylation_allele(SEXP infiles, SEXP regionChr, SEXP regionChrLen, SEXP regionStart,
				 SEXP seqstring, SEXP mode, SEXP returnZero) {
    /*
      mode == 0 : only C's in CpG context (+/- strands collapsed)
              1 : only C's in CpG context (+/- strands separate)
	      2 : all C's (+/- strands separate)
    */

    // validate arguments
    _verify_parameters(infiles, regionChr, regionChrLen, regionStart, seqstring, mode, returnZero);

    // declare parameters
    SEXP regionChrFirst = STRING_ELT(regionChr, 0), strandPlus = Rf_mkChar("+"), strandMinus = Rf_mkChar("-"), strandAny = Rf_mkChar("*");
    const char *target_name = Rf_translateChar(regionChrFirst);
    const char *seq = Rf_translateChar(STRING_ELT(seqstring, 0));
    int i = 0, j = 0, tid = 0, mode_int = Rf_asInteger(mode), nbIn = Rf_length(infiles),
	start = Rf_asInteger(regionStart) - 1, seqlen = strlen(seq), end = 0;
    end = start + seqlen; // end: 0-based, exclusive
    bool keepZero = Rf_asLogical(returnZero);

    const char **inf = (const char**) R_Calloc(nbIn, char*);
    for(i=0; i<nbIn; i++)
	inf[i] = Rf_translateChar(STRING_ELT(infiles, i));

    int *cntPlusTR = NULL, *cntPlusTU = NULL, *cntPlusTA = NULL,
	*cntPlusMR = NULL, *cntPlusMU = NULL, *cntPlusMA = NULL,
	*cntMinusTR = NULL, *cntMinusTU = NULL, *cntMinusTA = NULL,
	*cntMinusMR = NULL, *cntMinusMU = NULL, *cntMinusMA = NULL;
    cntPlusTR = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusTU = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusTA = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusMR = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusMU = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntPlusMA = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusTR = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusTU = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusTA = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusMR = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusMU = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    cntMinusMA = (int*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, int);
    bool *outputPlus = NULL, *outputMinus = NULL;
    outputPlus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    outputMinus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    int nOutputPlus = 0, nOutputMinus = 0, nOutput = 0;
    int leftextension = (start < MAX_READ_LENGTH ? start : MAX_READ_LENGTH);
    methCountersAllele data;
    data.Tp[0] = cntPlusTR;
    data.Tp[1] = cntPlusTU;
    data.Tp[2] = cntPlusTA;
    data.Mp[0] = cntPlusMR;
    data.Mp[1] = cntPlusMU;
    data.Mp[2] = cntPlusMA;
    data.Tm[0] = cntMinusTR;
    data.Tm[1] = cntMinusTU;
    data.Tm[2] = cntMinusTA;
    data.Mm[0] = cntMinusMR;
    data.Mm[1] = cntMinusMU;
    data.Mm[2] = cntMinusMA;
    data.op = outputPlus;
    data.om = outputMinus;
    data.offset = start - leftextension;

    // scan seq and initialize counters
    _scanSeqForCG(mode_int, seq, seqlen, leftextension, outputPlus, outputMinus, &nOutputPlus, &nOutputMinus, &nOutput);

    // loop over infiles and add to counters
    bam1_t *hit = bam_init1();
    samfile_t *fin;
    bam_index_t *idx;

    for(i=0; i<nbIn; i++) {
	fin = _bam_tryopen(inf[i], "rb", NULL);
	idx = bam_index_load(inf[i]); // load BAM index
	if (idx == 0)
	    Rf_error("BAM index for '%s' unavailable\n", inf[i]);


	// get target id
	tid = 0;
	while(strcmp(fin->header->target_name[tid], target_name) && tid+1<fin->header->n_targets)
	    tid++;

	if(strcmp(fin->header->target_name[tid], target_name))
	    Rf_error("could not find target '%s' in bam header of '%s'.\n", target_name, inf[i]);


	// call addHitToCounts on all alignments in region
	bam_fetch(fin->x.bam, idx, tid, start, end, &data, &addHitToCountsAllele);


	// clean bam file objects
	bam_index_destroy(idx);
	samclose(fin);
    }

    bam_destroy1(hit);


    // allocate result objects
    if(!keepZero) {
	// re-count the number of output elements (<= nOutput)
	j = 0;
	
	if((mode_int == 2) || (mode_int == 1)) {
	    // all C's (+/- strands separate) OR
	    // only C's in CpG context (+/- strands separate)
	    for(i=leftextension; i<seqlen+leftextension; i++) {
		if((outputPlus[i] && (cntPlusTR[i]>0 || cntPlusTU[i]>0 || cntPlusTA[i]>0)) ||
		   (outputMinus[i] && (cntMinusTR[i]>0 || cntMinusTU[i]>0 || cntMinusTA[i]>0)))
		    j++;
	    }
	} else if(mode_int == 0) {
	    // only C's in CpG context (+/- strands collapsed)
	    for(i=leftextension; i<seqlen-1+leftextension; i++) {
		if(outputPlus[i] && (cntPlusTR[i]>0 || cntPlusTU[i]>0 || cntPlusTA[i]>0 || cntMinusTR[i+1]>0 || cntMinusTU[i+1]>0 || cntMinusTA[i+1]>0))
		    j++;
	    }
	}
	nOutput = j;
    }

    SEXP resChr, resPos, resStrand, resTR, resTU, resTA, resMR, resMU, resMA, res, resNames;
    PROTECT(resChr = Rf_allocVector(STRSXP, nOutput));
    PROTECT(resPos = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resStrand = Rf_allocVector(STRSXP, nOutput));
    PROTECT(resTR = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resTU = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resTA = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resMR = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resMU = Rf_allocVector(INTSXP, nOutput));
    PROTECT(resMA = Rf_allocVector(INTSXP, nOutput));
    PROTECT(res = Rf_allocVector(VECSXP, 9));
    PROTECT(resNames = Rf_allocVector(STRSXP, 9));

    int *resPos_int = INTEGER(resPos),
	*resTR_int = INTEGER(resTR), *resTU_int = INTEGER(resTU), *resTA_int = INTEGER(resTA),
	*resMR_int = INTEGER(resMR), *resMU_int = INTEGER(resMU), *resMA_int = INTEGER(resMA);


    // fill in results objects
    j = 0;

    if((mode_int == 2) || (mode_int == 1)) {
	// all C's (+/- strands separate) OR
	// only C's in CpG context (+/- strands separate)

	for(i=leftextension; i<seqlen+leftextension; i++) {
	    if(outputPlus[i] && (keepZero || cntPlusTR[i]>0 || cntPlusTU[i]>0 || cntPlusTA[i]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandPlus);
		resTR_int[j] = cntPlusTR[i];
		resTU_int[j] = cntPlusTU[i];
		resTA_int[j] = cntPlusTA[i];
		resMR_int[j] = cntPlusMR[i];
		resMU_int[j] = cntPlusMU[i];
		resMA_int[j] = cntPlusMA[i];
		j++;
	    }
	    if(outputMinus[i] && (keepZero || cntMinusTR[i]>0 || cntMinusTU[i]>0 || cntMinusTA[i]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandMinus);
		resTR_int[j] = cntMinusTR[i];
		resTU_int[j] = cntMinusTU[i];
		resTA_int[j] = cntMinusTA[i];
		resMR_int[j] = cntMinusMR[i];
		resMU_int[j] = cntMinusMU[i];
		resMA_int[j] = cntMinusMA[i];
		j++;
	    }
	}

    } else if(mode_int == 0) {
	// only C's in CpG context (+/- strands collapsed)
	for(i=leftextension; i<seqlen-1+leftextension; i++) {
	    if(outputPlus[i] && (keepZero ||
				 cntPlusTR[i]>0 || cntPlusTU[i]>0 || cntPlusTA[i]>0 ||
				 cntMinusTR[i+1]>0 || cntMinusTU[i+1]>0 || cntMinusTA[i+1]>0)) {
		SET_STRING_ELT(resChr, j, regionChrFirst);
		resPos_int[j] = i + data.offset + 1;
		SET_STRING_ELT(resStrand, j, strandAny);
		resTR_int[j] = cntPlusTR[i] + cntMinusTR[i+1];
		resTU_int[j] = cntPlusTU[i] + cntMinusTU[i+1];
		resTA_int[j] = cntPlusTA[i] + cntMinusTA[i+1];
		resMR_int[j] = cntPlusMR[i] + cntMinusMR[i+1];
		resMU_int[j] = cntPlusMU[i] + cntMinusMU[i+1];
		resMA_int[j] = cntPlusMA[i] + cntMinusMA[i+1];
		j++;
	    }
	}
    }

    SET_VECTOR_ELT(res, 0, resChr);
    SET_VECTOR_ELT(res, 1, resPos);
    SET_VECTOR_ELT(res, 2, resStrand);
    SET_VECTOR_ELT(res, 3, resTR);
    SET_VECTOR_ELT(res, 4, resMR);
    SET_VECTOR_ELT(res, 5, resTU);
    SET_VECTOR_ELT(res, 6, resMU);
    SET_VECTOR_ELT(res, 7, resTA);
    SET_VECTOR_ELT(res, 8, resMA);

    SET_STRING_ELT(resNames, 0, Rf_mkChar("chr"));
    SET_STRING_ELT(resNames, 1, Rf_mkChar("position"));
    SET_STRING_ELT(resNames, 2, Rf_mkChar("strand"));
    SET_STRING_ELT(resNames, 3, Rf_mkChar("TR"));
    SET_STRING_ELT(resNames, 4, Rf_mkChar("MR"));
    SET_STRING_ELT(resNames, 5, Rf_mkChar("TU"));
    SET_STRING_ELT(resNames, 6, Rf_mkChar("MU"));
    SET_STRING_ELT(resNames, 7, Rf_mkChar("TA"));
    SET_STRING_ELT(resNames, 8, Rf_mkChar("MA"));
    Rf_setAttrib(res, R_NamesSymbol, resNames);

    // clean up
    R_Free(inf);
    R_Free(cntPlusTR);
    R_Free(cntPlusTU);
    R_Free(cntPlusTA);
    R_Free(cntPlusMR);
    R_Free(cntPlusMU);
    R_Free(cntPlusMA);
    R_Free(cntMinusTR);
    R_Free(cntMinusTU);
    R_Free(cntMinusTA);
    R_Free(cntMinusMR);
    R_Free(cntMinusMU);
    R_Free(cntMinusMA);
    R_Free(outputPlus);
    R_Free(outputMinus);

    // return
    UNPROTECT(11);
    return(res);
}

/*!
  @function  quantify_methylation_singleAlignments
  @abstract  parse bis-seq alignments and quantify methylation states (report results for individual reads)
  @param  infiles        character vector with one or several bam file names (results will be pooled)
  @param  regionChr      character(1) with target sequence (chromosome) name
  @param  regionChrLen   integer(1) with the length of the target sequence (chromosome)
  @param  regionStart    integer(1) with position on target sequence (chromosome) to start quantification of methylation states
                         ('regionEnd' is defined by the regionStart + strlen(seqstring))
  @param  seqstring      character(1) with the reference sequence [regionStart, ...] on regionChr
  @param  mode           analysis mode:
                             mode == 0 : *NOT ALLOWED HERE* only C's in CpG context (+/- strands collapsed)
                                     1 : only C's in CpG context (+/- strands separate)
                                     2 : all C's (+/- strands separate)
                                     3 : *NOT ALLOWED HERE* SNP detection, only C's in CpG context (+/- strands separate)

  @return list containing elements:
           aid    : integer vector with unique alignment identifiers
           Cid    : integer vector with unique C identifiers (C's in seqstring)
           meth   : integer vector with 1 (methylated) or 0 (unmethylated)
           // Cinfo  : character vector with annotation information for unique Cid's
 */
SEXP quantify_methylation_singleAlignments(SEXP infiles, SEXP regionChr, SEXP regionChrLen, SEXP regionStart,
					   SEXP seqstring, SEXP mode) {
    // validate arguments
    _verify_parameters(infiles, regionChr, regionChrLen, regionStart, seqstring, mode, NULL);

    // declare parameters
    const char *target_name = Rf_translateChar(STRING_ELT(regionChr, 0));
    const char *seq = Rf_translateChar(STRING_ELT(seqstring, 0));
    int i = 0, j = 0, tid = 0, mode_int = Rf_asInteger(mode), nbIn = Rf_length(infiles),
	start = Rf_asInteger(regionStart) - 1, seqlen = strlen(seq), end = 0;
    end = start + seqlen; // end: 0-based, exclusive

    if(mode_int != 1 && mode_int != 2)
	Rf_error("'mode' (%d) must be 1 or 2 for quantify_methylation_singleAlignments.\n", mode_int);

    const char **inf = (const char**) R_Calloc(nbIn, char*);
    for(i=0; i<nbIn; i++)
	inf[i] = Rf_translateChar(STRING_ELT(infiles, i));

    bool *outputPlus = NULL, *outputMinus = NULL;
    outputPlus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    outputMinus = (bool*) R_Calloc(seqlen + 2*MAX_READ_LENGTH, bool);
    int nOutputPlus = 0, nOutputMinus = 0, nOutput = 0;
    int leftextension = (start < MAX_READ_LENGTH ? start : MAX_READ_LENGTH);
    methCountersSingleAlignments data;
    data.op = outputPlus;
    data.om = outputMinus;
    data.offset = start - leftextension;

    // scan seq and initialize counters
    _scanSeqForCG(mode_int, seq, seqlen, leftextension, outputPlus, outputMinus, &nOutputPlus, &nOutputMinus, &nOutput);

    // loop over infiles and add to counters
    bam1_t *hit = bam_init1();
    samfile_t *fin;
    bam_index_t *idx;

    for(i=0; i<nbIn; i++) {
	fin = _bam_tryopen(inf[i], "rb", NULL);
	idx = bam_index_load(inf[i]); // load BAM index
	if (idx == 0)
	    Rf_error("BAM index for '%s' unavailable\n", inf[i]);


	// get target id
	tid = 0;
	while(strcmp(fin->header->target_name[tid], target_name) && tid+1<fin->header->n_targets)
	    tid++;

	if(strcmp(fin->header->target_name[tid], target_name))
	    Rf_error("could not find target '%s' in bam header of '%s'.\n", target_name, inf[i]);


	// call addHitToCountsSingleAlignments on all alignments in region
	bam_fetch(fin->x.bam, idx, tid, start, end, &data, &addHitToCountsSingleAlignments);


	// clean bam file objects
	bam_index_destroy(idx);
	samclose(fin);
    }

    bam_destroy1(hit);


    // allocate result objects
    int resLength = data.aid.size();
    if(resLength != data.Cid.size() || resLength != data.meth.size())
	Rf_error("result structures are incompatible (%d/%d/%d); fatal error in quantify_methylation_singleAlignments\n",
		 resLength, data.Cid.size(), data.meth.size());
    SEXP resAid, resCid, resStrand, resMeth, res, resNames, strandPlus = Rf_mkChar("+"), strandMinus = Rf_mkChar("-");
    //PROTECT(resAid    = Rf_allocVector(INTSXP, resLength));
    PROTECT(resAid    = Rf_allocVector(STRSXP, resLength));
    PROTECT(resCid    = Rf_allocVector(INTSXP, resLength));
    PROTECT(resStrand = Rf_allocVector(STRSXP, resLength));
    PROTECT(resMeth   = Rf_allocVector(INTSXP, resLength));
    PROTECT(res       = Rf_allocVector(VECSXP, 4));
    PROTECT(resNames  = Rf_allocVector(STRSXP, 4));

    // fill in results objects
    //memcpy( INTEGER(resAid),  data.aid.data(),  sizeof(int)*resLength );
    memcpy( INTEGER(resCid),  data.Cid.data(),  sizeof(int)*resLength );
    memcpy( INTEGER(resMeth), data.meth.data(), sizeof(int)*resLength );
    for(i=0; i<resLength; i++) {
	SET_STRING_ELT(resAid,    i, Rf_mkChar(data.aid[i].c_str()));
	SET_STRING_ELT(resStrand, i, data.strand[i]=='-' ? strandMinus : strandPlus);
    }
    
    // transform C (DNA base) identifiers to genomic positions
    int *Cid = INTEGER(resCid);
    for(i=0; i<resLength; i++)
	Cid[i] += data.offset + 1;

    // combine results into list
    SET_VECTOR_ELT(res, 0, resAid);
    SET_VECTOR_ELT(res, 1, resCid);
    SET_VECTOR_ELT(res, 2, resStrand);
    SET_VECTOR_ELT(res, 3, resMeth);

    SET_STRING_ELT(resNames, 0, Rf_mkChar("aid"));
    SET_STRING_ELT(resNames, 1, Rf_mkChar("Cid"));
    SET_STRING_ELT(resNames, 2, Rf_mkChar("strand"));
    SET_STRING_ELT(resNames, 3, Rf_mkChar("meth"));
    Rf_setAttrib(res, R_NamesSymbol, resNames);

    // clean up
    R_Free(inf);
    R_Free(outputPlus);
    R_Free(outputMinus);

    // return
    UNPROTECT(6);
    return(res);
}

