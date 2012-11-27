/*!
  @header
  
  Counts the alignments in regions, which fulfill specific criteria.
  
  @author:    Anita Lerch
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "count_alignments.h"
#include <stdlib.h>

#define SMART_SHIFT -1000000 // for half insert size shift towards the mate read

/*! @typedef
  @abstract Structure to provid the data to the bam_fetch() functions.
  @field sumU         Current count of the fetch region (non-allelic or allelic Unknown)
  @field sumR         Current count of the fetch region (allelic Reference)
  @field sumA         Current count of the fetch region (allelic Alternative)
  @field start        Start of the fetch region
  @field end          End of the fetch region
  @field strand       Strand of the fetch regions
  @field shift        Shift size of the reads
  @field readBitMask  select first/last read from a multi-read experiment
  @field selectReadPosition  weight of alignment on "s"tart or "e"nd
  @field allelic      allelic true(1) or false(0)
*/
typedef struct {
    int sumU;
    int sumR;
    int sumA;
    int start;
    int end;
    const char * strand;
    int shift;
    int readBitMask;
    char selectReadPosition;
    int allelic;
} regionInfoSums;


/*! @function
  @abstract  increase allele specific counts.
  @param     hit     the alignment
  @param     rinfo   regionInfoSums
  @return    0 if successful
 */
static int _sum_allelic(const bam1_t *hit, regionInfoSums *rinfo){
    static uint8_t *xv_ptr = 0;

    // get XV tag
    xv_ptr = bam_aux_get(hit,"XV");
    if(xv_ptr == 0)
        Rf_error("XV tag missing but needed for allele-specific counting");

    // increase count
    switch(bam_aux2A(xv_ptr)){
    case 'U':
        rinfo->sumU += 1;
        break;
    case 'R':
        rinfo->sumR += 1;
        break;
    case 'A':
        rinfo->sumA += 1;
        break;
    default:
        Rf_error("'%c' is not a valid XV tag value; should be one of 'U','R' or 'A'", bam_aux2A(xv_ptr));
    }

    return 0;
}


/*! @function
  @abstract  callback for bam_fetch(); sums alignments if correct strand and shifted biological start/end position overlap fetch region
  @param  hit   the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _addValidHitToSums(const bam1_t *hit, void *data){
    regionInfoSums *rinfo = (regionInfoSums*)data;
    static double shift = 0;
    static int pos = 0;

    // skip alignment if read1 or read2 flag is set (=paired-end) and if wrong readBitMask
    if((hit->core.flag & (BAM_FREAD1 + BAM_FREAD2)) && (hit->core.flag & rinfo->readBitMask) == 0)
        return 0;

    // skip alignment if region is not * and the strand of alignment or region is not the same
    if(strcmp(rinfo->strand, "*")
       && (((hit->core.flag & BAM_FREVERSE) == 0) != (strcmp(rinfo->strand, "+") == 0)))
        return 0;

    // set shift
    if(rinfo->shift == SMART_SHIFT){
        // the sign of isize needs to be examined to make sure that shift(x) == -shift(-x)
        // example: 13 -> 12/2 but also -13 -> -12/2
        if(hit->core.isize > 0)
            shift = ((double)hit->core.isize-1)/2;
        else
            shift = ((double)hit->core.isize+1)/2;
    }
    else if((hit->core.flag & BAM_FREVERSE) == 0)
        // alignment on plus-strand
        shift = (double)rinfo->shift;
    else
        // alignment on minus-strand
        shift = -(double)rinfo->shift;

    // calculate position
    if(((hit->core.flag & BAM_FREVERSE) == 0) == (rinfo->selectReadPosition == 's')) // XOR
	// plus-strand and startwithin OR minus-strand and endwithin
        // --> position on the left side of the read 
	pos = (int)((double)hit->core.pos + shift); // 0-based inclusive start
    else
	// plus-strand and endwithin OR minus-strand and startwithin
        // --> position on the right side of the read 
        pos = (int)((double)bam_calend(&hit->core, bam1_cigar(hit)) - 1 + shift); // 0-based exclusive end --> -1

    // check if in region
    if(rinfo->start <= pos && pos < rinfo->end){
        if(rinfo->allelic == 0)
            rinfo->sumU += 1;
        else
            _sum_allelic(hit, rinfo);
    }

    return 0;
}


/*! @function
  @abstract  verify the parameters of the count_alignments_non_allelic and count_alignments_allelic function
  @param  bamfile        Name of the bamfile

  @return       0 if successful
 */
int _verify_parameters(SEXP bamfile, SEXP tid,  SEXP start, SEXP end, SEXP strand,
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden){
    // check bamfile parameter
    if(!Rf_isString(bamfile) || Rf_length(bamfile) != 1)
        Rf_error("'bamfile' must be of type character(1)");

    // check region parameter
    if(!Rf_isInteger(tid))
        Rf_error("'tid' must be of type integer");
    if(!Rf_isInteger(start))
        Rf_error("'start' must be of type integer");
    if(!Rf_isInteger(end))
        Rf_error("'end' must be of type integer");
    if(!Rf_isString(strand))
        Rf_error("'strand' must be of type character");
    int num_regions = Rf_length(tid);
    if(num_regions != Rf_length(start) || num_regions != Rf_length(end) || num_regions != Rf_length(strand))
        Rf_error("'tid', 'start', 'end', 'strand' must have equal length");

    // check parameter overlap type
    if(!Rf_isString(selectReadPosition) || Rf_length(selectReadPosition) != 1)
        Rf_error("'selectReadPosition' must be of type character(1)");
    if(Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0] != 's' 
       && Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0] != 'e')
        Rf_error("The value of 'selectReadPosition' not supportet.");

    // check parameter readBitMask, shift, broaden
    if(!Rf_isInteger(readBitMask) || Rf_length(readBitMask) != 1)
        Rf_error("'readBitMask' must be of type integer(1)");
    if(!Rf_isInteger(shift) && Rf_length(shift) != 1)
        Rf_error("'shift' must be of type integer(1)");
    if(!Rf_isInteger(broaden) && Rf_length(broaden) != 1)
        Rf_error("'broaden' must be of type integer(1)");
    if(INTEGER(broaden)[0] < 0)
        Rf_error("'broaden' must be a positive value.");

    return 0;
}

/*! @function
  @abstract  Counts the alignments in regions, which fit the strand, overlap type and shift criteria.
  @param  bamfile   Name of the bamfile
  @param  regions       Coordinates of fetch regions
  @param  selectReadPosition  Type of the overlap criterion
  @param  shift         Shift size of the reads
  @return               Vector of the alignment counts
 */
SEXP count_alignments_non_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand,
                                  SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden){

    // check parameters
    _verify_parameters(bamfile, tid, start, end, strand, selectReadPosition, readBitMask, shift, broaden);
    
    // open bam file
    samfile_t *fin = 0;
    fin = samopen(Rf_translateChar(STRING_ELT(bamfile, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open BAM file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("BAM header missing or empty of file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (idx == 0){
        samclose(fin);
        Rf_error("failed to open BAM index file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }

    // initialise regionInfoSums
    regionInfoSums rinfo;
    rinfo.readBitMask = INTEGER(readBitMask)[0];
    rinfo.shift = INTEGER(shift)[0];
    rinfo.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rinfo.allelic = 0;
    
    // set shift for fetch to zero if smart shift
    int shift_f = abs(INTEGER(shift)[0]);
    if(INTEGER(shift)[0] == SMART_SHIFT)
        shift_f = 0;
    
    // select bam_fetch callback function
    bam_fetch_f fetch_func = _addValidHitToSums;
    
    // loop over query regions
    int num_regions = Rf_length(tid);
    SEXP count;
    PROTECT(count = allocVector(INTSXP, num_regions));
    for(int i = 0; i < num_regions; i++){
        // reset counter and setup region in rinfo
        rinfo.sumU = 0;
        rinfo.start = INTEGER(start)[i];
        rinfo.end = INTEGER(end)[i];
        rinfo.strand = Rf_translateChar(STRING_ELT(strand, i));
        // process alignments that overlap region
        bam_fetch(fin->x.bam, idx, INTEGER(tid)[i], 
                  INTEGER(start)[i]-shift_f-INTEGER(broaden)[0], // 0-based inclusive start
                  INTEGER(end)[i]+shift_f+INTEGER(broaden)[0],   // 0-based exclusive end
                  &rinfo, fetch_func);
	// copy result to return values
        INTEGER(count)[i] = rinfo.sumU;
    }
    
    // clean up
    samclose(fin);
    bam_index_destroy(idx);
    
    UNPROTECT(1);
    
    return count;
}

SEXP count_alignments_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand,
                                  SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden){

    // check parameters
    _verify_parameters(bamfile, tid, start, end, strand, selectReadPosition, readBitMask, shift, broaden);
    
    // open bam file
    samfile_t *fin = 0;
    fin = samopen(Rf_translateChar(STRING_ELT(bamfile, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open BAM file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("BAM header missing or empty of file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (idx == 0){
        samclose(fin);
        Rf_error("failed to open BAM index file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }

    // initialise regionInfoSums
    regionInfoSums rinfo;
    rinfo.readBitMask = INTEGER(readBitMask)[0];
    rinfo.shift = INTEGER(shift)[0];
    rinfo.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rinfo.allelic = 1;

    // set shift for fetch to zero if smart shift
    int shift_f = abs(INTEGER(shift)[0]);
    if(INTEGER(shift)[0] == SMART_SHIFT)
        shift_f = 0;
    
    // select bam_fetch callback function
    bam_fetch_f fetch_func = _addValidHitToSums;
    
    // loop over query regions
    int num_regions = Rf_length(tid);
    SEXP count, attrib, countU, countR, countA;
    PROTECT(countU = allocVector(INTSXP, num_regions));
    PROTECT(countR = allocVector(INTSXP, num_regions));
    PROTECT(countA = allocVector(INTSXP, num_regions));
    for(int i = 0; i < num_regions; i++){
        // reset counters and setup region in rinfo
	rinfo.sumU = 0;
	rinfo.sumR = 0;
	rinfo.sumA = 0;
        rinfo.start = INTEGER(start)[i];
        rinfo.end = INTEGER(end)[i];
        rinfo.strand = Rf_translateChar(STRING_ELT(strand, i));
        // process alignments that overlap region
        bam_fetch(fin->x.bam, idx, INTEGER(tid)[i], 
                  INTEGER(start)[i]-shift_f-INTEGER(broaden)[0], // 0-based inclusive start
                  INTEGER(end)[i]+shift_f+INTEGER(broaden)[0],   // 0-based exclusive end
                  &rinfo, fetch_func);
	// copy result to return values
	INTEGER(countU)[i] = rinfo.sumU;
	INTEGER(countR)[i] = rinfo.sumR;
	INTEGER(countA)[i] = rinfo.sumA;
    }

    // pack results into list
    PROTECT(count = allocVector(VECSXP, 3));
    PROTECT(attrib = allocVector(STRSXP, 3));

    SET_STRING_ELT(attrib, 0, mkChar("R"));
    SET_STRING_ELT(attrib, 1, mkChar("U"));
    SET_STRING_ELT(attrib, 2, mkChar("A"));
    SET_VECTOR_ELT(count, 0, countR);
    SET_VECTOR_ELT(count, 1, countU);
    SET_VECTOR_ELT(count, 2, countA);
    setAttrib(count, R_NamesSymbol, attrib);
 
    // clean up
    samclose(fin);
    bam_index_destroy(idx);
    
    UNPROTECT(5);
    
    return count;
}
