/*!
  @header
  
  Counts the alignments in regions, which fulfill specific criteria.
  
  @author:    Anita Lerch, Michael Stadler
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "count_alignments.h"
#include <stdlib.h>

#define SMART_SHIFT -1000000 // for half insert size shift towards the mate read
#define NO_ISIZE_FILTER -1   // disabled insert size-based alignment filtering

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
  @field skipSecondary  skip secondary alignments true(1) or false(0)
  @field selectReadPosition  weight of alignment on "s"tart or "e"nd
  @field allelic      allelic true(1) or false(0)
  @field includeSpliced  count spliced alignments true(1) or false(0)
  @field mapqMin      minimum mapping quality (MAPQ >= mapqMin)
  @field mapqMax      maximum mapping quality (MAPQ <= mapqMax)
  @field absIsizeMin  minimum absolute insert size (abs(ISIZE) >= absIsizeMin)
  @field absIsizeMax  maximum absolute isnert size (abs(ISIZE) <= absIsizeMax)
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
    int skipSecondary;
    char selectReadPosition;
    int allelic;
    int includeSpliced;
    uint8_t mapqMin;
    uint8_t mapqMax;
    int32_t absIsizeMin;
    int32_t absIsizeMax;
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

    // skip alignment if rinfo->includeSpliced == false and alignmend is spliced
    if(rinfo->includeSpliced == 0 && _isSpliced(hit) == 1)
        return 0;
    
    // skip alignment if mapping quality not in specified range
    if(hit->core.qual < rinfo->mapqMin || hit->core.qual > rinfo->mapqMax)
        return 0;
    
    // skip alignment if insert size not in specified range
    if((rinfo->absIsizeMin != NO_ISIZE_FILTER && abs(hit->core.isize) < rinfo->absIsizeMin) ||
       (rinfo->absIsizeMax != NO_ISIZE_FILTER && abs(hit->core.isize) > rinfo->absIsizeMax))
        return 0;
    
    // skip alignment if read1 or read2 flag is set (=paired-end) and if wrong readBitMask
    if((hit->core.flag & (BAM_FREAD1 + BAM_FREAD2)) && (hit->core.flag & rinfo->readBitMask) == 0)
        return 0;

    // skip alignment if secondary and skipSecondary==true
    if((hit->core.flag & BAM_FSECONDARY) && rinfo->skipSecondary)
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
                      SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                      SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){
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

    // check parameter readBitMask, shift, broaden, includeSpliced
    if(!Rf_isInteger(readBitMask) || Rf_length(readBitMask) != 1)
        Rf_error("'readBitMask' must be of type integer(1)");
    if(!Rf_isInteger(shift) && Rf_length(shift) != 1)
        Rf_error("'shift' must be of type integer(1)");
    if(!Rf_isInteger(broaden) && Rf_length(broaden) != 1)
        Rf_error("'broaden' must be of type integer(1)");
    if(INTEGER(broaden)[0] < 0)
        Rf_error("'broaden' must be a positive value.");
    if(!Rf_isLogical(includeSpliced) || 1 != Rf_length(includeSpliced))
        Rf_error("'includeSpliced' must be of type logical(1)");
        
    // check MAPQ parameters
    if(!Rf_isInteger(mapqMin) || Rf_length(mapqMin) !=1 || INTEGER(mapqMin)[0] < 0 || INTEGER(mapqMin)[0] > 255)
        Rf_error("'mapqMin' must be of type integer(1) and have a value between 0 and 255");
    if(!Rf_isInteger(mapqMax) || Rf_length(mapqMax) !=1 || INTEGER(mapqMax)[0] < 0 || INTEGER(mapqMax)[0] > 255)
        Rf_error("'mapqMax' must be of type integer(1) and have a value between 0 and 255");
    if(INTEGER(mapqMin)[0] > INTEGER(mapqMax)[0])
	Rf_error("'mapqMin' must not be greater than 'mapqMax'");
        
    // check TLEN parameters
    if(!Rf_isInteger(absIsizeMin) || Rf_length(absIsizeMin) !=1 || (INTEGER(absIsizeMin)[0] < 0 && INTEGER(absIsizeMin)[0] != NO_ISIZE_FILTER))
        Rf_error("'absIsizeMin' must be of type integer(1) and have a value greater than zero");
    if(!Rf_isInteger(absIsizeMax) || Rf_length(absIsizeMax) !=1 || (INTEGER(absIsizeMax)[0] < 0 && INTEGER(absIsizeMax)[0] != NO_ISIZE_FILTER))
        Rf_error("'absIsizeMax' must be of type integer(1) and have a value greater than zero");
    if(INTEGER(absIsizeMin)[0] != NO_ISIZE_FILTER && INTEGER(absIsizeMax)[0] != NO_ISIZE_FILTER && INTEGER(absIsizeMin)[0] > INTEGER(absIsizeMax)[0])
	Rf_error("'absIsizeMin' must not be greater than 'absIsizeMax'");

    return 0;
}

/*! @function
  @abstract  Counts the alignments in regions, which fit the strand, overlap type and shift criteria.
  @param  bamfile             Name of the bamfile
  @param  tid                 target region identifier
  @param  start               target region start
  @param  end                 target region end
  @param  strand              target region strand
  @param  selectReadPosition  alignment ancored at start/end
  @param  readBitMask         select first/second/any read in a paired-end experiment; select secondary alignments
  @param  shift               shift size
  @param  broaden             extend query region for bam_fetch to catch alignments with overlaps due to shifting
  @param  includeSpliced      also count spliced alignments
  @param  mapqMin             minimal mapping quality to count alignment (MAPQ >= mapqMin)
  @param  mapqMax             maximum mapping quality to count alignment (MAPQ <= mapqMax)
  @param  absIsizeMin         minimum absolute insert size (abs(ISIZE) >= absIsizeMin)
  @param  absIsizeMax         maximum absolute isnert size (abs(ISIZE) <= absIsizeMax)
  @return               Vector of the alignment counts
 */
SEXP count_alignments_non_allelic(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP strand,
                                  SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                                  SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){

    // check parameters
    _verify_parameters(bamfile, tid, start, end, strand, selectReadPosition, readBitMask, shift, broaden, includeSpliced,
		       mapqMin, mapqMax, absIsizeMin, absIsizeMax);
    
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
    rinfo.readBitMask = (INTEGER(readBitMask)[0] & (BAM_FREAD1 + BAM_FREAD2));
    rinfo.skipSecondary = ((INTEGER(readBitMask)[0] & BAM_FSECONDARY) ? 0 : 1);
    rinfo.shift = INTEGER(shift)[0];
    rinfo.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rinfo.allelic = 0;
    rinfo.includeSpliced = (Rf_asLogical(includeSpliced) ? 1 : 0);
    rinfo.mapqMin = (uint8_t)(INTEGER(mapqMin)[0]);
    rinfo.mapqMax = (uint8_t)(INTEGER(mapqMax)[0]);
    rinfo.absIsizeMin = (uint32_t)(INTEGER(absIsizeMin)[0]);
    rinfo.absIsizeMax = (uint32_t)(INTEGER(absIsizeMax)[0]);
    
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
                                  SEXP selectReadPosition, SEXP readBitMask, SEXP shift, SEXP broaden, SEXP includeSpliced,
                                  SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){

    // check parameters
    _verify_parameters(bamfile, tid, start, end, strand, selectReadPosition, readBitMask, shift, broaden, includeSpliced,
		       mapqMin, mapqMax, absIsizeMin, absIsizeMax);
    
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
    rinfo.readBitMask = (INTEGER(readBitMask)[0] & (BAM_FREAD1 + BAM_FREAD2));
    rinfo.skipSecondary = ((INTEGER(readBitMask)[0] & BAM_FSECONDARY) ? 0 : 1);
    rinfo.shift = INTEGER(shift)[0];
    rinfo.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rinfo.allelic = 1;
    rinfo.includeSpliced = (Rf_asLogical(includeSpliced) ? 1 : 0);
    rinfo.mapqMin = (uint8_t)(INTEGER(mapqMin)[0]);
    rinfo.mapqMax = (uint8_t)(INTEGER(mapqMax)[0]);
    rinfo.absIsizeMin = (uint32_t)(INTEGER(absIsizeMin)[0]);
    rinfo.absIsizeMax = (uint32_t)(INTEGER(absIsizeMax)[0]);

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
