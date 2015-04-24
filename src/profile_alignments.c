/*!
  @header
  
  Counts the alignments in regions, which fulfill specific criteria.
  
  @author:    Michael Stadler
  @date:      2011-08-17
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "profile_alignments.h"
#include <stdlib.h>

#define SMART_SHIFT -1000000 // for half insert size shift towards the mate read
#define NO_ISIZE_FILTER -1   // disabled insert size-based alignment filtering

/*! @typedef
  @abstract Structure to provid the data to the bam_fetch() functions.
  @field sumU         int[] over position with counts for the fetch region (non-allelic or allelic Unknown)
  @field sumR         int[] over position with counts for the fetch region (allelic Reference)
  @field sumA         int[] over position with counts for the fetch region (allelic Alternative)
  @field offset       int offset of anchor position within sumU/sumR/sumA arrays
  @field len          int length of sumU/sumR/sumA arrays
  @field start        start position of the fetch region
  @field end          end position of the fetch region
  @field ref          reference postition of the fetch region
  @field selstrand    strand of the fetch region (used to select alignments; have to be on the same strand)
  @field regstrand    strand of the fetch region (controls calculation of relative position: from the left(+,*) or the right(-)
  @field shift        shift size for alignment
  @field readBitMask  select first/last read from a multi-read experiment
  @field selectReadPosition  weight of alignment on "s"tart or "e"nd
  @field allelic      allelic true(1) or false(0)
  @field includeSpliced  count spliced alignments true(1) or false(0)
  @field mapqMin      minimum mapping quality (MAPQ >= mapqMin)
  @field mapqMax      maximum mapping quality (MAPQ <= mapqMax)
  @field absIsizeMin  minimum absolute insert size (abs(ISIZE) >= absIsizeMin)
  @field absIsizeMax  maximum absolute isnert size (abs(ISIZE) <= absIsizeMax)
*/
typedef struct {
    int *sumU;
    int *sumR;
    int *sumA;
    int offset;
    int len;
    int start;
    int end;
    int ref;
    const char *selstrand;
    const char *regstrand;
    int shift;
    int readBitMask;
    char selectReadPosition;
    int allelic;
    int includeSpliced;
    uint8_t mapqMin;
    uint8_t mapqMax;
    int32_t absIsizeMin;
    int32_t absIsizeMax;
} regionProfile;


/*! @function
  @abstract  increase allele specific counts.
  @param     hit     the alignment
  @param     rinfo   regionProfile
  @param     relpos  relative position within sumU/sumR/sumA
  @return    0 if successful
 */
static int _sum_allelic(const bam1_t *hit, regionProfile *rinfo, int relpos){
    static uint8_t *xv_ptr = 0;

    // get XV tag
    xv_ptr = bam_aux_get(hit,"XV");
    if(xv_ptr == 0)
        Rf_error("XV tag missing but needed for allele-specific counting");

    // increase count
    switch(bam_aux2A(xv_ptr)){
    case 'U':
        rinfo->sumU[relpos] += 1;
        break;
    case 'R':
        rinfo->sumR[relpos] += 1;
        break;
    case 'A':
        rinfo->sumA[relpos] += 1;
        break;
    default:
        Rf_error("'%c' is not a valid XV tag value; should be one of 'U','R' or 'A'", bam_aux2A(xv_ptr));
    }

    return 0;
}


/*! @function
  @abstract  callback for bam_fetch(); sums alignments by position if correct strand and shifted biological start/end position overlap fetch region
  @param  hit   the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _addValidHitToSums(const bam1_t *hit, void *data){
    regionProfile *rinfo = (regionProfile*)data;
    static double shift = 0;
    static int pos = 0, relpos = 0;

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

    // skip alignment if region is not * and the strand of alignment or region is not the same
    if(strcmp(rinfo->selstrand, "*")
       && (((hit->core.flag & BAM_FREVERSE) == 0) != (strcmp(rinfo->selstrand, "+") == 0)))
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

    // calculate absolute position
    if(((hit->core.flag & BAM_FREVERSE) == 0) == (rinfo->selectReadPosition == 's')) // XOR
	// plus-strand and startwithin OR minus-strand and endwithin
        // --> position on the left side of the read 
	pos = (int)((double)hit->core.pos + shift); // 0-based inclusive start
    else
	// plus-strand and endwithin OR minus-strand and startwithin
        // --> position on the right side of the read 
        pos = (int)((double)bam_calend(&hit->core, bam1_cigar(hit)) - 1 + shift); // 0-based exclusive end --> -1

    // calculate relative position
    if(strcmp(rinfo->regstrand, "-"))
	// plus-strand or unstranded region --> measure from the left
	relpos = pos - rinfo->ref + rinfo->offset;
    else
	// minus-strand region --> measure from the right
	relpos = rinfo->ref - pos + rinfo->offset;

    // check if in region
    if(rinfo->start <= pos && pos < rinfo->end &&
       relpos >= 0 && relpos < rinfo->len){
        if(rinfo->allelic == 0)
            rinfo->sumU[relpos] += 1;
        else
            _sum_allelic(hit, rinfo, relpos);
    }

    return 0;
}


/*! @function
  @abstract  verify the parameters of the profile_alignments_non_allelic and profile_alignments_allelic function
  @return       0 if successful
 */
int _verify_profile_parameters(SEXP bamfile, SEXP profileids, SEXP tid,  SEXP start, SEXP end, SEXP refpos,
                               SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                               SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                               SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){
    // check bamfile parameter
    if(!Rf_isString(bamfile) || Rf_length(bamfile) != 1)
        Rf_error("'bamfile' must be of type character(1)");

    // check region parameters
    if(!Rf_isInteger(profileids))
        Rf_error("'profileids' must be of type integer");
    if(!Rf_isInteger(tid))
        Rf_error("'tid' must be of type integer");
    if(!Rf_isInteger(start))
        Rf_error("'start' must be of type integer");
    if(!Rf_isInteger(end))
        Rf_error("'end' must be of type integer");
    if(!Rf_isInteger(refpos))
        Rf_error("'refpos' must be of type integer");
    if(!Rf_isString(selstrand))
        Rf_error("'selstrand' must be of type character");
    if(!Rf_isString(regstrand))
        Rf_error("'regstrand' must be of type character");
    int num_regions = Rf_length(profileids);
    if(num_regions != Rf_length(tid) || num_regions != Rf_length(start) || num_regions != Rf_length(end) ||
       num_regions != Rf_length(refpos) || num_regions != Rf_length(selstrand) || num_regions != Rf_length(regstrand))
        Rf_error("'tid', 'start', 'end', 'refpos', 'selstrand' and 'regstrand' must have equal length");

    // check selectReadPosition
    if(!Rf_isString(selectReadPosition) || Rf_length(selectReadPosition) != 1)
        Rf_error("'selectReadPosition' must be of type character(1)");
    if(Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0] != 's' 
       && Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0] != 'e')
        Rf_error("The value of 'selectReadPosition' not supportet.");

    // check parameters readBitMask, shift, broaden, maxWidth, includeSpliced
    if(!Rf_isInteger(readBitMask) || Rf_length(readBitMask) != 1)
        Rf_error("'readBitMask' must be of type integer(1)");
    if(!Rf_isInteger(shift) && Rf_length(shift) != 1)
        Rf_error("'shift' must be of type integer(1)");
    if(!Rf_isInteger(broaden) && Rf_length(broaden) != 1)
        Rf_error("'broaden' must be of type integer(1)");
    if(INTEGER(broaden)[0] < 0)
        Rf_error("'broaden' must be a positive value.");
    if(!Rf_isInteger(maxUp) && Rf_length(maxUp) != 1)
        Rf_error("'maxUp' must be of type integer(1)");
    if(!Rf_isInteger(maxDown) && Rf_length(maxDown) != 1)
        Rf_error("'maxDown' must be of type integer(1)");
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
  @abstract  Counts alignments by position in regions (orientation in regstrand), selecting on selstrand and shift
  @param  bamfile   name of the bamfile
  @param  profileids       profile integer identifier (monotoniously increasing)
                                - target regions (see tid, start, end, refpos, selstrand and regstrand) with identical
                                  profile integer identifiers will be combine into the same profile
                                - the number of profiles returned corresponds to profileids[len-1]-profileids[0]+1
                                  (which is the number of unique integer identifiers in profileids)
  @param  tid                 target region identifier
  @param  start               target region start
  @param  end                 target region end
  @param  refpos              target region anchor position
  @param  selstrand           target region strand (used to select alignments; have to be on the same strand)
  @param  regstrand           target region strand (controls how relative position is calculated, from left(+,*) or right(-) side)
  @param  selectReadPosition  alignment ancored at start/end
  @param  readBitMask         select first/second/any read in a paired-end experiment
  @param  shift               shift size
  @param  broaden             extend query region for bam_fetch to catch alignments with overlaps due to shifting
  @param  maxUp               maximal upstream length of region
  @param  maxDown             maximal downstream length of region
  @param  includeSpliced      also count spliced alignments
  @param  mapqMin             minimal mapping quality to count alignment (MAPQ >= mapqMin)
  @param  mapqMax             maximum mapping quality to count alignment (MAPQ <= mapqMax)
  @param  absIsizeMin         minimum absolute insert size (abs(ISIZE) >= absIsizeMin)
  @param  absIsizeMax         maximum absolute isnert size (abs(ISIZE) <= absIsizeMax)
  @return          vector of length maxWidth with alignment counts per relative position in regions
 */
SEXP profile_alignments_non_allelic(SEXP bamfile, SEXP profileids, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                    SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                    SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                    SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){

    // check parameters
    _verify_profile_parameters(bamfile, profileids, tid, start, end, refpos, selstrand, regstrand,
                               selectReadPosition, readBitMask, shift, broaden, maxUp, maxDown, includeSpliced,
			       mapqMin, mapqMax, absIsizeMin, absIsizeMax);
    
    // open bam file
    samfile_t *fin = 0;
    fin = samopen(Rf_translateChar(STRING_ELT(bamfile, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open BAM file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("BAM header missing or empty in file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (idx == 0){
        samclose(fin);
        Rf_error("failed to open BAM index for file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }

    // initialise regionProfile
    int i, *profU_c;
    int *profId = INTEGER(profileids);
    for(i=Rf_length(tid)-1; i>=0; i--)
        profId[i] = profId[i] - profId[0];
    int profIdNum = profId[Rf_length(tid)-1] + 1;
    int mu = INTEGER(maxUp)[0], md = INTEGER(maxDown)[0];
    int mw = mu+md+1;
    SEXP profU;
    PROTECT(profU = allocMatrix(INTSXP, mw, profIdNum));
    profU_c = INTEGER(profU);
    for(i=0; i<mw*profIdNum; i++)
        profU_c[i] = 0;

    regionProfile rprof;
    rprof.offset = mu;
    rprof.len = mw;
    rprof.shift = INTEGER(shift)[0];
    rprof.readBitMask = INTEGER(readBitMask)[0];
    rprof.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rprof.allelic = 0;
    rprof.includeSpliced = (Rf_asLogical(includeSpliced) ? 1 : 0);
    rprof.mapqMin = (uint8_t)(INTEGER(mapqMin)[0]);
    rprof.mapqMax = (uint8_t)(INTEGER(mapqMax)[0]);
    rprof.absIsizeMin = (uint32_t)(INTEGER(absIsizeMin)[0]);
    rprof.absIsizeMax = (uint32_t)(INTEGER(absIsizeMax)[0]);
   
    // set shift for fetch to zero if smart shift
    int shift_f = abs(INTEGER(shift)[0]);
    if(INTEGER(shift)[0] == SMART_SHIFT)
        shift_f = 0;
    
    // select bam_fetch callback function
    bam_fetch_f fetch_func = _addValidHitToSums;
    
    // loop over query regions
    for(i=0; i < Rf_length(tid); i++){
        // setup region in rprof
        rprof.sumU = profU_c + mw*profId[i];
        rprof.start = INTEGER(start)[i];
        rprof.end = INTEGER(end)[i];
        rprof.ref = INTEGER(refpos)[i];
        rprof.selstrand = Rf_translateChar(STRING_ELT(selstrand, i));
        rprof.regstrand = Rf_translateChar(STRING_ELT(regstrand, i));

        // process alignments that overlap region
        bam_fetch(fin->x.bam, idx, INTEGER(tid)[i], 
                  INTEGER(start)[i]-shift_f-INTEGER(broaden)[0], // 0-based inclusive start
                  INTEGER(end)[i]+shift_f+INTEGER(broaden)[0],   // 0-based exclusive end
                  &rprof, fetch_func);
    }
    
    // clean up
    samclose(fin);
    bam_index_destroy(idx);
    
    UNPROTECT(1);
    
    return profU;
}

SEXP profile_alignments_allelic(SEXP bamfile, SEXP profileids, SEXP tid, SEXP start, SEXP end, SEXP refpos,
                                SEXP selstrand, SEXP regstrand, SEXP selectReadPosition, SEXP readBitMask,
                                SEXP shift, SEXP broaden, SEXP maxUp, SEXP maxDown, SEXP includeSpliced,
                                SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax){
    // check parameters
    _verify_profile_parameters(bamfile, profileids, tid, start, end, refpos, selstrand, regstrand,
                               selectReadPosition, readBitMask, shift, broaden, maxUp, maxDown, includeSpliced,
			       mapqMin, mapqMax, absIsizeMin, absIsizeMax);
    
    // open bam file
    samfile_t *fin = 0;
    fin = samopen(Rf_translateChar(STRING_ELT(bamfile, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open BAM file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("BAM header missing or empty in file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (idx == 0){
        samclose(fin);
        Rf_error("failed to open BAM index for file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }

    // initialise regionProfile
    int i, *profU_c, *profR_c, *profA_c;
    int *profId = INTEGER(profileids);
    for(i=Rf_length(tid)-1; i>=0; i--)
        profId[i] = profId[i] - profId[0];
    int profIdNum = profId[Rf_length(tid)-1] + 1;
    int mu = INTEGER(maxUp)[0], md = INTEGER(maxDown)[0];
    int mw = mu+md+1;
    SEXP profU, profR, profA, prof, attrib;
    PROTECT(profU = allocMatrix(INTSXP, mw, profIdNum));
    PROTECT(profR = allocMatrix(INTSXP, mw, profIdNum));
    PROTECT(profA = allocMatrix(INTSXP, mw, profIdNum));
    profU_c = INTEGER(profU);
    profR_c = INTEGER(profR);
    profA_c = INTEGER(profA);
    for(i=0; i<mw*profIdNum; i++)
        profU_c[i] = profR_c[i] = profA_c[i] = 0;

    regionProfile rprof;
    rprof.offset = mu;
    rprof.len = mw;
    rprof.shift = INTEGER(shift)[0];
    rprof.readBitMask = INTEGER(readBitMask)[0];
    rprof.selectReadPosition = Rf_translateChar(STRING_ELT(selectReadPosition, 0))[0];
    rprof.allelic = 1;
    rprof.includeSpliced = (Rf_asLogical(includeSpliced) ? 1 : 0);
    rprof.mapqMin = (uint8_t)(INTEGER(mapqMin)[0]);
    rprof.mapqMax = (uint8_t)(INTEGER(mapqMax)[0]);
    rprof.absIsizeMin = (uint32_t)(INTEGER(absIsizeMin)[0]);
    rprof.absIsizeMax = (uint32_t)(INTEGER(absIsizeMax)[0]);

    // set shift for fetch to zero if smart shift
    int shift_f = abs(INTEGER(shift)[0]);
    if(INTEGER(shift)[0] == SMART_SHIFT)
        shift_f = 0;
    
    // select bam_fetch callback function
    bam_fetch_f fetch_func = _addValidHitToSums;
    
    // loop over query regions
    for(i=0; i < Rf_length(tid); i++){
        // setup region in rprof
        rprof.sumU = profU_c + mw*profId[i];
        rprof.sumR = profR_c + mw*profId[i];
        rprof.sumA = profA_c + mw*profId[i];
        rprof.start = INTEGER(start)[i];
        rprof.end = INTEGER(end)[i];
        rprof.ref = INTEGER(refpos)[i];
        rprof.selstrand = Rf_translateChar(STRING_ELT(selstrand, i));
        rprof.regstrand = Rf_translateChar(STRING_ELT(regstrand, i));

        // process alignments that overlap region
        bam_fetch(fin->x.bam, idx, INTEGER(tid)[i], 
                  INTEGER(start)[i]-shift_f-INTEGER(broaden)[0], // 0-based inclusive start
                  INTEGER(end)[i]+shift_f+INTEGER(broaden)[0],   // 0-based exclusive end
                  &rprof, fetch_func);
    }

    // pack results into list
    PROTECT(prof = allocVector(VECSXP, 3));
    PROTECT(attrib = allocVector(STRSXP, 3));

    SET_STRING_ELT(attrib, 0, mkChar("R"));
    SET_STRING_ELT(attrib, 1, mkChar("U"));
    SET_STRING_ELT(attrib, 2, mkChar("A"));
    SET_VECTOR_ELT(prof, 0, profR);
    SET_VECTOR_ELT(prof, 1, profU);
    SET_VECTOR_ELT(prof, 2, profA);
    setAttrib(prof, R_NamesSymbol, attrib);
 
    // clean up
    samclose(fin);
    bam_index_destroy(idx);
    
    UNPROTECT(5);
    
    return prof;
}
