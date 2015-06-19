#include <stdlib.h>
#include <math.h>
#include "zlib.h"
#include "export_wig.h"


#define NO_ISIZE_FILTER -1   // disabled insert size-based alignment filtering


typedef struct {
    int bs;              // binsize
    int32_t cTid;        // current target id
    int32_t cTlen;       // current target length
    int32_t cTbin;       // number of bins on current target
    long unsigned int *count; // bin counts for current target
    int shift;           // shift for alignments
    int paired;          // paired experiment
    const char * strand; // alignment strand
    int log2p1;          // output log2(x+1)?
    uint8_t mapqMin;     // minimum mapping quality (MAPQ >= mapqMin)
    uint8_t mapqMax;     // maximum mapping quality (MAPQ <= mapqMax)
    int32_t absIsizeMin; // minimum absolute insert size (abs(ISIZE) >= absIsizeMin)
    int32_t absIsizeMax; // maximum absolute isnert size (abs(ISIZE) <= absIsizeMax)
} targetCoverage;


// call-back function for bam_fetch()
static int _addHitToCoverage(const bam1_t *hit, void *data){
    targetCoverage *tcov = (targetCoverage*)data;
    static int32_t hitPos;
    static int32_t hitBin;
    
    // skip alignment if region is not * and the strand of alignment or region is not the same
    if(strcmp(tcov->strand, "*")
       && (((hit->core.flag & BAM_FREVERSE) == 0) != (strcmp(tcov->strand, "+") == 0)))
        return 0;
    
    // skip alignment if mapping quality not in specified range
    if(hit->core.qual < tcov->mapqMin || hit->core.qual > tcov->mapqMax)
        return 0;
    
    // skip alignment if insert size not in specified range
    if((tcov->absIsizeMin != NO_ISIZE_FILTER && abs(hit->core.isize) < tcov->absIsizeMin) ||
       (tcov->absIsizeMax != NO_ISIZE_FILTER && abs(hit->core.isize) > tcov->absIsizeMax))
        return 0;
    
    if(tcov->paired) {
        if ((hit->core.flag & BAM_FPROPER_PAIR) && // skip reads that are not aligned as a proper pair
            !(hit->core.flag & BAM_FREAD2)) {      // skip read2 of proper pairs
            // get (shifted) position of hit (bam positions and 'hitPos' are zero-based)
	        if (hit->core.flag & BAM_FREVERSE) // alignment on minus strand
	            hitPos = (int)((double)bam_calend(&(hit->core), bam1_cigar(hit)) -1 + ((double)hit->core.isize - hit->core.isize/abs(hit->core.isize)) /2);
            else                               // alignment on plus strand
                hitPos = (int)((double)hit->core.pos + ((double)hit->core.isize - hit->core.isize/abs(hit->core.isize)) /2);
        } else {
	    return 0;
	}
    } else {
        if (!(hit->core.flag & BAM_FUNMAP)) { // skip unmapped reads
            // get (shifted) position of hit (bam positions and 'hitPos' are zero-based)
	        if (hit->core.flag & BAM_FREVERSE) // alignment on minus strand
		        hitPos = (int32_t)(bam_calend(&(hit->core), bam1_cigar(hit))) - 1 - tcov->shift;
    	    else                               // alignment on plus strand
    	        hitPos = hit->core.pos + tcov->shift;
        } else {
	    return 0;
	}
    }

    hitBin = (int32_t)floor((double)hitPos /tcov->bs);
    if (hitBin < 0 || hitBin >= tcov->cTbin) // skip out of chromosome positions (should be only last partial bin)
	return 0;

    // add hit to current window
    tcov->count[ hitBin ]++;

    return 0;
}


void output_current_target(targetCoverage *tcov, int compr, double fact, gzFile gzfh, FILE *fh) {
    int j;
    if(compr) {
	if(tcov->log2p1) {
	    for(j=0; j < tcov->cTbin; j++) {
		gzprintf(gzfh, "%.2f\n", log2((double)(tcov->count[j]) * fact + 1));
	    }
	} else {
	    for(j=0; j < tcov->cTbin; j++) {
		gzprintf(gzfh, "%.2f\n", (double)(tcov->count[j]) * fact);
	    }
	}
    } else {
	if(tcov->log2p1) {
	    for(j=0; j < tcov->cTbin; j++) {
		fprintf(fh, "%.2f\n", log2((double)(tcov->count[j]) * fact + 1));
	    }
	} else {
	    for(j=0; j < tcov->cTbin; j++) {
		fprintf(fh, "%.2f\n", (double)(tcov->count[j]) * fact);
	    }
	}
    }
}


void start_new_target(targetCoverage *tcov, bam_header_t *bh, int compr, gzFile gzfh, FILE *fh) {
    // assumes that tcov->cTid has been set to the new value
    tcov->cTlen = (int32_t)(bh->target_len[tcov->cTid]);
    tcov->cTbin = (int32_t)(floor((double)(tcov->cTlen) /tcov->bs)); // drop last partial bin on chromosome for compatibility with bigWig format
    if(compr)
	gzprintf(gzfh, "fixedStep chrom=%s start=1 step=%d span=%d\n", bh->target_name[tcov->cTid], tcov->bs, tcov->bs);
    else
	fprintf(fh, "fixedStep chrom=%s start=1 step=%d span=%d\n", bh->target_name[tcov->cTid], tcov->bs, tcov->bs);
    if(tcov->count != NULL)
        R_Free(tcov->count);
    tcov->count = (unsigned long int*) R_Calloc(tcov->cTbin, unsigned long int);
}


/* from one or several input bam files, produce a single, one-track wig file */
SEXP bamfile_to_wig(SEXP _bam_in, SEXP _wig_out, SEXP _paired, SEXP _binsize, SEXP _shift,
                    SEXP _strand, SEXP _norm_factor, SEXP _tracknames, SEXP _log2p1,
                    SEXP _colors, SEXP _compress, SEXP mapqMin, SEXP mapqMax, SEXP absIsizeMin, SEXP absIsizeMax) {
    // validate parameters
    if (!Rf_isString(_bam_in))
        Rf_error("'_bam_in' must be a character vector");
    if (!Rf_isString(_wig_out) || 1 != Rf_length(_wig_out))
        Rf_error("'_wig_out' must be a character(1)");
    if (!Rf_isLogical(_paired) || 1 != Rf_length(_paired))
        Rf_error("'_paired' must be logical(1)");
    if (!Rf_isInteger(_binsize) || 1 != Rf_length(_binsize))
        Rf_error("'_binsize' must be an integer(1)");
    if (!Rf_isInteger(_shift) || 1 != Rf_length(_shift))
        Rf_error("'_shift' must be an integer(1)");
    if (!Rf_isString(_strand) || 1 != Rf_length(_strand))
        Rf_error("'_strand' must be a character(1)");
    if (!Rf_isNumeric(_norm_factor) || 1 != Rf_length(_norm_factor))
        Rf_error("'_norm_factor' must be a numerical(1)");
    if (!Rf_isString(_tracknames) || 1 != Rf_length(_tracknames))
        Rf_error("'_tracknames' must be a character(1)");
    if (!Rf_isLogical(_log2p1) || 1 != Rf_length(_log2p1))
        Rf_error("'_log2p1' must be logical(1)");
    if (!Rf_isString(_colors) || 1 != Rf_length(_colors))
        Rf_error("'_colors' must be a character(1)");
    if (!Rf_isLogical(_compress) || 1 != Rf_length(_compress))
	Rf_error("'_compress' must be a logical(1)");
    if (!Rf_isInteger(mapqMin) || Rf_length(mapqMin) !=1 || INTEGER(mapqMin)[0] < 0 || INTEGER(mapqMin)[0] > 255)
        Rf_error("'mapqMin' must be of type integer(1) and have a value between 0 and 255");
    if (!Rf_isInteger(mapqMax) || Rf_length(mapqMax) !=1 || INTEGER(mapqMax)[0] < 0 || INTEGER(mapqMax)[0] > 255)
        Rf_error("'mapqMax' must be of type integer(1) and have a value between 0 and 255");
    if(INTEGER(mapqMin)[0] > INTEGER(mapqMax)[0])
	Rf_error("'mapqMin' must not be greater than 'mapqMax'");
    if(!Rf_isInteger(absIsizeMin) || Rf_length(absIsizeMin) !=1 || (INTEGER(absIsizeMin)[0] < 0 && INTEGER(absIsizeMin)[0] != NO_ISIZE_FILTER))
        Rf_error("'absIsizeMin' must be of type integer(1) and have a value greater than zero");
    if(!Rf_isInteger(absIsizeMax) || Rf_length(absIsizeMax) !=1 || (INTEGER(absIsizeMax)[0] < 0 && INTEGER(absIsizeMax)[0] != NO_ISIZE_FILTER))
        Rf_error("'absIsizeMax' must be of type integer(1) and have a value greater than zero");
    if(INTEGER(absIsizeMin)[0] != NO_ISIZE_FILTER && INTEGER(absIsizeMax)[0] != NO_ISIZE_FILTER && INTEGER(absIsizeMin)[0] > INTEGER(absIsizeMax)[0])
	Rf_error("'absIsizeMin' must not be greater than 'absIsizeMax'");
   

    // declare internal variables
    int t, i, n = Rf_length(_bam_in);
    double norm_factor = REAL(_norm_factor)[0];
    int compress = Rf_asLogical(_compress);
    const char **bam_in = (const char**) R_Calloc(n, char*);
    const char *wig_out = Rf_translateChar(STRING_ELT(_wig_out, 0));
    const char *tracknames = Rf_translateChar(STRING_ELT(_tracknames, 0));
    const char  *colors = Rf_translateChar(STRING_ELT(_colors, 0));
    for(i=0; i<n; i++)
        bam_in[i] = Rf_translateChar(STRING_ELT(_bam_in, i));

    targetCoverage tcov;
    tcov.bs = INTEGER(_binsize)[0];
    tcov.shift = INTEGER(_shift)[0];
    tcov.paired = Rf_asLogical(_paired);
    tcov.count = NULL;
    tcov.strand = Rf_translateChar(STRING_ELT(_strand, 0));
    tcov.log2p1 = Rf_asLogical(_log2p1);
    tcov.mapqMin = (uint8_t)(INTEGER(mapqMin)[0]);
    tcov.mapqMax = (uint8_t)(INTEGER(mapqMax)[0]);
    tcov.absIsizeMin = (uint32_t)(INTEGER(absIsizeMin)[0]);
    tcov.absIsizeMax = (uint32_t)(INTEGER(absIsizeMax)[0]);
 

    // open bam input files
    samfile_t **fin = (samfile_t**) R_Calloc(n, samfile_t*);
    bam_index_t **idx = (bam_index_t**) R_Calloc(n, bam_index_t*);
    bam1_t *hit = bam_init1();
    for(i=0; i<n; i++) {
        fin[i] = _bam_tryopen(bam_in[i], "rb", NULL);
    	if (fin[i]->header == 0) {
	    samclose(fin[i]);
	    Rf_error("invalid bam file header for %s", bam_in[i]);
	}
        idx[i] = bam_index_load(bam_in[i]);
	if (idx == 0) // index is unavailable
	    Rf_error("bam index unavailable for %s", bam_in[i]);
    }


    // open output file
    gzFile gzfout = NULL;
    FILE *fout = NULL;
    if(compress) {
	gzfout = gzopen(wig_out, "wb9");
	if(gzfout == NULL)
	    Rf_error("could not create compressed output file: %s",wig_out);
    } else {
	fout = fopen(wig_out, "w");
	if (fout == NULL)
	    Rf_error("could not create output file: %s", wig_out);
    }
    

    // output track header
    if(compress)
	gzprintf(gzfout, "track type=wiggle_0 name='%s' description='%s%s' visibility=full color=%s altColor=%s priority=100 autoscale=off gridDefault=on maxHeightPixels=64:64:11 graphType=bar yLineMark=0.0 yLineOnOff=off windowingFunction=maximum smoothingWindow=off\n", tracknames, tracknames, norm_factor==1.0 ? "" : " (scaled)", colors, colors);
    else
	fprintf(fout, "track type=wiggle_0 name='%s' description='%s%s' visibility=full color=%s altColor=%s priority=100 autoscale=off gridDefault=on maxHeightPixels=64:64:11 graphType=bar yLineMark=0.0 yLineOnOff=off windowingFunction=maximum smoothingWindow=off\n", tracknames, tracknames, norm_factor==1.0 ? "" : " (scaled)", colors, colors);


    // select bam_fetch callback function
    //bam_fetch_f fetch_func = _addHitToCoverage;
 
    // loop over targets
    for(t=0; t<fin[0]->header->n_targets; t++) {
        // start new target
        tcov.cTid = t; // is always first in a sorted bam file
        start_new_target(&tcov, fin[0]->header, compress, gzfout, fout);
        
        // loop over input bam files for current target --> sum coverage
        for(i=0; i<n; i++)
            bam_fetch(fin[i]->x.bam, idx[i], t, 0, tcov.cTlen, &tcov, (bam_fetch_f)_addHitToCoverage);
	//            bam_fetch(fin[i]->x.bam, idx[i], t, 0, tcov.cTlen, &tcov, fetch_func);

        // output current target
        output_current_target(&tcov, compress, norm_factor, gzfout, fout);
    }

    // close output file
    if(compress)
	gzclose(gzfout);
    else
	fclose(fout);


    // clean up
    bam_destroy1(hit);
    for(i=0; i<n; i++) {
        bam_index_destroy(idx[i]);
    	samclose(fin[i]);
    }
    R_Free(tcov.count);
    R_Free(bam_in);
    R_Free(fin);
    R_Free(idx);

    return R_NilValue;
}

