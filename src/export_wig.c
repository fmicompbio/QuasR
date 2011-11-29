/*!
  @header

  TODO.
 
  @author:    Michael Stadler
  @date:      2011-11-24
  @copyright: Friedrich Miescher Institute for Biomedical Research, Switzerland
  @license: GPLv3
 */

#include "quasr.h"
#include <stdlib.h>
#include <math.h>


/*
  remark: could support compression by linking to the zlibbioc package
  see: http://www.bioconductor.org/packages/2.9/bioc/vignettes/zlibbioc/inst/doc/UsingZlibbioc.pdf
*/
SEXP bamfile_to_wig(SEXP bam_in, SEXP wig_out, SEXP width, SEXP shift, SEXP maxHits, SEXP normFactor, SEXP trackname, SEXP color, SEXP append) {
    if (!isString(bam_in)    || 1 != length(bam_in))
        error("'bam_in' must be character(1)");
    if (!isString(wig_out)   || 1 != length(wig_out))
        error("'wig_out' must be character(1)");
    if (!isInteger(width)    || 1 != length(width))
        error("'width' must be integer(1)");
    if (!isInteger(shift)    || 1 != length(shift))
        error("'shift' must be integer(1)");
    if (!isInteger(maxHits)  || 1 != length(maxHits))
        error("'maxHits' must be integer(1)");
    if (!isReal(normFactor)  || 1 != length(normFactor))
        error("'normFactor' must be numeric(1)");
    if (!isString(trackname) || 1 != length(trackname))
        error("'trackname' must be character(1)");
    if (!isString(color) || 1 != length(color))
        error("'color' must be character(1)");
    if (!isLogical(append)   || 1 != length(append))
        error("'append' must be logical(1)");
    
    int i=0, w=asInteger(width), s=asInteger(shift), m=asInteger(maxHits), beg=0, end=0x7fffffff;
    int32_t ih=0, currTid=0, currTlen=0, currTbin=0, hitPos=0;
    double *currCount=NULL, f=asReal(normFactor);
    int a=asLogical(append);
    const char *bi=translateChar(asChar(bam_in));
    const char *wo=translateChar(asChar(wig_out));
    const char *tn=translateChar(asChar(trackname));
    const char *co=translateChar(asChar(color));
    bam1_t *hit = bam_init1();

    if(w == NA_INTEGER) error("invalid '%s' value", "width");
    if(s == NA_INTEGER) error("invalid '%s' value", "shift");
    if(m == NA_INTEGER) error("invalid '%s' value", "maxHits");
    if(f == NA_REAL)    error("invalid '%s' value", "normFactor");
    if(a == NA_LOGICAL) error("'%s' must be TRUE or FALSE", "append");

    samfile_t *fin = _bam_tryopen(bi, "rb", NULL);
    if (fin->header == 0) {
        samclose(fin);
        error("invalid BAM file header");
    }

    bam_index_t *idx = bam_index_load(bi); // load BAM index
    if (idx == 0) // index is unavailable
        error("BAM index unavailable");

    FILE *fout;
    if (a)
	fout = fopen(wo, "a");
    else
	fout = fopen(wo, "w");
    if (fout == NULL)
	error("could not create outfile: %s\n", wo);

    // output track header
    fprintf(fout, "track type=wiggle_0 name='%s' description='%s' visibility=full color='%s' altColor='%s' priority=100 autoscale=off gridDefault=on maxHeightPixels=128:128:11 graphType=bar yLineMark=0.0 yLineOnOff=off windowingFunction=maximum smoothingWindow=off\n", tn, tn, co, co);

    // start a new target
    currTid = 0;                                 // is always first in a sorted bam file
    currTlen = fin->header->target_len[currTid]; // length of the current target
    currTbin = ceil((double)currTlen /w);        // number of bins on current target
    currCount = (double*)calloc(currTbin,sizeof(double)); // sum of alignment weights per window on current target
    fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", fin->header->target_name[currTid], w, w);

    while( samread(fin, hit) > 0) {
	if (hit->core.flag & BAM_FUNMAP) // skip unmapped reads
	    continue;

	// get alignment inverse weight (aux tag "IH")
	ih = get_inverse_weight(hit);
	if (ih > m || ih==0) // skip alignments of reads with more than maxHits (m) hits
	    continue;

	// get (shifted) position of hit (bam positions and 'hitPos' are zero-based)
	if (hit->core.flag & BAM_FREVERSE) // alignment on minus strand
	    hitPos = bam_calend(&(hit->core), bam1_cigar(hit)) - s;
	else                              // alignment on plus strand
	    hitPos = hit->core.pos + s;
	if (hitPos < 0)
	    hitPos = 0;
	if (hitPos >= currTlen)
	    hitPos = currTlen-1;

	if (currTid != hit->core.tid) {
	    // output former target
	    for(i=0; i<currTbin; i++)
		fprintf(fout, "%.1f\n", currCount[i]);
	    // start a new target
	    currTid = hit->core.tid;
	    currTlen = fin->header->target_len[currTid];
	    currTbin = ceil((double)currTlen /w);
	    fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", fin->header->target_name[currTid], w, w);
	    free(currCount);
	    currCount = (double*)calloc(currTbin,sizeof(double));
	}

	// add hit to current window
	currCount[ (int)floor((double)hitPos /w) ] += f/ih;
    }

    // output last target
    for(i=0; i<currTbin; i++)
	fprintf(fout, "%.1f\n", currCount[i]);

    free(currCount);
    fclose(fout);
    bam_index_destroy(idx);
    samclose(fin);
    bam_destroy1(hit);

    return wig_out;
}
