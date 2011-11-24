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


/*
  remark: could support compression by linking to the zlibbioc package
  see: http://www.bioconductor.org/packages/2.9/bioc/vignettes/zlibbioc/inst/doc/UsingZlibbioc.pdf
*/
SEXP bamfile_to_wig(SEXP bam_in, SEXP wig_out, SEXP width) {
    if (!IS_CHARACTER(bam_in) || 1 != LENGTH(bam_in))
        Rf_error("'bam_in' must be character(1)");
    if (!IS_CHARACTER(wig_out) || 1 != LENGTH(wig_out))
        Rf_error("'wig_out' must be character(1)");
    if (!IS_INTEGER(width) || 1 != LENGTH(width))
        Rf_error("'width' must be integer(1)");
    
    int w=asInteger(width), beg=0, end=0x7fffffff;
    int32_t ih=0, currTid=0, currPos=0;
    double currCount=0.0;  
    bam1_t *hit = bam_init1();

    samfile_t *fin = _bam_tryopen(translateChar(asChar(bam_in)), "rb", NULL);
    if (fin->header == 0) {
        samclose(fin);
        error("invalid header");
    }

    bam_index_t *idx = bam_index_load(translateChar(asChar(bam_in))); // load BAM index
    if (idx == 0) // index is unavailable
        error("BAM index unavailable");

    FILE *fout = fopen(translateChar(asChar(wig_out), "w");
    if (fout == NULL)
	error("could not create outfile: %s\n", wig_out);

    // start a new target
    currTid = 0;     // is always first in a sorted bam file
    currPos = w-1;   // end of current window (zero-based)
    currCount = 0.0; // sum of alignment weight in current window
    fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", fin->header->target_name[currTid], w, w);

    while( samread(fin, hit) > 0) {
	if (hit->core.flag & BAM_FUNMAP) // skip unmapped reads
	    continue;

	// get alignment inverse weight (aux tag "IH")
	ih = _get_inverse_weight(hit);
	/*
	  int32_t hit->core.pos : 0-based leftmost coordinate
	  int32_t hit->core.tid : chromosome ID, defined by bam_header_t
	  uint32_t bam_calend() : calculate the rightmost coordinate of an alignment on the reference genome
	*/

	if (currTid != hit->core.tid) {
	    // output last window on former target
	    fprintf(fout, "%.1f\n", currCount);
	    // start a new target
	    currTid = hit->core.tid;
	    currPos = w-1;
	    fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", fin->header->target_name[currTid], w, w);
	    // output empty windows
	    while ((currPos += w) < hit->core.pos)
		fprintf(fout, "0\n");
	    // start new window
	    currCount = 0.0;

	} else	if (currPos < hit->core.pos) {
	    // output last window
	    fprintf(fout, "%.1f\n", currCount);
	    // output empty windows
	    while ((currPos += w) < hit->core.pos)
		fprintf(fout, "0\n");
	    // start new window
	    currCount = 0.0;
	}

	// add hit to current window
	currCount += 1.0/ih;
    }
    // output very last window
    fprintf(fout, "%.1f\n", currCount);

    fclose(fout);
    bam_index_destroy(idx);
    samclose(fin);
    bam_destroy1(hit);

    return wig_out;
}
