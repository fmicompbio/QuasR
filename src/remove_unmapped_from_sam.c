#include "remove_unmapped_from_sam.h"

/*! @function
  @abstract  Walk through samfile and write all mapped alignments and the corresponding mate alignments to the output file. This function can also be used for a single-end alignments.
  @param  inSam  character(1) Input name of the samfile
  @param  outBam  character(1) Output name of the bamfile
  @param  paired  logical(1) TRUE if paired-end and FALSE if single-end
  @return    character vector of the output filename
 */
SEXP remove_unmapped_from_sam_and_convert_to_bam(SEXP inSam, SEXP outBam)
{
    // check parameters
    if(!Rf_isString(inSam) || 1 != Rf_length(inSam)){
        Rf_error("'inSam' must be character(1)");
    }
    
    if(!Rf_isString(outBam) || 1 < Rf_length(outBam)){
        Rf_error("'outBam' must be character(1).");
    }

    samfile_t *fin = samopen(Rf_translateChar(STRING_ELT(inSam, 0)), "r", NULL);
    if (fin == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", Rf_translateChar(STRING_ELT(inSam, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("SAM/BAM header missing or empty file: '%s'", Rf_translateChar(STRING_ELT(inSam, 0)));
    }

    samfile_t *fout = samopen(Rf_translateChar(STRING_ELT(outBam, 0)), "wb", fin->header);
    if(fout == 0)
        Rf_error("Error in opening the output file %s", Rf_translateChar(STRING_ELT(outBam, 0)));

    int r, count = 0;
    bam1_t *aln = bam_init1();

    while (0 <= (r = samread(fin, aln))){
        //check if mapped or corresponding mate is mapped
        if((aln->core.flag & BAM_FUNMAP) == 0 || ((aln->core.flag & BAM_FPAIRED) != 0  && (aln->core.flag & BAM_FMUNMAP) == 0))
            samwrite(fout, aln);
        count++;
    }

    // clean up
    bam_destroy1(aln);
    samclose(fin);
    samclose(fout);

    return outBam;
}
