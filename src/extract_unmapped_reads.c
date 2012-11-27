#include "extract_unmapped_reads.h"
#include "utilities.h"

// get sequence string from an alignment
// function copied from io_sam.c of Rsamtools
char * _bamseq(const bam1_t * bam, int reverse_comp)
{
    static const char key[] = {
        '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    };

    const uint32_t len = bam->core.l_qseq;
    const unsigned char *seq = bam1_seq(bam);
    char *s = R_Calloc(len + 1, char);
    for (uint32_t i = 0; i < len; ++i)
        s[i] = key[bam1_seqi(seq, i)];
    if (reverse_comp){
        _complement(s, len);
        _reverse(s, len);
    }
    s[len] = '\0';
    return s;
}

// get quality string from an alignment
// function copied from io_sam.c of Rsamtools
char * _bamqual(const bam1_t * bam, int reverse)
{
    const uint32_t len = bam->core.l_qseq;
    const unsigned char *bamq = bam1_qual(bam);
    char *s = R_Calloc(len + 1, char);
    for (uint32_t i = 0; i < len; ++i)
        s[i] = bamq[i] + 33;
    if (reverse)
        _reverse(s, len);
    s[len] = '\0';
    return s;
}

// write id and sequence of an alignment into a file as fastq-format
void _write_fastq_seq(FILE * fp, const bam1_t * bam, int reverse_comp){
    char * seq = _bamseq(bam, reverse_comp);
    char * qual = _bamqual(bam, reverse_comp);
    fprintf(fp, "@%s\n%s\n+\n%s\n", bam1_qname(bam), seq, qual);
    R_Free(qual);
    R_Free(seq);
}

// write id and sequence of an alignment into a file as fasta-format
void _write_fasta_seq(FILE * fp, const bam1_t * bam, int reverse_comp){
    char * seq = _bamseq(bam, reverse_comp);
    fprintf(fp, ">%s\n%s\n", bam1_qname(bam), seq);
    R_Free(seq);
}

// get unmapped single-end reads and write them into a file
int _extract_unmapped_single_reads(samfile_t *fin, FILE *fout, int fastq){

    int r, count = 0;
    bam1_t *aln = bam_init1();
    while (0 <= (r = samread(fin, aln))) {
        // check if unmapped 
        if((aln->core.flag & BAM_FUNMAP) != 0){
            // check if fastq format
            if(fastq){
                _write_fastq_seq(fout, aln, 0);
            }else{
                _write_fasta_seq(fout, aln, 0);
            }
        }
        count++;
    }
    bam_destroy1(aln);
    return r >= -1 ? count : -1 * count;
}

// get unmapped paired-end reads and write them into a file
int _extract_unmapped_paired_reads(samfile_t *fin, FILE *fout1, FILE *fout2, int fastq, int rc_read2)
{
    int r, count = 0;
    bam1_t *aln1 = bam_init1();
    bam1_t *aln2 = bam_init1();
    while (0 <= (r = samread(fin, aln1))) {
        //check if unmapped and mate is unmapped too
        if((aln1->core.flag & BAM_FUNMAP) != 0 && (aln1->core.flag & BAM_FMUNMAP) != 0){
            // TODO possibly check flag BAM_FPAIRED too
            // read mate, mate must follow as next alignment
            if(0 <= (r = samread(fin, aln2))){
                // check if mate sequence
                if(((aln2->core.flag & BAM_FUNMAP) != 0) && ((aln2->core.flag & BAM_FMUNMAP) != 0)
                    && ((aln1->core.flag & BAM_FREAD1) != 0) && ((aln2->core.flag & BAM_FREAD2) != 0)){
                    // check if fastq format and write unmapped pairs
                    if(fastq){
                        _write_fastq_seq(fout1, aln1, 0);
                        _write_fastq_seq(fout2, aln2, rc_read2);
                    }else{
                        _write_fasta_seq(fout1, aln1, 0);
                        _write_fasta_seq(fout2, aln2, rc_read2);
                    }
                }else{
                    Rf_error("The order of unmapped paired-end sequences in bamfile is inconsistent at %i-th alignment.\n", count);
                }
            }else{ // read of mate not succesfull, stop the loop
                if(r >= -1)
                    Rf_error("The order of unmapped paired-end sequences in bamfile is inconsistent at EOF.\n");
                else
                    Rf_error("Reading failed after %i-th alignment.\n", count);
                break;
            }
            count++;
        }
        count++;
    }
    bam_destroy1(aln1);
    bam_destroy1(aln2);
    return r >= -1 ? count : -1 * count;
}

/*! @function
  @abstract Extract unmapped reads from a input bamfile and write the sequence to a fasta or fastq file.
  @param  inBam  character(1) Input name of the bamfile
  @param  outBam  character(1) Output name of the bamfile
  @param  fastq  logical(1) TRUE if output sequence should be in the fastq format or FALSE if in fasta format
  @param  rcRead2  logical(1) TRUE if second read should be reverse complemented
  @return    outBam  character(1) Output name of the bamfile
 */
SEXP extract_unmapped_reads(SEXP inBam, SEXP outFile, SEXP fastq, SEXP rcRead2)
{
    // check parameters
    if(!Rf_isString(inBam) || 1 != Rf_length(inBam)){
        Rf_error("'inBam' must be character(1)");
    }
    
    if(!Rf_isString(outFile) || 2 < Rf_length(outFile)){
        Rf_error("'outFile' must be character(1) if single-end and character(2) if paired-end.");
    }
    
    if(!Rf_isLogical(fastq) || 1 != Rf_length(fastq)){
        Rf_error("'fastq' must be logical(1)");
    }

    if(!Rf_isLogical(rcRead2) || 1 != Rf_length(rcRead2)){
        Rf_error("'rcRead2' must be logical(1)");
    }

    samfile_t *fin = samopen(Rf_translateChar(STRING_ELT(inBam, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", Rf_translateChar(STRING_ELT(inBam, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("SAM/BAM header missing or empty file: '%s'", Rf_translateChar(STRING_ELT(inBam, 0)));
    }

    FILE * fout = NULL;
    fout = fopen(Rf_translateChar(STRING_ELT(outFile, 0)), "wb");
    if(fout == NULL)
        Rf_error("Error in opening the output file %s", Rf_translateChar(STRING_ELT(outFile, 0)));
                        
    // execute
    int res;
    if(Rf_length(outFile) == 2){
        // paired-end
        FILE * fout_mate = NULL;
        fout_mate = fopen(Rf_translateChar(STRING_ELT(outFile, 1)), "wb");
        if(fout == NULL)
            Rf_error("Error in opening the output file %s", Rf_translateChar(STRING_ELT(outFile, 1)));
        res = _extract_unmapped_paired_reads(fin, fout, fout_mate, Rf_asLogical(fastq), Rf_asLogical(rcRead2));
        fclose(fout_mate);
    }else{
        // single-end
        res = _extract_unmapped_single_reads(fin, fout, Rf_asLogical(fastq));
    }

    // check res
    if(res < 0)
	Rf_error("Error while extracting unmapped reads (return value: %d)", res);
    
    // clean up
    samclose(fin);
    fclose(fout);

    return outFile;
}

