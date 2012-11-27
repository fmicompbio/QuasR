#include <Rdefines.h>
#include "samtools/sam.h"
#include <stdio.h>

SEXP extract_unmapped_reads(SEXP inBam, SEXP outFile, SEXP fastq, SEXP rcRead2);

// src copy from io_sam.c of Rsamtools
char * _bamseq(const bam1_t * bam, int reverse_comp);
char * _bamqual(const bam1_t * bam, int reverse);

void _write_fastq_seq(FILE * fp, const bam1_t * bam, int reverse_comp);
void _write_fasta_seq(FILE * fp, const bam1_t * bam, int reverse_comp);
int _extract_unmapped_single_reads(samfile_t *fin, FILE *fout, int fastq);
int _extract_unmapped_paired_reads(samfile_t *fin, FILE *fout1, FILE *fout2, int fastq, int rc_read2);
