/*

bam_cat -- efficiently concatenates bam files

bam_cat can be used to concatenate BAM files. Under special
circumstances, it can be used as an alternative to 'samtools merge' to
concatenate multiple sorted files into a single sorted file. For this
to work each file must be sorted, and the sorted files must be given
as command line arguments in order such that the final read in file i
is less than or equal to the first read in file i+1.

This code is derived from the bam_reheader function in samtools 0.1.8
and modified to perform concatenation by Chris Saunders on behalf of
Illumina.


########## License:

The MIT License

Original SAMtools work copyright (c) 2008-2009 Genome Research Ltd.
Modified SAMtools work copyright (c) 2010 Illumina, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/


/*
makefile:
"""
CC=gcc
CFLAGS+=-g -Wall -O2 -D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -I$(SAMTOOLS_DIR)
LDFLAGS+=-L$(SAMTOOLS_DIR)
LDLIBS+=-lbam -lz

all:bam_cat
"""
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "samtools/bgzf.h"
#include "samtools/bam.h"
#include "cat_bam.h"

#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28

// copied bam_cat from samtools
int _bam_cat(int nfn, char * const *fn, const bam_header_t *h, const char* outbam)
{
    BGZF *fp;
    FILE* fp_file;
    uint8_t *buf;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;
    
    // fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    fp = bgzf_open(outbam, "w");
    if (fp == 0) {
        Rf_error("[%s] ERROR: fail to open output file '%s'.\n", __func__, outbam);
        return 1;
    }
    if (h) bam_header_write(fp, h);
    
    buf = (uint8_t*) malloc(BUF_SIZE);
    for(i = 0; i < nfn; ++i){
        BGZF *in;
        bam_header_t *old;
        int len,j;
        
        in = strcmp(fn[i], "-")? bam_open(fn[i], "r") : bam_dopen(fileno(stdin), "r");
        if (in == 0) {
            Rf_error("[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            return -1;
        }
        if (in->open_mode != 'r') return -1;
        
        old = bam_header_read(in);
		if (h == 0 && i == 0) bam_header_write(fp, old);
        
        if (in->block_offset < in->block_length) {
            bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset);
            bgzf_flush(fp);
        }
        
        j=0;
#ifdef _USE_KNETFILE
        fp_file=fp->x.fpw;
        while ((len = knet_read(in->x.fpr, buf, BUF_SIZE)) > 0) {
#else  
        fp_file=fp->file;
        while (!feof(in->file) && (len = fread(buf, 1, BUF_SIZE, in->file)) > 0) {
#endif
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    Rf_error("[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i]);
                    return -1;
                }
                fwrite(ebuf, 1, len, fp_file);
                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) fwrite(ebuf, 1, es, fp_file);
                len-= es;
                memcpy(ebuf,buf+len,es);
                fwrite(buf, 1, len, fp_file);
            }
            j=1;
        }

        /* check final gzip block */
        {
            const uint8_t gzip1=ebuf[0];
            const uint8_t gzip2=ebuf[1];
            const uint32_t isize=*((uint32_t*)(ebuf+es-4));
            if(((gzip1!=GZIPID1) || (gzip2!=GZIPID2)) || (isize!=0)) {
                Rprintf("[%s] WARNING: Unexpected block structure in file '%s'.", __func__, fn[i]);
                Rf_error(" Possible output corruption.\n");
                fwrite(ebuf, 1, es, fp_file);
            }
        }
        bam_header_destroy(old);
        bgzf_close(in);
    }
    free(buf);
    bgzf_close(fp);
    return 0;
}



SEXP cat_bam(SEXP inbam, SEXP outbam)
{
    bam_header_t *h = 0;
	int ret = 0;
	
	// check parameter inbam
    if(!Rf_isString(inbam)){
        Rf_error("'inbam' must be character()");
    }
    // number of input bamfiles
    int nfin = Rf_length(inbam);
    // convert/copy STRSXP to a const char **
    const char ** fin = (const char **)R_Calloc(nfin, const char *);
    for(int i = 0; i < nfin; i++){
    	fin[i] = Rf_translateChar(STRING_ELT(inbam, i));
    }

    // check parameter outbam
    if(!Rf_isString(outbam) || 1 != Rf_length(outbam)){
        Rf_error("'outbam' must be character(1)");
    }
    const char * fout = Rf_translateChar(STRING_ELT(outbam, 0));

	// execute bam_cat copied from samtools
    ret = _bam_cat(nfin, (char * const *)fin, h, fout);

	R_Free(fin);

	return Rf_ScalarInteger(ret);
}
