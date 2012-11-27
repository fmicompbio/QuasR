#include <stdio.h>
#include <math.h>
#include <Rinternals.h>
#include "samtools/sam.h"


void _reverse(char *buf, int len);
void _complement(char *buf, int len);
samfile_t * _bam_tryopen(const char *filename, const char *filemode, void *aux);
SEXP _getListElement(SEXP list, const char *str);
