#include <Rdefines.h>
#include "samtools/sam.h"

SEXP nucleotide_alignment_frequencies(SEXP bamfile, SEXP refSequence, SEXP refChr, SEXP refStart, SEXP mmDist, SEXP chunkSize);
