#include "htslib/kseq.h"
// include Boolean.h early, will define TRUE/FALSE enum prefent Rdefines.h from defining them as int constants
#include <R_ext/Boolean.h>
#include <Rdefines.h>
#include <zlib.h>
#include <stdio.h>

/* for mingw32 long long storage and printf %llu  */
/* includes <stdint.h> for long long I64u declaration
   and has the printf %I64u format specifier */
#include <inttypes.h>  

SEXP convert_reads_id_bis_rc(SEXP inFile, SEXP outFile, SEXP fromToCharacter, SEXP reverseComplement);


