#include <R_ext/Rdynload.h>
#include "quasr.h"

static const R_CallMethodDef callMethods[] = {
    /* weight_alignments.c */
    {".weight_alignments", (DL_FUNC) &weight_alignments, 4},
    {".getAllowedMaxHits", (DL_FUNC) &get_allowed_max_hits, 0},
    /* count_alignments.c */
    {".count_alignments", (DL_FUNC) &count_alignments, 5},
    {".seqname", (DL_FUNC) &seqname, 1},
    {NULL, NULL, 0}
};


void R_init_QuasR(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_QuasR(DllInfo *info){}

