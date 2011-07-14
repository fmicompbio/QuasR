#include <R_ext/Rdynload.h>
#include "quasr.h"

static const 
R_CallMethodDef callMethods[] = {
  /* count_alignments.c */
	{".countAlignments", (DL_FUNC) &count_alignments, 3},
	{".getAllowedMaxHits", (DL_FUNC) &get_allowed_max_hits, 0},
  {NULL, NULL, 0}
};


void
R_init_QuasR(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void
R_unload_QuasR(DllInfo *info)
{
}

