// prevent remapping e.g. Ralloc, which causes conflicts under windows
#define STRICT_R_HEADERS 1

#include "merge_reorder_sam.h"
#include "quantify_methylation.h"
#include "count_junctions.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <R_ext/Rdynload.h>
#include "split_sam_chr.h"
#include "cat_bam.h"
#include "idxstats_bam.h"
#include "convert_reads_id_bis_rc.h"
#include "extract_unmapped_reads.h"
#include "remove_unmapped_from_sam.h"
#include "count_alignments.h"
#include "profile_alignments.h"
#include "count_alignments_subregions.h"
#include "export_wig.h"
#include "nucleotide_alignment_frequencies.h"

static const R_CallMethodDef callMethods[] = {
    /* split_sam_chr.c */
    {"splitSamChr", (DL_FUNC) &split_sam_chr, 2},
    /* cat_bam.c */
    {"catBam", (DL_FUNC) &cat_bam, 2},
    /* idxstats_bam.c */
    {"idxstatsBam", (DL_FUNC) &idxstats_bam, 1},
    /* merge_reorder_sam.c */
    {"mergeReorderSam", (DL_FUNC) &merge_reorder_sam, 4},
    /* convert_bisulfite_reads.c */
    {"convertReadsIdBisRc", (DL_FUNC) &convert_reads_id_bis_rc, 4},
    /* extract_unmapped_reads.c */
    {"extractUnmappedReads", (DL_FUNC) &extract_unmapped_reads, 4},
    /* remove_unmapped_from_sam.c */
    {"removeUnmappedFromSamAndConvertToBam", (DL_FUNC) &remove_unmapped_from_sam_and_convert_to_bam, 2},
    /* count_alignments.c */
    {"countAlignmentsNonAllelic", (DL_FUNC) &count_alignments_non_allelic, 14},
    {"countAlignmentsAllelic", (DL_FUNC) &count_alignments_allelic, 14},
    /* count_junctions.cpp */
    {"countJunctions", (DL_FUNC) &count_junctions, 7},
    /* profile_alignments.c */
    {"profileAlignmentsNonAllelic", (DL_FUNC) &profile_alignments_non_allelic, 19},
    {"profileAlignmentsAllelic", (DL_FUNC) &profile_alignments_allelic, 19},
    /* count_alignments_subregions.c */
    {"countAlignmentsSubregions", (DL_FUNC) &count_alignments_subregions, 10},
    /* quantify_methylation.cpp */
    {"quantifyMethylation", (DL_FUNC) &quantify_methylation, 9},
    {"detectSNVs", (DL_FUNC) &detect_SNVs, 8},
    {"quantifyMethylationAllele", (DL_FUNC) &quantify_methylation_allele, 9},
    {"quantifyMethylationSingleAlignments", (DL_FUNC) &quantify_methylation_singleAlignments, 8},
    /* export_wig.c */
    {"bamfileToWig", (DL_FUNC) &bamfile_to_wig, 15},
    /* nucleotide_alignment_frequencies.c */
    {"nucleotideAlignmentFrequencies", (DL_FUNC) &nucleotide_alignment_frequencies, 6},
    {NULL, NULL, 0}
};


void R_init_QuasR(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_QuasR(DllInfo *info){}

#ifdef __cplusplus
}  /*    extern "C"   */
#endif
