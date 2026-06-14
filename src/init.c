#include <stddef.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

extern SEXP flounder_rust_capabilities(void);
extern SEXP flounder_pod5_find(SEXP path);
extern SEXP flounder_pod5_verify(SEXP path);
extern SEXP flounder_pod5_file_info(SEXP path);
extern SEXP flounder_pod5_folder_info(SEXP path);
extern SEXP flounder_pod5_manifest(SEXP path);
extern SEXP flounder_pod5_compare(SEXP left, SEXP right);
extern SEXP flounder_pod5_subdivide_plan(SEXP path, SEXP strategy, SEXP files_per_chunk, SEXP seconds_per_chunk, SEXP reads_per_chunk);
extern SEXP flounder_bam_summary(SEXP path, SEXP sample_records, SEXP prefer_index, SEXP include_mapq_hist, SEXP include_flags, SEXP allow_incomplete);
extern SEXP flounder_bam_verify(SEXP path);
extern SEXP flounder_bam_validate(SEXP path, SEXP max_errors, SEXP max_warnings, SEXP header_only, SEXP records, SEXP fail_fast, SEXP include_warnings);
extern SEXP flounder_bam_check_eof(SEXP path);
extern SEXP flounder_bam_check_index(SEXP path, SEXP require, SEXP prefer_csi);
extern SEXP flounder_bam_check_map(SEXP path, SEXP sample_records, SEXP prefer_index);
extern SEXP flounder_bam_check_sort(SEXP path, SEXP sample_records, SEXP strict);
extern SEXP flounder_bam_check_tag(SEXP path, SEXP tag, SEXP sample_records, SEXP full_scan, SEXP require_type, SEXP count_hits);

static const R_CallMethodDef CallEntries[] = {
    {"flounder_rust_capabilities", (DL_FUNC) &flounder_rust_capabilities, 0},
    {"flounder_pod5_find", (DL_FUNC) &flounder_pod5_find, 1},
    {"flounder_pod5_verify", (DL_FUNC) &flounder_pod5_verify, 1},
    {"flounder_pod5_file_info", (DL_FUNC) &flounder_pod5_file_info, 1},
    {"flounder_pod5_folder_info", (DL_FUNC) &flounder_pod5_folder_info, 1},
    {"flounder_pod5_manifest", (DL_FUNC) &flounder_pod5_manifest, 1},
    {"flounder_pod5_compare", (DL_FUNC) &flounder_pod5_compare, 2},
    {"flounder_pod5_subdivide_plan", (DL_FUNC) &flounder_pod5_subdivide_plan, 5},
    {"flounder_bam_summary", (DL_FUNC) &flounder_bam_summary, 6},
    {"flounder_bam_verify", (DL_FUNC) &flounder_bam_verify, 1},
    {"flounder_bam_validate", (DL_FUNC) &flounder_bam_validate, 7},
    {"flounder_bam_check_eof", (DL_FUNC) &flounder_bam_check_eof, 1},
    {"flounder_bam_check_index", (DL_FUNC) &flounder_bam_check_index, 3},
    {"flounder_bam_check_map", (DL_FUNC) &flounder_bam_check_map, 3},
    {"flounder_bam_check_sort", (DL_FUNC) &flounder_bam_check_sort, 3},
    {"flounder_bam_check_tag", (DL_FUNC) &flounder_bam_check_tag, 6},
    {NULL, NULL, 0}
};

attribute_visible void R_init_floundeR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
