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

static const R_CallMethodDef CallEntries[] = {
    {"flounder_rust_capabilities", (DL_FUNC) &flounder_rust_capabilities, 0},
    {"flounder_pod5_find", (DL_FUNC) &flounder_pod5_find, 1},
    {"flounder_pod5_verify", (DL_FUNC) &flounder_pod5_verify, 1},
    {"flounder_pod5_file_info", (DL_FUNC) &flounder_pod5_file_info, 1},
    {"flounder_pod5_folder_info", (DL_FUNC) &flounder_pod5_folder_info, 1},
    {"flounder_pod5_manifest", (DL_FUNC) &flounder_pod5_manifest, 1},
    {"flounder_pod5_compare", (DL_FUNC) &flounder_pod5_compare, 2},
    {NULL, NULL, 0}
};

attribute_visible void R_init_floundeR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
