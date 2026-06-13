#include <stddef.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

extern SEXP flounder_rust_capabilities(void);
extern SEXP flounder_pod5_find(SEXP path);

static const R_CallMethodDef CallEntries[] = {
    {"flounder_rust_capabilities", (DL_FUNC) &flounder_rust_capabilities, 0},
    {"flounder_pod5_find", (DL_FUNC) &flounder_pod5_find, 1},
    {NULL, NULL, 0}
};

attribute_visible void R_init_floundeR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
