#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare C routines */
extern SEXP _ZC_Test(SEXP _X, SEXP _Y, SEXP _Param);
extern SEXP _CGZ_Test(SEXP _X, SEXP _Y, SEXP _Param);
extern SEXP _FZW_Test(SEXP _X, SEXP _Param);

extern SEXP RPtest0(SEXP X_, SEXP Y_, SEXP DIM_);
extern SEXP RPtest(SEXP X_, SEXP Y_, SEXP Pk_, SEXP DIM_);

extern SEXP GCtest_(SEXP X_, SEXP Y_, SEXP PSI_, SEXP DIM_);

/* Register routines */
static const R_CallMethodDef CallEntries[] = {
    {"_ZC_Test",  (DL_FUNC) &_ZC_Test,  3},
    {"_CGZ_Test", (DL_FUNC) &_CGZ_Test, 3},
    {"_FZW_Test", (DL_FUNC) &_FZW_Test, 2},

    {"RPtest0",   (DL_FUNC) &RPtest0,   3},
    {"RPtest",    (DL_FUNC) &RPtest,    4},

    {"GCtest_",   (DL_FUNC) &GCtest_,   4},

    {NULL, NULL, 0}
};

void R_init_hdcrtRepro(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}