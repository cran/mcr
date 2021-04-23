#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void PaBaLargeData(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP calcAngleMat(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"PaBaLargeData", (DL_FUNC) &PaBaLargeData, 10},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"calcAngleMat", (DL_FUNC) &calcAngleMat, 3},
    {NULL, NULL, 0}
};

void R_init_mcr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
