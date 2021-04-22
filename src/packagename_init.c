#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void fnFilterAlphaC(void *, void *, void *, void *, void *);
extern void fnFilterWC(void *, void *, void *, void *, void *);
extern void fnGarchC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void fnQDiffC(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
      {"fnFilterAlphaC", (DL_FUNC) &fnFilterAlphaC, 5},
      {"fnFilterWC",     (DL_FUNC) &fnFilterWC,     5},
      {"fnGarchC",       (DL_FUNC) &fnGarchC,       8},
      {"fnQDiffC",       (DL_FUNC) &fnQDiffC,       5},
      {NULL, NULL, 0}
};

void R_init_bayesGARCH(DllInfo *dll)
{
      R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
      R_useDynamicSymbols(dll, FALSE);
}