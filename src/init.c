#include "pomp_defines.h"
#include <R_ext/Rdynload.h>

extern SEXP do_dmeasure3(SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);

static const R_CallMethodDef callMethods[] = {
  {"do_dmeasure3", (DL_FUNC) &do_dmeasure3, 7},
  {NULL, NULL, 0}
};

void R_init_pomp (DllInfo *info) {
  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
