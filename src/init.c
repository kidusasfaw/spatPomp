#include <R_ext/Rdynload.h>
#include "spatpomp3_defines.h"

static const R_CallMethodDef callMethods[] = {
  {"do_dmeasure3", (DL_FUNC) &do_dmeasure3, 7},
  {NULL, NULL, 0}
};

void R_init_spatpomp3 (DllInfo *info) {
  // Register routines
  sp_load_stack_incr = (psp_load_stack_incr)R_GetCCallable("pomp", "load_stack_incr");
  sp_load_stack_decr = (psp_load_stack_decr)R_GetCCallable("pomp", "load_stack_decr");
  sp_pomp_fun_handler = (psp_pomp_fun_handler)R_GetCCallable("pomp", "pomp_fun_handler");

  sp_lookup_in_table = (psp_lookup_in_table)R_GetCCallable("pomp", "lookup_in_table");
  sp_table_lookup = (psp_table_lookup)R_GetCCallable("pomp", "table_lookup");
  sp_make_covariate_table = (psp_make_covariate_table)R_GetCCallable("pomp", "make_covariate_table");


  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
