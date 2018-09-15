#include <R_ext/Rdynload.h>
#include "spatpomp3_defines.h"

static const R_CallMethodDef callMethods[] = {
  {"do_unit_dmeasure", (DL_FUNC) &do_unit_dmeasure, 8},
  {"iif_computations", (DL_FUNC) &iif_computations, 10},
  {"hippie_computations", (DL_FUNC) &hippie_computations, 10},
  {NULL, NULL, 0}
};

void R_init_spatpomp3 (DllInfo *info) {
  // Register routines

  lsi = (load_stack_incr_t *) R_GetCCallable("pomp", "load_stack_incr");
  //sp_load_stack_incr = (psp_load_stack_incr)R_GetCCallable("pomp", "load_stack_incr");


  lsd = (load_stack_decr_t *) R_GetCCallable("pomp", "load_stack_decr");
  //sp_load_stack_decr = (psp_load_stack_decr)R_GetCCallable("pomp", "load_stack_decr");


  pfh = (pomp_fun_handler_t *) R_GetCCallable("pomp", "pomp_fun_handler");
  //sp_pomp_fun_handler = (psp_pomp_fun_handler)R_GetCCallable("pomp", "pomp_fun_handler");


  lit = (lookup_in_table_t *) R_GetCCallable("pomp", "lookup_in_table");
  //sp_lookup_in_table = (psp_lookup_in_table)R_GetCCallable("pomp", "lookup_in_table");


  tl = (table_lookup_t *) R_GetCCallable("pomp", "table_lookup");
  //sp_table_lookup = (psp_table_lookup)R_GetCCallable("pomp", "table_lookup");


  mct = (make_covariate_table_t *) R_GetCCallable("pomp", "make_covariate_table");
  //sp_make_covariate_table = (psp_make_covariate_table)R_GetCCallable("pomp", "make_covariate_table");


  spu = (set_pomp_userdata_t *) R_GetCCallable("pomp", "set_pomp_userdata");


  upu = (unset_pomp_userdata_t *) R_GetCCallable("pomp", "unset_pomp_userdata");


  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
