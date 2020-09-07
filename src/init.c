#include <R_ext/Rdynload.h>
#include "spatPomp_defines.h"

static const R_CallMethodDef callMethods[] = {
  {"do_dunit_measure", (DL_FUNC) &do_dunit_measure, 8},
  {"do_runit_measure", (DL_FUNC) &do_runit_measure, 6},
  {"abf_computations", (DL_FUNC) &abf_computations, 5},
  {"abfir_resample", (DL_FUNC) &abfir_resample, 5},
  {"girf_computations", (DL_FUNC) &girf_computations, 9},
  {"hippie_computations", (DL_FUNC) &hippie_computations, 10},
  {"bpfilter_computations", (DL_FUNC) &bpfilter_computations, 5},
  {"do_fcst_samp_var", (DL_FUNC) &do_fcst_samp_var, 6},
  {"do_v_to_theta", (DL_FUNC) &do_v_to_theta, 7},
  {"do_theta_to_v", (DL_FUNC) &do_theta_to_v, 6},
  {"do_theta_to_e", (DL_FUNC) &do_theta_to_e, 6},
  {NULL, NULL, 0}
};

void R_init_spatPomp (DllInfo *info) {
  // Register routines

  lsi = (load_stack_incr_t *) R_GetCCallable("pomp", "load_stack_incr");
  //sp_load_stack_incr = (psp_load_stack_incr)R_GetCCallable("pomp", "load_stack_incr");


  lsd = (load_stack_decr_t *) R_GetCCallable("pomp", "load_stack_decr");
  //sp_load_stack_decr = (psp_load_stack_decr)R_GetCCallable("pomp", "load_stack_decr");


  pfh = (pomp_fun_handler_t *) R_GetCCallable("pomp", "pomp_fun_handler");
  //sp_pomp_fun_handler = (psp_pomp_fun_handler)R_GetCCallable("pomp", "pomp_fun_handler");

  gcn = (get_covariate_names_t *) R_GetCCallable("pomp", "get_covariate_names");

  //DO I NEED THIS LINE BELOW??
  //lit = (lookup_in_table_t *) R_GetCCallable("pomp", "lookup_in_table");
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
