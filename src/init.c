#include <R_ext/Rdynload.h>
#include "spatPomp_defines.h"

SEXP randwalk_perturbation_spatPomp(SEXP params, SEXP rw_sd) {
  return(randwalk_perturbation_pomp(params, rw_sd));
}

SEXP lookup_in_table_spatPomp(SEXP covar, SEXP t) {
  return(lookup_in_table_pomp(covar, t));
}

static const R_CallMethodDef callMethods[] = {
  {"randwalk_perturbation_spatPomp", (DL_FUNC) &randwalk_perturbation_spatPomp,2},
  {"lookup_in_table_spatPomp", (DL_FUNC) &lookup_in_table_spatPomp,2},
  {"do_dunit_measure", (DL_FUNC) &do_dunit_measure, 8},
  {"do_runit_measure", (DL_FUNC) &do_runit_measure, 6},
  {"abf_computations", (DL_FUNC) &abf_computations, 5},
  {"abfir_resample", (DL_FUNC) &abfir_resample, 5},
  {"girf_computations", (DL_FUNC) &girf_computations, 9},
  {"bpfilter_computations", (DL_FUNC) &bpfilter_computations, 6},
  {"do_fcst_samp_var", (DL_FUNC) &do_fcst_samp_var, 6},
  {"do_v_to_theta", (DL_FUNC) &do_v_to_theta, 7},
  {"do_theta_to_v", (DL_FUNC) &do_theta_to_v, 6},
  {"do_theta_to_e", (DL_FUNC) &do_theta_to_e, 6},
  {NULL, NULL, 0}
};

SEXP (*randwalk_perturbation_pomp)(SEXP,SEXP);
SEXP (*lookup_in_table_pomp)(SEXP,SEXP);
SEXP (*gcn)(SEXP);
SEXP (*lsd)(SEXP);
SEXP (*lsi)(SEXP);
make_covariate_table_t * mct;
SEXP (*pfh)(SEXP, SEXP, pompfunmode *, SEXP, SEXP, SEXP, SEXP);
void (*spu)(SEXP);
void (*upu)(void);
void (*tl)(lookup_table_t *, double, double *);


void R_init_spatPomp (DllInfo *info) {
  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
  randwalk_perturbation_pomp = (SEXP(*) (SEXP, SEXP)) R_GetCCallable("pomp","randwalk_perturbation");
  lookup_in_table_pomp = (SEXP(*) (SEXP, SEXP)) R_GetCCallable("pomp","lookup_in_table");
  lsi = (load_stack_incr_t *) R_GetCCallable("pomp", "load_stack_incr");
  lsd = (load_stack_decr_t *) R_GetCCallable("pomp", "load_stack_decr");
  pfh = (pomp_fun_handler_t *) R_GetCCallable("pomp", "pomp_fun_handler");
  gcn = (get_covariate_names_t *) R_GetCCallable("pomp", "get_covariate_names");
  tl = (table_lookup_t *) R_GetCCallable("pomp", "table_lookup");
  mct = (make_covariate_table_t *) R_GetCCallable("pomp", "make_covariate_table");
  spu = (set_pomp_userdata_t *) R_GetCCallable("pomp", "set_pomp_userdata");
  upu = (unset_pomp_userdata_t *) R_GetCCallable("pomp", "unset_pomp_userdata");
}
