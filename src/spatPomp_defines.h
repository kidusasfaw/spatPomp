// -*- C++ -*-

#ifndef _SPATPOMP_DEFINES_H_
#define _SPATPOMP_DEFINES_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "pomp_defines.h"

typedef void spatPomp_unit_measure_model_density (double *lik, const double *y, const double *x, const double *p, int give_log,
                                         const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
                                         int ncovars, const double *covars, double t, int u);

typedef void spatPomp_unit_measure_model_simulator (double *y, const double *x, const double *p,
                                                  const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
                                                  int ncovars, const double *covars, double t, int u);
typedef void spatPomp_unit_measure_mean (double *y, const double *x, const double *p,
                                                    const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
                                                    int ncovars, const double *covars, double t, int u);
typedef void spatPomp_unit_mmeasure (double *y, const double *x, const double *p, const double *v,
                                     const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
                                     int ncovars, const double *covars, double t, int u);
typedef void spatPomp_unit_measure_var (double *v, const double *x, const double *p,
                                        const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
                                        int ncovars, const double *covars, double t, int u);

load_stack_incr_t *lsi;
load_stack_decr_t *lsd;
pomp_fun_handler_t *pfh;
get_covariate_names_t *gcn;
table_lookup_t *tl;
make_covariate_table_t *mct;
set_pomp_userdata_t *spu;
unset_pomp_userdata_t *upu;
pomp_onestep_sim *pos;

// dunit_measure.c
extern SEXP do_dunit_measure(SEXP object, SEXP y, SEXP x, SEXP times, SEXP units, SEXP params, SEXP log, SEXP gnsi);
// runit_measure.c
extern SEXP do_runit_measure(SEXP object, SEXP x, SEXP times, SEXP units, SEXP params, SEXP gnsi);
// abf.c
extern SEXP abf_computations(SEXP x, SEXP params, SEXP Np, SEXP trackancestry, SEXP weights);
// abfir.c
extern SEXP abfir_resample(SEXP x, SEXP Np, SEXP weights, SEXP gps, SEXP tol);
// bpfilter.c
extern SEXP bpfilter_computations(SEXP x, SEXP params, SEXP Np, SEXP trackancestry, SEXP doparRS, SEXP resamp_weights);
// girf.c
extern SEXP girf_computations(SEXP x, SEXP params, SEXP Np, SEXP trackancestry, SEXP doparRS, SEXP weights, SEXP gps, SEXP fsv, SEXP tol);
// fcstsampvar.c
extern SEXP do_fcst_samp_var(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// v_to_theta.c
extern SEXP do_v_to_theta(SEXP object, SEXP X, SEXP vc, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// theta_to_v.c
extern SEXP do_theta_to_v(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// theta_to_e.c
extern SEXP do_theta_to_e(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// resample.c
extern void nosort_resamp(int nw, double *w, int np, int *p, int offset);
extern SEXP systematic_resampling(SEXP weights);

#endif
