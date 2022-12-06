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

extern load_stack_incr_t *lsi;
extern load_stack_decr_t *lsd;
extern pomp_fun_handler_t *pfh;
extern get_covariate_names_t *gcn;
extern table_lookup_t *tl;
extern make_covariate_table_t *mct;
extern set_pomp_userdata_t *spu;
extern unset_pomp_userdata_t *upu;
extern pomp_onestep_sim *pos;
extern SEXP(*randwalk_perturbation_pomp)(SEXP,SEXP);
extern SEXP (*lookup_in_table_pomp)(SEXP,SEXP);

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
//hippie.c
extern SEXP iabf_computations(SEXP x, SEXP params, SEXP Np, SEXP rw_sd, SEXP predmean, SEXP predvar, SEXP filtmean, SEXP trackancestry, SEXP onepar, SEXP weights);
// fcstsampvar.c
extern SEXP do_fcst_samp_var(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// munit_measure.c
extern SEXP do_munit_measure(SEXP object, SEXP X, SEXP vc, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// vunit_measure.c
extern SEXP do_vunit_measure(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// eunit_measure.c
extern SEXP do_eunit_measure(SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi);
// resample.c
extern void nosort_resamp(int nw, double *w, int np, int *p, int offset);
extern SEXP systematic_resampling(SEXP weights);

#endif
