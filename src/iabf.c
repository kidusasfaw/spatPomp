// -*- C++ -*-

#include "spatPomp_defines.h"
#include <Rdefines.h>

// examines resampling weights for filtering failure
// computes log likelihood and effective sample size
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// resampling weights are used in filtering mean computation.
// if length(resamp_weights) == 1, an unweighted average is computed.
// returns all of the above in a named list
SEXP iabf_computations (SEXP x, SEXP params, SEXP Np,
			   SEXP rw_sd,
			   SEXP predmean, SEXP predvar,
			   SEXP filtmean, SEXP trackancestry, SEXP onepar,
			   SEXP resamp_weights)
{
  int nprotect = 0;
  SEXP pm = R_NilValue, pv = R_NilValue, fm = R_NilValue, anc = R_NilValue;
  SEXP rw_names, ess, fail, loglik;
  SEXP newstates = R_NilValue, newparams = R_NilValue;
  SEXP retval, retvalnames;
  const char *dimnm[2] = {"variable","rep"};
  double *xx = 0, *xp = 0;
  int *xanc = 0;
  SEXP dimX, dimP, Xnames, Pnames, pindex;
  int *dim, lv, np;
  int nvars, npars = 0, nrw = 0, nreps;
  int do_rw, do_pm, do_pv, do_fm, do_ta, do_par_resamp;
  //double sum, sumsq, vsq, ws, w, toler;
  int j, k, l;

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  xx = REAL(x);
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;

  PROTECT(dimP = GET_DIM(params)); nprotect++;
  dim = INTEGER(dimP);
  npars = dim[0];
  if (nreps % dim[1] != 0)
    errorcall(R_NilValue,"ncol('states') should be a multiple of ncol('params')"); // # nocov
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  //np = *(INTEGER(AS_INTEGER(Np)));       // number of particles to resample
  np = 1;       // number of particles to resample for iif is just one

  nrw = LENGTH(rw_sd);	     // number of parameters that are variable
  do_rw = nrw > 0;	     // do random walk in parameters?
  if (do_rw) {
    // names of parameters undergoing random walk
    PROTECT(rw_names = GET_NAMES(rw_sd)); nprotect++;
  }

  do_pm = *(LOGICAL(AS_LOGICAL(predmean))); // calculate prediction means?
  do_pv = *(LOGICAL(AS_LOGICAL(predvar)));  // calculate prediction variances?
  do_fm = *(LOGICAL(AS_LOGICAL(filtmean))); // calculate filtering means?
  do_ta = *(LOGICAL(AS_LOGICAL(trackancestry))); // track ancestry?
  do_par_resamp = *(LOGICAL(AS_LOGICAL(onepar))); // are all cols of 'params' the same?
  do_par_resamp = !do_par_resamp || do_rw; // should we do parameter resampling?

  if (do_par_resamp) {
    if (dim[1] != nreps)
      errorcall(R_NilValue,"ncol('states') should be equal to ncol('params')"); // # nocov
  }

  PROTECT(ess = NEW_NUMERIC(1)); nprotect++; // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); nprotect++; // log likelihood
  PROTECT(fail = NEW_LOGICAL(1)); nprotect++;	// particle failure?

  if (do_rw) {
    // indices of parameters undergoing random walk
    PROTECT(pindex = matchnames(Pnames,rw_names,"parameters")); nprotect++;
    xp = REAL(params);
    lv = nvars+nrw;
  } else {
    lv = nvars;
  }

  if (do_pm || do_pv) {
    PROTECT(pm = NEW_NUMERIC(lv)); nprotect++;
  }

  if (do_pv) {
    PROTECT(pv = NEW_NUMERIC(lv)); nprotect++;
  }

  if (do_fm) {
    if (do_rw) {
      PROTECT(fm = NEW_NUMERIC(nvars+npars)); nprotect++;
    } else {
      PROTECT(fm = NEW_NUMERIC(nvars)); nprotect++;
    }
  }

  if (do_ta) {
    PROTECT(anc = NEW_INTEGER(np)); nprotect++;
    xanc = INTEGER(anc);
  }
  GetRNGstate();

  // do the following whether or not there is "filtering failure"
  int xdim[2];
  int copies = *(INTEGER(AS_INTEGER(Np)));       // number of particles to copy after systematic resampling one particle
  int sample[copies];
  double *ss = 0, *st = 0, *ps = 0, *pt = 0;

  // create storage for new states
  xdim[0] = nvars; xdim[1] = copies;
  PROTECT(newstates = makearray(2,xdim)); nprotect++;
  setrownames(newstates,Xnames,2);
  fixdimnames(newstates,dimnm,2);
  ss = REAL(x);
  st = REAL(newstates);

  // create storage for new parameters
  if (do_par_resamp) {
    xdim[0] = npars; xdim[1] = copies;
    PROTECT(newparams = makearray(2,xdim)); nprotect++;
    setrownames(newparams,Pnames,2);
    fixdimnames(newparams,dimnm,2);
    ps = REAL(params);
    pt = REAL(newparams);
  }

  // resample
  nosort_resamp(nreps,REAL(resamp_weights),np,sample,0);
  for(l = 1; l < copies; l++){ // replicate the resampled particle number `copies` number of times
    sample[l] = sample[l-1];
  }

  for (k = 0; k < copies; k++) { // copy the particles
    for (j = 0, xx = ss+nvars*sample[k]; j < nvars; j++, st++, xx++) *st = *xx;
    if (do_par_resamp) {
     for (j = 0, xp = ps+npars*sample[k]; j < npars; j++, pt++, xp++) *pt = *xp;
    }
    if (do_ta) xanc[k] = sample[k]+1;
  }

  PutRNGstate();

  PROTECT(retval = NEW_LIST(3)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(3)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("states"));
  SET_STRING_ELT(retvalnames,1,mkChar("params"));
  SET_STRING_ELT(retvalnames,2,mkChar("ancestry"));
  SET_NAMES(retval,retvalnames);
  SET_ELEMENT(retval,0,newstates);
  SET_ELEMENT(retval,1,newparams);
  if (do_ta) {
    SET_ELEMENT(retval,2,anc);
  }

  UNPROTECT(nprotect);
  return(retval);
}
