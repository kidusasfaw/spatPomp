// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>
#include <string.h>

#include "spatPomp_defines.h"
#include "pomp.h"

static R_INLINE SEXP ret_array (int n, int nreps, int ntimes, SEXP names) {
  int dim[3] = {n, nreps, ntimes};
  const char *dimnm[3] = {"variable", "rep", "time"};
  SEXP Y;

  PROTECT(Y = makearray(3,dim));
  setrownames(Y,names,3);
  fixdimnames(Y,dimnm,3);

  UNPROTECT(1);
  return Y;
}

SEXP do_runit_measure (SEXP object, SEXP x, SEXP times, SEXP units, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nvars, npars, ncovars, nreps, nrepsx, nrepsp;
  int nobs = 0;
  SEXP Snames, Pnames, Cnames, Onames = R_NilValue;
  SEXP fn, args;
  SEXP pompfun;
  SEXP Y = R_NilValue;
  int *dim;
  lookup_table_t covariate_table;
  SEXP cvec;
  double *cov;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1)
    errorcall(R_NilValue,"length('times') = 0, no work to do.");

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1];

  if (ntimes != dim[2])
    errorcall(R_NilValue,"length of 'times' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[1];

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    errorcall(R_NilValue,"larger number of replicates is not a multiple of smaller.");

  PROTECT(pompfun = GET_SLOT(object,install("runit_measure"))); nprotect++;
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = (*gcn)(GET_SLOT(object,install("covar")))); nprotect++;
  PROTECT(Onames = GET_SLOT(pompfun,install("obsnames"))); nprotect++;

  // set up the covariate table
  covariate_table = (*mct)(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  // extract the user-defined function
  PROTECT(fn = (*pfh)(pompfun,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {

  case Rfun: {
  }

    break;

  case native: case regNative: {
    double *yt = 0, *xp, *pp;
    int *unit = INTEGER(units);
    double *time = REAL(times), *xs = REAL(x), *ps = REAL(params);
    int *oidx, *sidx, *pidx, *cidx;
    spatPomp_unit_measure_model_simulator *ff = NULL;
    int j, k;

    nobs = LENGTH(Onames);
    // extract observable, state, parameter covariate indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfun,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    // create array to store results
    PROTECT(Y = ret_array(nobs,nreps,ntimes,Onames)); nprotect++;
    yt = REAL(Y);

    GetRNGstate();

    for (k = 0; k < ntimes; k++, time++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      // interpolate the covar functions for the covariates
      (*tl)(&covariate_table,*time,cov);

      for (j = 0; j < nreps; j++, yt += nobs) { // loop over replicates

        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];

        (*ff)(yt,xp,pp,oidx,sidx,pidx,cidx,ncovars,cov,*time,*unit);

      }
    }

    PutRNGstate();

  }

    break;

  default: {
    nobs = LENGTH(Onames);
    int dim[3] = {nobs, nreps, ntimes};
    const char *dimnm[3] = {"variable","rep","time"};
    double *yt = 0;
    int i, n = nobs*nreps*ntimes;

    PROTECT(Y = makearray(3,dim)); nprotect++;
    setrownames(Y,Onames,3);
    fixdimnames(Y,dimnm,3);

    for (i = 0, yt = REAL(Y); i < n; i++, yt++) *yt = R_NaReal;

    warningcall(R_NilValue,"'runit_measure' unspecified: NAs generated.");
  }

  }

  UNPROTECT(nprotect);
  return Y;
}
