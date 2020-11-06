// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "spatPomp_defines.h"
#include "pomp.h"

static R_INLINE SEXP ret_array (int nreps, int ntimes) {
  int dim[2] = {nreps, ntimes};
  const char *dimnm[2] = {"rep","time"};
  SEXP F;
  PROTECT(F = makearray(2,dim));
  fixdimnames(F,dimnm,2);
  UNPROTECT(1);
  return F;
}

SEXP do_dunit_measure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP units, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nvars, npars, ncovars, nreps, nrepsx, nrepsp, nobs;
  SEXP Snames, Pnames, Cnames, Onames;
  SEXP cvec, pompfun;
  SEXP fn, args;
  SEXP F;
  int *dim;
  lookup_table_t covariate_table;
  double *cov;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1) errorcall(R_NilValue,"length('times') = 0, no work to do.");

  PROTECT(y = as_matrix(y)); nprotect++;
  dim = INTEGER(GET_DIM(y));
  nobs = dim[0];

  if (ntimes != dim[1]) errorcall(R_NilValue,"length of 'times' and 2nd dimension of 'y' do not agree.");

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

  PROTECT(Onames = GET_ROWNAMES(GET_DIMNAMES(y))); nprotect++;
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = (*gcn)(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = (*mct)(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dunit_measure"))); nprotect++;
  PROTECT(fn = (*pfh)(pompfun,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // create array to store results
  PROTECT(F = ret_array(nreps,ntimes)); nprotect++;

  switch (mode) {

  case Rfun: {
  }

    break;

  case native: case regNative: {
    int *oidx, *sidx, *pidx, *cidx;
    int give_log;
    spatPomp_unit_measure_model_density *ff = NULL;
    double *yp = REAL(y), *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    int *unit = INTEGER(units);
    double *ft = REAL(F);
    double *xp, *pp;
    int j, k;

    // extract state, parameter, covariate, observable indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfun,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    give_log = *(INTEGER(AS_INTEGER(log)));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    (*spu)(args);

    for (k = 0; k < ntimes; k++, time++, yp += nobs) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      // interpolate the covar functions for the covariates
      (*tl)(&covariate_table,*time,cov);

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];

        (*ff)(ft,yp,xp,pp,give_log,oidx,sidx,pidx,cidx,ncovars,cov,*time,*unit);
      }
    }

    (*upu)();

  }

    break;

  default: {
    double *ft = REAL(F);
    int j, k;

    for (k = 0; k < ntimes; k++) { // loop over times
      for (j = 0; j < nreps; j++, ft++) { // loop over replicates
        *ft = R_NaReal;
      }
    }

    warningcall(R_NilValue,"'dunit_measure' unspecified: likelihood undefined.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}



