#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "spatPomp_defines.h"
#include "pomp.h"

SEXP do_munit_measure(SEXP object, SEXP X, SEXP vc, SEXP Np, SEXP times, SEXP params, SEXP gnsi){
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nunits, nvars, npars, ncovars, nparticles, nreps, nrepsx, nrepsp;
  SEXP Snames, Pnames, Cnames, Onames;
  SEXP cvec, pompfun;
  SEXP fn, args;
  SEXP mparams;
  SEXP x;
  SEXP unitnames;
  int *dim;
  lookup_table_t covariate_table;
  double *cov;
  PROTECT(Np = AS_INTEGER(Np)); nprotect++;
  nparticles = *INTEGER(Np);
  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1) errorcall(R_NilValue,"length('times') = 0, no work to do.");


  PROTECT(x = as_state_array(X)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1];

  if (ntimes != dim[2])
    errorcall(R_NilValue,"length of 'times' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = duplicate(params)); nprotect++;
  PROTECT(mparams = duplicate(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[2];

  nreps = (nrepsp > nparticles) ? nrepsp : nparticles;

  if ((nreps % nrepsp != 0))
    errorcall(R_NilValue,"larger number of replicates is not a multiple of smaller.");


  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("munit_measure"))); nprotect++;

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = (*gcn)(GET_SLOT(object,install("covar")))); nprotect++;
  PROTECT(Onames = GET_SLOT(pompfun,install("obsnames"))); nprotect++;

  PROTECT(fn = (*pfh)(pompfun,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // set up the covariate table
  covariate_table = (*mct)(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  PROTECT(unitnames = GET_SLOT(object,install("unit_names"))); nprotect++;
  nunits = length(unitnames);

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // create array to store results
  // PROTECT(F = ret_array(npars, nunits, nreps, ntimes)); nprotect++;
  switch (mode) {

  case native: case regNative: {
    int *oidx, *sidx, *pidx, *cidx;
    spatPomp_unit_mmeasure *ff = NULL;
    double *xs = REAL(x), *ps = REAL(params), *time = REAL(times), *v = REAL(vc);
    double *ft = REAL(mparams);
    double *xp, *pp;
    int i, j, k;

    // extract state, parameter, covariate, observable indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfun,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    for (k = 0; k < ntimes; k++, time++) { // loop over times
      // interpolate the covar functions for the covariates
      (*tl)(&covariate_table,*time,cov);
      R_CheckUserInterrupt();	// check for user interrupt
      for (j = 0; j < nreps; j++) { // loop over replicates
        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)+nrepsp*k];
        for(i = 0; i < nunits; i++, ft+=npars, v++){
          (*ff)(ft,xp,pp,v,oidx,sidx,pidx,cidx,ncovars,cov,*time,i);
        }
      }

    }

  }

    break;

  default: {
    double *ft = REAL(mparams);
    int j, k;

    for (k = 0; k < ntimes; k++) { // loop over times
      for (j = 0; j < nreps; j++, ft++) { // loop over replicates
        *ft = R_NaReal;
      }
    }

    warningcall(R_NilValue,"'munit_measure' unspecified.");

  }

  }
  // create array to store variances for each combination of unit, particle and lookahead

  UNPROTECT(nprotect);
  return mparams;
}
