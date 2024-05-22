#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <R_ext/Rdynload.h>
#include "spatPomp_defines.h"
#include "pomp.h"

static R_INLINE SEXP ret_array (int nunits, int nreps, int ntimes) {
  int dim[3] = {nunits, nreps, ntimes};
  const char *dimnm[3] = {"unit","rep","time"};
  SEXP F;
  PROTECT(F = makearray(3,dim));
  fixdimnames(F,dimnm,3);
  UNPROTECT(1);
  return F;
}

SEXP do_eunit_measure (SEXP object, SEXP X, SEXP Np, SEXP times, SEXP params, SEXP gnsi){
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nunits, nvars, npars, ncovars, nreps, nrepsx, nrepsp;
  SEXP Snames, Pnames, Cnames, Onames;
  SEXP cvec, pompfunthetatoe;
  SEXP fnthetatoe, args;
  SEXP F;
  SEXP x;
  SEXP unitnames;
  int *dim;
  lookup_table_t covariate_table;
  double *cov;
  PROTECT(Np = AS_INTEGER(Np)); nprotect++;
  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1) errorcall(R_NilValue,"length('times') = 0, no work to do.");


  PROTECT(x = as_state_array(X)); nprotect++;
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


  // extract the user-defined function
  PROTECT(pompfunthetatoe = GET_SLOT(object,install("eunit_measure"))); nprotect++;


  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = (*gcn)(GET_SLOT(object,install("covar")))); nprotect++;
  PROTECT(Onames = GET_SLOT(pompfunthetatoe,install("obsnames"))); nprotect++;

  PROTECT(fnthetatoe = (*pfh)(pompfunthetatoe,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // set up the covariate table
  covariate_table = (*mct)(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  PROTECT(unitnames = GET_SLOT(object,install("unit_names"))); nprotect++;
  nunits = length(unitnames);

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // create array to store results
  PROTECT(F = ret_array(nunits, nreps, ntimes)); nprotect++;
  switch (mode) {

  case native: case regNative: {
    int *oidx, *sidx, *pidx, *cidx;
    spatPomp_unit_measure_var *ffthetatoe = NULL;
    double *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    double *ft = REAL(F);
    double *xp, *pp;
    int i, j, k;

    // extract state, parameter, covariate, observable indices
    sidx = INTEGER(GET_SLOT(pompfunthetatoe,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfunthetatoe,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfunthetatoe,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfunthetatoe,install("covarindex")));

    // address of native routine
    *((void **) (&ffthetatoe)) = R_ExternalPtrAddr(fnthetatoe);

    for (k = 0; k < ntimes; k++, time++) { // loop over times
      // interpolate the covar functions for the covariates
      (*tl)(&covariate_table,*time,cov);
      R_CheckUserInterrupt();	// check for user interrupt
      for (j = 0; j < nreps; j++) { // loop over replicates
        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];
        for(i = 0; i < nunits; i++, ft++){
          (*ffthetatoe)(ft,xp,pp,oidx,sidx,pidx,cidx,ncovars,cov,*time,i);
        }
      }
    }
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
    warningcall(R_NilValue,"'eunit_measure' unspecified.");
  }
  }
  // create array to store variances for each combination of unit, particle and lookahead
  UNPROTECT(nprotect);
  return F;
}
