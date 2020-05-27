// -*- C++ -*-

#include "spatPomp_defines.h"
#include <Rdefines.h>

SEXP bpfilter_computations (SEXP x, SEXP params, SEXP Np,
			   SEXP trackancestry,
			   SEXP resamp_weights)
{
  int nprotect = 0;
  SEXP anc = R_NilValue;
  SEXP newstates = R_NilValue;
  SEXP retval, retvalnames;
  const char *dimnm[2] = {"variable","rep"};
  double *xw = 0, *xx = 0;
  int *xanc = 0;
  SEXP dimX, dimP, Xnames, Pnames;
  int *dim, np;
  int nvars, npars = 0, nreps;
  int do_ta;
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

  np = *(INTEGER(AS_INTEGER(Np)));

  do_ta = *(LOGICAL(AS_LOGICAL(trackancestry))); // track ancestry?

  xw = REAL(resamp_weights);

  if (do_ta) {
    PROTECT(anc = NEW_INTEGER(np)); nprotect++;
    xanc = INTEGER(anc);
  }
  GetRNGstate();

  int xdim[2];
  int sample[np];
  double *ss = 0, *st = 0;

  // create storage for new states
  xdim[0] = nvars; xdim[1] = np;
  PROTECT(newstates = makearray(2,xdim)); nprotect++;
  setrownames(newstates,Xnames,2);
  fixdimnames(newstates,dimnm,2);
  ss = REAL(x);
  st = REAL(newstates);


  // resample
  nosort_resamp(nreps,REAL(resamp_weights),np,sample,0);
  for (k = 0; k < np; k++) { // copy the particles
    for (j = 0, xx = ss+nvars*sample[k]; j < nvars; j++, st++, xx++) *st = *xx;
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
  SET_ELEMENT(retval,1,params);

  if (do_ta) {
    SET_ELEMENT(retval,2,anc);
  }

  UNPROTECT(nprotect);
  return(retval);
}
