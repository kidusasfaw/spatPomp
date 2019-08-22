// -*- C++ -*-
// resamples states and guide function after checking weights for filtering failure.

#include "spatPomp_defines.h"
#include <Rdefines.h>

SEXP asifir_resample (SEXP x, SEXP Np, SEXP weights, SEXP gps, SEXP tol) {
  int nprotect = 0, all_fail = 0, j, k;
  SEXP newstates = R_NilValue, gfs = R_NilValue, dimX, newdim, Xnames;
  SEXP fail, retval, retvalnames;
  const char *dimnm[2] = {"variable","rep"};
  int *dim, np, nvars, nreps, nlost;
  double toler;

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  np = *(INTEGER(AS_INTEGER(Np))); 

  // identify lost particles (weights less than toler)
  double *xw = REAL(weights);
  toler = *(REAL(tol));		
  for (k = 0, nlost = 0; k < nreps; k++) { 
    if (xw[k] < toler) nlost++;
  }
  if (nlost >= nreps) all_fail = 1; // all particles are lost
  PROTECT(fail = NEW_LOGICAL(1)); nprotect++;  
  *(LOGICAL(fail)) = all_fail;

  GetRNGstate();

  if (!all_fail) { // resample the particles unless we have filtering failure
    int xdim[2], sample[np];
    double *g, *xx, *xf, *xp, *gp, *gf;

    // create storage for resampled states and guide function 
    xdim[0] = nvars; xdim[1] = np;
    PROTECT(newstates = makearray(2,xdim)); nprotect++;
    setrownames(newstates,Xnames,2);
    fixdimnames(newstates,dimnm,2);
    PROTECT(gfs = NEW_NUMERIC(np)); nprotect++;
    xp = REAL(x); 
    xf = REAL(newstates);
    gp = REAL(gps);
    gf = REAL(gfs);

    // resample
    nosort_resamp(nreps,REAL(weights),np,sample,0);
    for (k = 0, g = gp + sample[k]; k < np; k++, gf++) { 
      g = gp + sample[k];
      *gf = *g;
      for (j = 0, xx = xp+nvars*sample[k]; j < nvars; j++, xf++, xx++)
        *xf = *xx;
    }

  } else { // don't resample: just drop 3rd dimension in x prior to return
    PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
    dim = INTEGER(newdim);
    dim[0] = nvars; dim[1] = nreps;
    SET_DIM(x,newdim);
    setrownames(x,Xnames,2);
    fixdimnames(x,dimnm,2);
  }

  PutRNGstate();

  PROTECT(retval = NEW_LIST(3)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(3)); nprotect++;

  SET_STRING_ELT(retvalnames,0,mkChar("fail"));
  SET_ELEMENT(retval,0,fail);
  SET_STRING_ELT(retvalnames,1,mkChar("states"));
  SET_ELEMENT(retval,1, all_fail ? x : newstates);
  SET_STRING_ELT(retvalnames,2,mkChar("filterguides"));
  SET_ELEMENT(retval,2, all_fail ? gps : gfs);
  SET_NAMES(retval,retvalnames);

  UNPROTECT(nprotect);
  return(retval);
}
