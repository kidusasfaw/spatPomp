// -*- C++ -*-

#ifndef _SPATPOMP3_DEFINES_H_
#define _SPATPOMP3_DEFINES_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "pomp_defines.h"

typedef SEXP(*psp_load_stack_incr)();
typedef SEXP(*psp_load_stack_decr)();
typedef SEXP(*psp_pomp_fun_handler)();

typedef SEXP(*psp_lookup_in_table)();
typedef void(*psp_table_lookup)();
typedef struct lookup_table(*psp_make_covariate_table)();

psp_load_stack_incr sp_load_stack_incr;
psp_load_stack_decr sp_load_stack_decr;
psp_pomp_fun_handler sp_pomp_fun_handler;

psp_lookup_in_table sp_lookup_in_table;
psp_table_lookup sp_table_lookup;
psp_make_covariate_table sp_make_covariate_table;

//dmeasure3.c
extern SEXP do_dmeasure3(SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);


#endif
