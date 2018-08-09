// dear emacs, please treat this as -*- C++ -*-

#include "spatpomp3_defines.h"

SEXP __pomp_userdata;
#define USERDATA  (__pomp_userdata)

void set_pomp_userdata (SEXP userdata) {
  USERDATA = userdata;
}

void unset_pomp_userdata (void) {
  USERDATA = R_NilValue;
}
