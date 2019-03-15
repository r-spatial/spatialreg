/* Copyright 2019 by Roger S. Bivand. */

#include "spreg.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*static const R_CMethodDef CEntries[]  = {
    {NULL, NULL, 0}
};*/


static R_CallMethodDef CallEntries[] = {
    {"lmin21", (DL_FUNC) &lmin21, 4},
    {"lmin22", (DL_FUNC) &lmin22, 5},
    {"lmin23", (DL_FUNC) &lmin23, 6},
    {"lmin3", (DL_FUNC) &lmin3, 6},
    {"lmin3S", (DL_FUNC) &lmin3S, 7},
    {"mom_calc_int2", (DL_FUNC) &mom_calc_int2, 5},
    {NULL, NULL, 0}
};


void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_spdep(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

}



