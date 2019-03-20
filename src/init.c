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
    {"listw2dsT", (DL_FUNC) &listw2dsT, 4},
    {"listw2dgR", (DL_FUNC) &listw2dgR, 4},
    {"opt_error_free", (DL_FUNC) &opt_error_free, 1},
    {"hess_error_free", (DL_FUNC) &hess_error_free, 1},
    {"hess_lag_free", (DL_FUNC) &hess_lag_free, 1},
    {"opt_error_init", (DL_FUNC) &opt_error_init, 0},
    {"hess_error_init", (DL_FUNC) &hess_error_init, 0},
    {"hess_lag_init", (DL_FUNC) &hess_lag_init, 0},
    {"R_ml_sse_env", (DL_FUNC) &R_ml_sse_env, 2},
    {"R_ml1_sse_env", (DL_FUNC) &R_ml1_sse_env, 3},
    {"R_ml2_sse_env", (DL_FUNC) &R_ml2_sse_env, 3},
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



