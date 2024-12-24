/* Copyright 2019-2024 by Roger S. Bivand. */

#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#ifndef R_NO_REMAP
# define R_NO_REMAP
#endif
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <R_ext/Applic.h>
/* #include <R_ext/Linpack.h> */
#include <R_ext/Utils.h>
#define ROFFSET 1

SEXP lmin21(SEXP nb, SEXP y, SEXP cy, SEXP card);
SEXP lmin22(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta);
SEXP lmin23(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta, SEXP tol);
SEXP lmin3(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP beta, SEXP tol);
SEXP lmin3S(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP card, SEXP beta, SEXP tol);
SEXP mom_calc_int2(SEXP is, SEXP m, SEXP nb, SEXP weights, SEXP card);
SEXP listw2dsT(SEXP nbs, SEXP wts, SEXP card, SEXP ncard2);
SEXP listw2dgR(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);

SEXP opt_error_free(SEXP ptr); 
SEXP hess_error_free(SEXP ptr); 
SEXP hess_lag_free(SEXP ptr); 
SEXP opt_error_init(void); /* -Wstrict-prototypes */
SEXP hess_error_init(void); /* -Wstrict-prototypes */
SEXP hess_lag_init(void); /* -Wstrict-prototypes */
SEXP R_ml_sse_env(SEXP env, SEXP coef); 
SEXP R_ml1_sse_env(SEXP env, SEXP lambda, SEXP beta); 
SEXP R_ml2_sse_env(SEXP env, SEXP rho, SEXP beta); 
SEXP mom_calc_int2(SEXP is, SEXP m, SEXP nb, SEXP weights, SEXP card); 

void opt_error_set(SEXP env); 
void hess_error_set(SEXP env); 
void hess_lag_set(SEXP env); 


