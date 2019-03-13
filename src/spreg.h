/* Copyright 2019 by Roger S. Bivand. */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>
/* #include <R_ext/Linpack.h> */
#include <R_ext/Utils.h>
#define ROFFSET 1

SEXP lmin21(SEXP nb, SEXP y, SEXP cy, SEXP card);
SEXP lmin22(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta);
SEXP lmin23(SEXP nb, SEXP y, SEXP cy, SEXP card, SEXP beta, SEXP tol);
SEXP lmin3(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP beta, SEXP tol);
SEXP lmin3S(SEXP nb, SEXP ev1, SEXP ev1_lag, SEXP n_nei, SEXP card, SEXP beta, SEXP tol);

