# Spatial regression model Jacobian computations

These functions are made available in the package namespace for other
developers, and are not intended for users. They provide a shared
infrastructure for setting up data for Jacobian computation, and then
for caclulating the Jacobian, either exactly or approximately, in
maximum likelihood fitting of spatial regression models. The techniques
used are the exact eigenvalue, Cholesky decompositions (Matrix, spam),
and LU ones, with Chebyshev and Monte Carlo approximations; moments use
the methods due to Martin and Smirnov/Anselin.

## Usage

``` r
do_ldet(coef, env, which=1)
jacobianSetup(method, env, con, pre_eig=NULL, trs=NULL, interval=NULL, which=1)
cheb_setup(env, q=5, which=1)
mcdet_setup(env, p=16, m=30, which=1)
eigen_setup(env, which=1)
eigen_pre_setup(env, pre_eig, which=1)
spam_setup(env, pivot="MMD", which=1)
spam_update_setup(env, in_coef=0.1, pivot="MMD", which=1)
Matrix_setup(env, Imult, super=as.logical(NA), which=1)
Matrix_J_setup(env, super=FALSE, which=1)
LU_setup(env, which=1)
LU_prepermutate_setup(env, coef=0.1, order=FALSE, which=1)
moments_setup(env, trs=NULL, m, p, type="MC", correct=TRUE, trunc=TRUE, eq7=TRUE, which=1)
SE_classic_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_whichMin_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_interp_setup(env, SE_method="LU", p=16, m=30, nrho=200,
 interval=c(-1,0.999), which=1)
can.be.simmed(listw)
```

## Arguments

- coef:

  spatial coefficient value

- env:

  environment containing pre-computed objects, fixed after assignment in
  setup functions

- which:

  default 1; if 2, use second listw object

- method:

  string value, used by `jacobianSetup` to choose method

- con:

  control list passed from model fitting function and parsed in
  `jacobianSetup` to set environment variables for method-specific setup

- pre_eig:

  pre-computed eigenvalues of length n

- q:

  Chebyshev approximation order; default in calling spdep functions is
  5, here it cannot be missing and does not have a default

- p:

  Monte Carlo approximation number of random normal variables; default
  calling spdep functions is 16, here it cannot be missing and does not
  have a default

- m:

  Monte Carlo approximation number of series terms; default in calling
  spdep functions is 30, here it cannot be missing and does not have a
  default; `m` serves the same purpose in the moments method

- pivot:

  default “MMD”, may also be “RCM” for Cholesky decompisition using spam

- in_coef:

  fill-in initiation coefficient value, default 0.1

- Imult:

  see
  [`Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html);
  numeric scalar which defaults to zero. The matrix that is decomposed
  is A+m\*I where m is the value of Imult and I is the identity matrix
  of order ncol(A). Default in calling spdep functions is 2, here it
  cannot be missing and does not have a default, but is rescaled for
  binary weights matrices in proportion to the maximim row sum in those
  calling functions

- super:

  see
  [`Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html);
  logical scalar indicating is a supernodal decomposition should be
  created. The alternative is a simplicial decomposition. Default in
  calling spdep functions is FALSE for “Matrix_J” and `as.logical(NA)`
  for “Matrix”. Setting it to NA leaves the choice to a CHOLMOD-internal
  heuristic

- order:

  default FALSE; used in LU_prepermutate, note warnings given for `lu`
  method

- trs:

  A numeric vector of `m` traces, as from `trW`

- type:

  moments trace type, see
  [`trW`](https://r-spatial.github.io/spatialreg/reference/trW.md)

- correct:

  default TRUE: use Smirnov correction term, see
  [`trW`](https://r-spatial.github.io/spatialreg/reference/trW.md)

- trunc:

  default TRUE: truncate Smirnov correction term, see
  [`trW`](https://r-spatial.github.io/spatialreg/reference/trW.md)

- eq7:

  default TRUE; use equation 7 in Smirnov and Anselin (2009), if FALSE
  no unit root correction

- SE_method:

  default “LU”, alternatively “MC”; underlying lndet method to use for
  generating SE toolbox emulation grid

- nrho:

  default 200, number of lndet values in first stage SE toolbox
  emulation grid

- interval:

  default c(-1,0.999) if interval argument NULL, bounds for SE toolbox
  emulation grid

- interpn:

  default 2000, number of lndet values to interpolate in second stage SE
  toolbox emulation grid

- SElndet:

  default NULL, used to pass a pre-computed two-column matrix of
  coefficient values and corresponding interpolated lndet values

- listw:

  a spatial weights object

## Details

Since environments are containers in the R workspace passed by reference
rather than by value, they are useful for passing objects to functions
called in numerical optimisation, here for the maximum likelihood
estimation of spatial regression models. This technique can save a
little time on each function call, balanced against the need to access
the objects in the environment inside the function. The environment
should contain a `family` string object either “SAR”, “CAR” or “SMA”
(used in `do_ldet` to choose spatial moving average in `spautolm`, and
these specific objects before calling the set-up functions:

- eigen:

  Classical Ord eigenvalue computations - either:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  verbose

  :   logical scalar: legacy report print control, for historical
      reasons only

  or:

  pre_eig

  :   pre-computed eigenvalues

  and assigns to the environment:

  eig

  :   a vector of eigenvalues

  eig.range

  :   the search interval for the spatial coefficient

  method

  :   string: “eigen”

- Matrix:

  Sparse matrix pre-computed Cholesky decomposition with fast updating:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  and assigns to the environment:

  csrw

  :   sparse spatial weights matrix

  nW

  :   negative sparse spatial weights matrix

  pChol

  :   a “CHMfactor” from factorising `csrw` with
      [`Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html)

  nChol

  :   a “CHMfactor” from factorising `nW` with
      [`Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html)

  method

  :   string: “Matrix”

- Matrix_J:

  Standard Cholesky decomposition without updating:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  n

  :   number of spatial objects

  and assigns to the environment:

  csrw

  :   sparse spatial weights matrix

  I

  :   sparse identity matrix

  super

  :   the value of the `super` argument

  method

  :   string: “Matrix_J”

- spam:

  Standard Cholesky decomposition without updating:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  n

  :   number of spatial objects

  and assigns to the environment:

  csrw

  :   sparse spatial weights matrix

  I

  :   sparse identity matrix

  pivot

  :   string — pivot method

  method

  :   string: “spam”

- spam_update:

  Pre-computed Cholesky decomposition with updating:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  n

  :   number of spatial objects

  and assigns to the environment:

  csrw

  :   sparse spatial weights matrix

  I

  :   sparse identity matrix

  csrwchol

  :   A Cholesky decomposition for updating

  method

  :   string: “spam”

- LU:

  Standard LU decomposition without updating:

  listw

  :   A listw spatial weights object

  n

  :   number of spatial objects

  and assigns to the environment:

  W

  :   sparse spatial weights matrix

  I

  :   sparse identity matrix

  method

  :   string: “LU”

- LU_prepermutate:

  Standard LU decomposition with updating (pre-computed fill-reducing
  permutation):

  listw

  :   A listw spatial weights object

  n

  :   number of spatial objects

  and assigns to the environment:

  W

  :   sparse spatial weights matrix

  lu_order

  :   order argument to lu

  pq

  :   2-column matrix for row and column permutation for fill-reduction

  I

  :   sparse identity matrix

  method

  :   string: “LU”

- MC:

  Monte Carlo approximation:

  listw

  :   A listw spatial weights object

  and assigns to the environment:

  clx

  :   list of Monte Carlo approximation terms (the first two simulated
      traces are replaced by their analytical equivalents)

  W

  :   sparse spatial weights matrix

  method

  :   string: “MC”

- cheb:

  Chebyshev approximation:

  listw

  :   A listw spatial weights object

  and assigns to the environment:

  trT

  :   vector of Chebyshev approximation terms

  W

  :   sparse spatial weights matrix

  method

  :   string: “Chebyshev”

- moments:

  moments approximation:

  listw

  :   A listw spatial weights object

  can.sim

  :   logical scalar: can the spatial weights be made symmetric by
      similarity

  and assigns to the environment:

  trs

  :   vector of traces, possibly approximated

  q12

  :   integer vector of length 2, unit roots terms, ignored until 0.5-52

  eq7

  :   logical scalar: use equation 7

  correct

  :   logical scalar: use Smirnov correction term

  trunc

  :   logical scalar: truncate Smirnov correction term

  method

  :   string: “moments”

- SE_classic:

  :

  listw

  :   A listw spatial weights object

  n

  :   number of spatial objects

  and assigns to the environment:

  detval

  :   two column matrix of lndet grid values

  method

  :   string: “SE_classic”

  SE_method

  :   string: “LU” or “MC”

- SE_whichMin:

  :

  listw

  :   A listw spatial weights object

  n

  :   number of spatial objects

  and assigns to the environment:

  detval

  :   two column matrix of lndet grid values

  method

  :   string: “SE_whichMin”

  SE_method

  :   string: “LU” or “MC”

- SE_interp:

  :

  listw

  :   A listw spatial weights object

  n

  :   number of spatial objects

  and assigns to the environment:

  fit

  :   fitted spline object from which to predict lndet values

  method

  :   string: “SE_interp”

  SE_method

  :   string: “LU” or “MC”

Some set-up functions may also assign `similar` to the environment if
the weights were made symmetric by similarity.

Three set-up functions emulate the behaviour of the Spatial Econometrics
toolbox (March 2010) maximum likelihood lndet grid performance. The
toolbox lndet functions compute a smaller number of lndet values for a
grid of coefficient values (spacing 0.01), and then interpolate to a
finer grid of values (spacing 0.001). “SE_classic”, which is an
implementation of the SE toolbox code, for example in f_sar.m, appears
to have selected a row in the grid matrix one below the correct row when
the candidate coefficient value was between 0.005 and 0.01-fuzz, always
rounding the row index down. A possible alternative is to choose the
index that is closest to the candidate coefficient value
(“SE_whichMin”). Another alternative is to fit a spline model to the
first stage coarser grid, and pass this fitted model to the log
likelihood function to make a point prediction using the candidate
coefficient value, rather than finding the grid index (“SE_interp”).

## Value

`do_ldet` returns the value of the Jacobian for the calculation method
recorded in the environment argument, and for the Monte Carlo
approximation, returns a measure of the spread of the approximation as
an “sd” attribute; the remaining functions modify the environment in
place as a side effect and return nothing.

## References

LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC
Press, Boca Raton, pp. 77–110.

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the
Jacobian in Gaussian spatial autoregressive models: An illustrated
comparison of available methods. *Geographical Analysis*, 45(2),
150-179.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`spautolm`](https://r-spatial.github.io/spatialreg/reference/spautolm.md),
[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`errorsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html)

## Examples

``` r
data(boston, package="spData")
#require("spdep", quietly=TRUE)
lw <- spdep::nb2listw(boston.soi)
can.sim <- can.be.simmed(lw)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
eigen_setup(env)
get("similar", envir=env)
#> [1] TRUE
do_ldet(0.5, env)
#> [1] -18.26702
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
eigen_pre_setup(env, pre_eig=eigenw(similar.listw(lw)))
do_ldet(0.5, env)
#> [1] -18.26702
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
Matrix_setup(env, Imult=2, super=FALSE)
get("similar", envir=env)
#> [1] TRUE
do_ldet(0.5, env)
#> [1] -18.26702
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
spam_setup(env)
get("similar", envir=env)
#> [1] TRUE
do_ldet(0.5, env)
#> [1] -18.26702
#> attr(,"logarithm")
#> [1] TRUE
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_setup(env)
get("similar", envir=env)
#> [1] FALSE
do_ldet(0.5, env)
#> [1] -18.26702
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_prepermutate_setup(env)
get("similar", envir=env)
#> [1] FALSE
do_ldet(0.5, env)
#> [1] -18.26702
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
cheb_setup(env, q=5)
get("similar", envir=env)
#> [1] FALSE
do_ldet(0.5, env)
#> [1] -18.26176
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
set.seed(12345)
mcdet_setup(env, p=16, m=30)
get("similar", envir=env)
#> [1] FALSE
do_ldet(0.5, env)
#> [1] -18.38606
#> attr(,"sd")
#> [1] 0.2045107
rm(env)
```
