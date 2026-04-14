# Spatial simultaneous autoregressive model estimation by maximum likelihood

The `lagsarlm` function provides Maximum likelihood estimation of
spatial simultaneous autoregressive lag and spatial Durbin (mixed)
models of the form:

\$\$y = \rho W y + X \beta + \varepsilon\$\$

where \\\rho\\ is found by
[`optimize()`](https://rdrr.io/r/stats/optimize.html) first, and
\\\beta\\ and other parameters by generalized least squares subsequently
(one-dimensional search using optim performs badly on some platforms).
In the spatial Durbin (mixed) model, the spatially lagged independent
variables are added to X. Note that interpretation of the fitted
coefficients should use impact measures, because of the feedback loops
induced by the data generation process for this model. With one of the
sparse matrix methods, larger numbers of observations can be handled,
but the `interval=` argument may need be set when the weights are not
row-standardised.

Maximum likelihood estimation of spatial simultaneous autoregressive
error models of the form:

\$\$y = X \beta + u, u = \lambda W u + \varepsilon\$\$

where \\\lambda\\ is found by
[`optimize()`](https://rdrr.io/r/stats/optimize.html) first, and
\\\beta\\ and other parameters by generalized least squares
subsequently. With one of the sparse matrix methods, larger numbers of
observations can be handled, but the `interval=` argument may need be
set when the weights are not row-standardised. When `etype` is “emixed”,
a so-called spatial Durbin error model is fitted.

Maximum likelihood estimation of spatial simultaneous autoregressive
“SAC/SARAR” models of the form:

\$\$y = \rho W1 y + X \beta + u, u = \lambda W2 u + \varepsilon\$\$

where \\\rho\\ and \\\lambda\\ are found by `nlminb` or
[`optim()`](https://rdrr.io/r/stats/optim.html) first, and \\\beta\\ and
other parameters by generalized least squares subsequently.

## Usage

``` r
lagsarlm(formula, data = list(), listw, na.action, Durbin, type,
 method="eigen", quiet=NULL, zero.policy=NULL, interval=NULL,
 tol.solve=.Machine$double.eps, trs=NULL, control=list())
errorsarlm(formula, data=list(), listw, na.action, weights=NULL,
 Durbin, etype, method="eigen", quiet=NULL, zero.policy=NULL,
 interval = NULL, tol.solve=.Machine$double.eps, trs=NULL, control=list())
sacsarlm(formula, data = list(), listw, listw2 = NULL, na.action, Durbin, type,
 method="eigen", quiet=NULL, zero.policy=NULL, tol.solve=.Machine$double.eps,
 llprof=NULL, interval1=NULL, interval2=NULL, trs1=NULL, trs2=NULL,
 control = list())
# S3 method for class 'Sarlm'
summary(object, correlation = FALSE, Nagelkerke = FALSE,
 Hausman=FALSE, adj.se=FALSE, ...)
# S3 method for class 'Sarlm'
print(x, ...)
# S3 method for class 'summary.Sarlm'
print(x, digits = max(5, .Options$digits - 3),
 signif.stars = FALSE, ...)
# S3 method for class 'Sarlm'
residuals(object, ...)
# S3 method for class 'Sarlm'
deviance(object, ...)
# S3 method for class 'Sarlm'
coef(object, ...)
# S3 method for class 'Sarlm'
vcov(object, ...)
# S3 method for class 'Sarlm'
fitted(object, ...)
# S3 method for class 'Sarlm'
nobs(object, ...)
# S3 method for class 'Sarlm'
set_coef(model, coefs, ...)
```

## Arguments

- formula:

  a symbolic description of the model to be fit. The details of model
  specification are given for [`lm()`](https://rdrr.io/r/stats/lm.html)

- data:

  an optional data frame containing the variables in the model. By
  default the variables are taken from the environment which the
  function is called.

- listw, listw2:

  a `listw` object created for example by `nb2listw`; if `nb2listw` not
  given, set to the same spatial weights as the `listw` argument

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the weights list will be subsetted to remove NAs in the
  data. It may be necessary to set zero.policy to TRUE because this
  subsetting may create no-neighbour observations. Note that only
  weights lists created without using the glist argument to `nb2listw`
  may be subsetted.

- weights:

  an optional vector of weights to be used in the fitting process.
  Non-NULL weights can be used to indicate that different observations
  have different variances (with the values in weights being inversely
  proportional to the variances); or equivalently, when the elements of
  weights are positive integers w_i, that each response y_i is the mean
  of w_i unit-weight observations (including the case that there are w_i
  observations equal to y_i and the data have been summarized) -
  [`lm`](https://rdrr.io/r/stats/lm.html)

- Durbin:

  default FALSE (spatial lag, error or SARAR model); if TRUE, full
  spatial Durbin model (SDM, SDEM or GNM); if a formula object, the
  subset of explanatory variables to lag. From version 1.3-7, the
  presence of factors (categorical variables) in the Durbin term will
  give a warning, as it is as yet unknown how spatial lags of
  categorical variables should be interpreted.

- type:

  (use the ‘Durbin=’ argument - retained for backwards compatibility
  only) default "lag", may be set to "mixed"; when "mixed", the lagged
  intercept is dropped for spatial weights style "W", that is
  row-standardised weights, but otherwise included; “Durbin” may be used
  instead of “mixed”

- etype:

  (use the ‘Durbin=’ argument - retained for backwards compatibility
  only) default "error", may be set to "emixed" to include the spatially
  lagged independent variables added to X; when "emixed", the lagged
  intercept is dropped for spatial weights style "W", that is
  row-standardised weights, but otherwise included

- method:

  "eigen" (default) - the Jacobian is computed as the product of (1 -
  rho\*eigenvalue) using `eigenw`, and "spam" or "Matrix_J" for strictly
  symmetric weights lists of styles "B" and "C", or made symmetric by
  similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S",
  using code from the spam or Matrix packages to calculate the
  determinant; “Matrix” and “spam_update” provide updating Cholesky
  decomposition methods; "LU" provides an alternative sparse matrix
  decomposition approach. In addition, there are "Chebyshev" and Monte
  Carlo "MC" approximate log-determinant methods; the
  Smirnov/Anselin (2009) trace approximation is available as "moments".
  Three methods: "SE_classic", "SE_whichMin", and "SE_interp" are
  provided experimentally, the first to attempt to emulate the behaviour
  of Spatial Econometrics toolbox ML fitting functions. All use grids of
  log determinant values, and the latter two attempt to ameliorate some
  features of "SE_classic".

- quiet:

  default NULL, use !verbose global option value; if FALSE, reports
  function values during optimization.

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE (default) assign
  NA - causing `lagsarlm()` to terminate with an error

- interval:

  default is NULL, search interval for autoregressive parameter

- tol.solve:

  the tolerance for detecting linear dependencies in the columns of
  matrices to be inverted - passed to
  [`solve()`](https://rdrr.io/r/base/solve.html) (default=1.0e-10). This
  may be used if necessary to extract coefficient standard errors (for
  instance lowering to 1e-12), but errors in
  [`solve()`](https://rdrr.io/r/base/solve.html) may constitute
  indications of poorly scaled variables: if the variables have scales
  differing much from the autoregressive coefficient, the values in this
  matrix may be very different in scale, and inverting such a matrix is
  analytically possible by definition, but numerically unstable;
  rescaling the RHS variables alleviates this better than setting
  tol.solve to a very small value

- llprof:

  default NULL, can either be an integer, to divide the feasible ranges
  into a grid of points, or a two-column matrix of spatial coefficient
  values, at which to evaluate the likelihood function

- trs1, trs2:

  default NULL, if given, vectors for each weights object of powered
  spatial weights matrix traces output by `trW`; when given, used in
  some Jacobian methods

- interval1, interval2:

  default is NULL, search intervals for each weights object for
  autoregressive parameters

- trs:

  default NULL, if given, a vector of powered spatial weights matrix
  traces output by `trW`; when given, insert the asymptotic analytical
  values into the numerical Hessian instead of the approximated values;
  may be used to get around some problems raised when the numerical
  Hessian is poorly conditioned, generating NaNs in subsequent
  operations; the use of trs is recommended

- control:

  list of extra control arguments - see section below

- object, model:

  `Sarlm` object from `lagsarlm`, `errorsarlm` or `sacsarlm`

- coefs:

  numerical, replacement coefficients to be inserted into a fitted model

- correlation:

  logical; if 'TRUE', the correlation matrix of the estimated parameters
  including sigma is returned and printed (default=FALSE)

- Nagelkerke:

  if TRUE, the Nagelkerke pseudo R-squared is reported

- Hausman:

  if TRUE, the results of the Hausman test for error models are reported

- adj.se:

  if TRUE, adjust the coefficient standard errors for the number of
  fitted coefficients

- x:

  `Sarlm` object from `lagsarlm`, `errorsarlm` or `sacsarlm` in
  `print.Sarlm`, summary object from `summary.Sarlm` for
  `print.summary.Sarlm`

- digits:

  the number of significant digits to use when printing

- signif.stars:

  logical. If TRUE, "significance stars" are printed for each
  coefficient.

- ...:

  further arguments passed to or from other methods

## Details

The asymptotic standard error of \\\rho\\ is only computed when
method=“eigen”, because the full matrix operations involved would be
costly for large n typically associated with the choice of method="spam"
or "Matrix". The same applies to the coefficient covariance matrix.
Taken as the asymptotic matrix from the literature, it is typically
badly scaled, and with the elements involving \\\rho\\ (lag model) or
\\\lambda\\ (error model) being very small, while other parts of the
matrix can be very large (often many orders of magnitude in difference).
It often happens that the `tol.solve` argument needs to be set to a
smaller value than the default, or the RHS variables can be centred or
reduced in range.

Versions of the package from 0.4-38 include numerical Hessian values
where asymptotic standard errors are not available. This change has been
introduced to permit the simulation of distributions for impact
measures. The warnings made above with regard to variable scaling also
apply in this case.

Note that the fitted() function for the output object assumes that the
response variable may be reconstructed as the sum of the trend, the
signal, and the noise (residuals). Since the values of the response
variable are known, their spatial lags are used to calculate signal
components (Cressie 1993, p. 564). This differs from other software,
including GeoDa, which does not use knowledge of the response variable
in making predictions for the fitting data. Refer to the help page of
[`predict.Sarlm`](https://r-spatial.github.io/spatialreg/reference/predict.sarlm.md)
for discussions and references.

Because numerical optimisation is used to find the values of lambda and
rho in `sacsarlm`, care needs to be shown. It has been found that the
surface of the 2D likelihood function often forms a “banana trench” from
(low rho, high lambda) through (high rho, high lambda) to (high rho, low
lambda) values. In addition, sometimes the banana has optima towards
both ends, one local, the other global, and conseqently the choice of
the starting point for the final optimization becomes crucial. The
default approach is not to use just (0, 0) as a starting point, nor the
(rho, lambda) values from `gstsls`, which lie in a central part of the
“trench”, but either four values at (low rho, high lambda), (0, 0),
(high rho, high lambda), and (high rho, low lambda), and to use the best
of these start points for the final optimization. Optionally, nine
points can be used spanning the whole (lower, upper) space.

From version 1.4.1, functions for models including spatially lagged
independent variables warn on fitting if any of the right-hand side
variables are factors. This is because the interpretation of
coefficients that are not slopes is unclear when the variable is not
interpretable on an unbounded line, such as factors. Factor variable
names are shown with the suffix “(F)”, others “dy/dx” in output from
impact methods. A discussion can be found at
<https://github.com/rsbivand/eqc25_talk>.

## Control arguments

- tol.opt::

  the desired accuracy of the optimization - passed to
  [`optimize()`](https://rdrr.io/r/stats/optimize.html) (default=square
  root of double precision machine tolerance, a larger root may be used
  needed, see help(boston) for an example)

- returnHcov::

  (error model) default TRUE, return the Vo matrix for a spatial Hausman
  test

- pWOrder::

  (error model) default 250, if returnHcov=TRUE and the method is not
  “eigen”, pass this order to `powerWeights` as the power series maximum
  limit

- fdHess::

  default NULL, then set to (method != "eigen") internally; use `fdHess`
  to compute an approximate Hessian using finite differences when using
  sparse matrix methods; used to make a coefficient covariance matrix
  when the number of observations is large; may be turned off to save
  resources if need be. If TRUE, use of compiled sum of squared error
  (SSE) calculation is controlled by `compiled_sse` below

- optimHess::

  default FALSE, use `fdHess` from nlme, if TRUE, use `optim` to
  calculate Hessian at optimum

- optimHessMethod::

  default “optimHess”, may be “nlm” or one of the `optim` methods

- compiled_sse::

  default FALSE for `lagsarlm` but TRUE for `errorsarlm` as sum of
  squared error (SSE) calculation is more demanding; logical value used
  in the log likelihood function to choose compiled code for computing
  SSE

- Imult::

  default 2; used for preparing the Cholesky decompositions for updating
  in the Jacobian function

- super::

  if NULL (default), set to FALSE to use a simplicial decomposition for
  the sparse Cholesky decomposition and method “Matrix_J”, set to
  `as.logical(NA)` for method “Matrix”, if TRUE, use a supernodal
  decomposition

- cheb_q::

  default 5; highest power of the approximating polynomial for the
  Chebyshev approximation

- MC_p::

  default 16; number of random variates

- MC_m::

  default 30; number of products of random variates matrix and spatial
  weights matrix

- spamPivot::

  default “MMD”, alternative “RCM”

- in_coef:

  default 0.1, coefficient value for initial Cholesky decomposition in
  “spam_update”

- type:

  default “MC”, used with method “moments”; alternatives “mult” and
  “moments”, for use if `trs` is missing,
  [`trW`](https://r-spatial.github.io/spatialreg/reference/trW.md)

- correct:

  default TRUE, used with method “moments” to compute the
  Smirnov/Anselin correction term

- trunc:

  default TRUE, used with method “moments” to truncate the
  Smirnov/Anselin correction term

- SE_method:

  default “LU”, may be “MC”

- nrho:

  default 200, as in SE toolbox; the size of the first stage lndet grid;
  it may be reduced to for example 40

- interpn:

  default 2000, as in SE toolbox; the size of the second stage lndet
  grid

- small_asy:

  default TRUE; if the method is not “eigen”, use asymmetric covariances
  rather than numerical Hessian ones if n \<= small

- small:

  default 1500; threshold number of observations for asymmetric
  covariances when the method is not “eigen”

- SElndet:

  default NULL, may be used to pass a pre-computed SE toolbox style
  matrix of coefficients and their lndet values to the "SE_classic" and
  "SE_whichMin" methods

- LU_order:

  default FALSE; used in “LU_prepermutate”, note warnings given for `lu`
  method

- pre_eig:

  default NULL; may be used to pass a pre-computed vector of eigenvalues

- return_impacts:

  default TRUE; may be set FALSE to avoid problems calculating impacts
  with aliased variables

- OrdVsign:

  default 1; used to set the sign of the final component to negative if
  -1 (alpha times ((sigma squared) squared) in Ord (1975) equation B.1).

- opt_method::

  default “nlminb”, may be set to “L-BFGS-B” to use box-constrained
  optimisation in `optim`

- opt_control::

  default [`list()`](https://rdrr.io/r/base/list.html), a control list
  to pass to `nlminb` or `optim`

- pars::

  default `NULL`, for which five trial starting values spanning the
  lower/upper range are tried and the best selected, starting values of
  \\\rho\\ and \\\lambda\\

- npars:

  default integer `4L`, four trial points; if not default value, nine
  trial points

- pre_eig1, pre_eig2:

  default NULL; may be used to pass pre-computed vectors of eigenvalues

## References

Cliff, A. D., Ord, J. K. 1981 *Spatial processes*, Pion; Ord, J. K. 1975
Estimation methods for models of spatial interaction, *Journal of the
American Statistical Association*, 70, 120-126; Anselin, L. 1988
*Spatial econometrics: methods and models.* (Dordrecht: Kluwer);
Anselin, L. 1995 SpaceStat, a software program for the analysis of
spatial data, version 1.80. Regional Research Institute, West Virginia
University, Morgantown, WV; Anselin L, Bera AK (1998) Spatial dependence
in linear regression models with an introduction to spatial
econometrics. In: Ullah A, Giles DEA (eds) Handbook of applied economic
statistics. Marcel Dekker, New York, pp. 237-289; Nagelkerke NJD (1991)
A note on a general definition of the coefficient of determination.
Biometrika 78: 691-692; Cressie, N. A. C. 1993 *Statistics for spatial
data*, Wiley, New York; LeSage J and RK Pace (2009) Introduction to
Spatial Econometrics. CRC Press, Boca Raton.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of
Estimation Methods for Spatial Econometrics. *Journal of Statistical
Software*, 63(18), 1-36.
[doi:10.18637/jss.v063.i18](https://doi.org/10.18637/jss.v063.i18) .

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the
Jacobian in Gaussian spatial autoregressive models: An illustrated
comparison of available methods. *Geographical Analysis*, 45(2),
150-179.

## Author

Roger Bivand <Roger.Bivand@nhh.no>, with thanks to Andrew Bernat for
contributions to the asymptotic standard error code.

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`impacts`](https://r-spatial.github.io/spatialreg/reference/impacts.md)

## Examples

``` r
data(oldcol, package="spdep")
listw <- spdep::nb2listw(COL.nb, style="W")
ev <- eigenw(listw)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", quiet=FALSE, control=list(pre_eig=ev, OrdVsign=1))
#> 
#> Spatial lag model
#> Jacobian calculated using neighbourhood matrix eigenvalues
#> 
#> rho:  -0.5674437     function value:  -202.2909 
#> rho:  0.03126655     function value:  -186.749 
#> rho:  0.4012898  function value:  -182.419 
#> rho:  0.6138418  function value:  -183.5636 
#> rho:  0.4157662  function value:  -182.398 
#> rho:  0.4295565  function value:  -182.3905 
#> rho:  0.4311288  function value:  -182.3904 
#> rho:  0.4310273  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
(x <- summary(COL.lag.eig, correlation=TRUE))
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "eigen", quiet = FALSE, control = list(pre_eig = ev, 
#>         OrdVsign = 1))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079251   7.177347  6.2808 3.369e-10
#> INC         -1.031616   0.305143 -3.3808 0.0007229
#> HOVAL       -0.265926   0.088499 -3.0049 0.0026570
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Asymptotic standard error: 0.11768
#>     z-value: 3.6626, p-value: 0.00024962
#> Wald statistic: 13.415, p-value: 0.00024962
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: 0.31955, p-value: 0.57188
#> 
#>  Correlation of coefficients 
#>             sigma rho   (Intercept) INC  
#> rho         -0.14                        
#> (Intercept)  0.12 -0.83                  
#> INC         -0.05  0.35 -0.61            
#> HOVAL       -0.01  0.08 -0.25       -0.44
#> 
coef(x)
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept) 45.0792505 7.17734654  6.280768 3.369041e-10
#> INC         -1.0316157 0.30514297 -3.380762 7.228517e-04
#> HOVAL       -0.2659263 0.08849862 -3.004863 2.657002e-03
# \dontrun{
COL.lag.eig$fdHess
#> [1] FALSE
COL.lag.eig$resvar
#>                    sigma           rho (Intercept)         INC         HOVAL
#> sigma       379.77510363 -0.3236306318  16.3015079 -0.29590801 -0.0202478459
#> rho          -0.32363063  0.0138487533  -0.6975717  0.01266245  0.0008664428
#> (Intercept)  16.30150787 -0.6975716740  51.5143034 -1.32602704 -0.1616379856
#> INC          -0.29590801  0.0126624511  -1.3260270  0.09311223 -0.0117959715
#> HOVAL        -0.02024785  0.0008664428  -0.1616380 -0.01179597  0.0078320057
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", control=list(pre_eig=ev, OrdVsign=-1))
summary(COL.lag.eigb)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "eigen", control = list(pre_eig = ev, OrdVsign = -1))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079251   9.617835  4.6870 2.772e-06
#> INC         -1.031616   0.326524 -3.1594  0.001581
#> HOVAL       -0.265926   0.088855 -2.9928  0.002764
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Asymptotic standard error: 0.17322
#>     z-value: 2.4884, p-value: 0.012833
#> Wald statistic: 6.1919, p-value: 0.012833
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: -0.93825, p-value: 1
#> 
COL.lag.eigb$fdHess
#> [1] FALSE
COL.lag.eigb$resvar
#>                    sigma          rho (Intercept)         INC        HOVAL
#> sigma       388.59742932 -0.701154266  35.3176451 -0.64109248 -0.043867490
#> rho          -0.70115427  0.030003688  -1.5113074  0.02743353  0.001877171
#> (Intercept)  35.31764509 -1.511307358  92.5027556 -2.07005707 -0.212549096
#> INC          -0.64109248  0.027433533  -2.0700571  0.10661800 -0.010871824
#> HOVAL        -0.04386749  0.001877171  -0.2125491 -0.01087182  0.007895242
# force numerical Hessian
COL.lag.eig1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25))
summary(COL.lag.eig1)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "Matrix", control = list(small = 25))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (numerical Hessian approximate standard errors) 
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 45.079249   7.870857  5.7274 1.02e-08
#> INC         -1.031616   0.328409 -3.1413 0.001682
#> HOVAL       -0.265926   0.088216 -3.0145 0.002574
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Approximate (numerical Hessian) standard error: 0.12361
#>     z-value: 3.487, p-value: 0.0004885
#> Wald statistic: 12.159, p-value: 0.0004885
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> 
COL.lag.eig1$fdHess
#>                       rho (Intercept)         INC         HOVAL
#> rho          0.0152792403  -0.8344615  0.02005365  0.0002838387
#> (Intercept) -0.8344615335  61.9503925 -1.78346075 -0.1334797491
#> INC          0.0200536508  -1.7834607  0.10785231 -0.0122138975
#> HOVAL        0.0002838387  -0.1334797 -0.01221390  0.0077819895
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25), trs=trMatc)
summary(COL.lag.eig1a)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "Matrix", trs = trMatc, control = list(small = 25))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (numerical Hessian approximate standard errors) 
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 45.079249   7.937000  5.6796 1.35e-08
#> INC         -1.031616   0.329329 -3.1325 0.001733
#> HOVAL       -0.265926   0.088222 -3.0143 0.002576
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Approximate (numerical Hessian) standard error: 0.12502
#>     z-value: 3.4476, p-value: 0.00056553
#> Wald statistic: 11.886, p-value: 0.00056553
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> 
COL.lag.eig1a$fdHess
#>                   sigma2           rho (Intercept)         INC         HOVAL
#> sigma2      380.74788212 -0.3652578483  19.9480308 -0.47938159 -0.0067837895
#> rho          -0.36525785  0.0156300573  -0.8536131  0.02051362  0.0002902909
#> (Intercept)  19.94803083 -0.8536130500  62.9959621 -1.80853220 -0.1338484559
#> INC          -0.47938159  0.0205136229  -1.8085322  0.10845751 -0.0122072022
#> HOVAL        -0.00678379  0.0002902909  -0.1338485 -0.01220720  0.0077831884
COL.lag.eig$resvar[2,2]
#> [1] 0.01384875
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb$resvar[2,2]
#> [1] 0.03000369
# force numerical Hessian
COL.lag.eig1$fdHess[1,1]
#> [1] 0.01527924
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a$fdHess[2,2]
#> [1] 0.01563006
# }
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix", quiet=FALSE))
#> 
#> Spatial lag model
#> Jacobian calculated using sparse matrix Cholesky decomposition
#> rho:  -0.2364499     function value:  -192.9523 
#> rho:  0.2354499  function value:  -183.542 
#> rho:  0.5271001  function value:  -182.7039 
#> rho:  0.4455543  function value:  -182.3974 
#> rho:  0.4267907  function value:  -182.391 
#> rho:  0.4311986  function value:  -182.3904 
#> rho:  0.4310114  function value:  -182.3904 
#> rho:  0.4310231  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310902  function value:  -182.3904 
#> rho:  0.4310488  function value:  -182.3904 
#> rho:  0.431033   function value:  -182.3904 
#> rho:  0.431027   function value:  -182.3904 
#> rho:  0.4310247  function value:  -182.3904 
#> rho:  0.4310238  function value:  -182.3904 
#> rho:  0.4310234  function value:  -182.3904 
#> rho:  0.4310233  function value:  -182.3904 
#> rho:  0.4310233  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> Computing eigenvalues ...
#> 
#>    user  system elapsed 
#>   0.157   0.000   0.158 
summary(COL.lag.M)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "Matrix", quiet = FALSE)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079249   7.177346  6.2808 3.369e-10
#> INC         -1.031616   0.305143 -3.3808 0.0007229
#> HOVAL       -0.265926   0.088499 -3.0049 0.0026570
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Asymptotic standard error: 0.11768
#>     z-value: 3.6626, p-value: 0.00024962
#> Wald statistic: 13.415, p-value: 0.00024962
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: 0.31954, p-value: 0.57188
#> 
impacts(COL.lag.M, listw=listw)
#> Impact measures (lag, exact):
#>                 Direct   Indirect      Total
#> INC dy/dx   -1.0860220 -0.7270849 -1.8131068
#> HOVAL dy/dx -0.2799509 -0.1874254 -0.4673763
# \dontrun{
system.time(COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="spam", quiet=FALSE))
#> 
#> Spatial lag model
#> Jacobian calculated using sparse matrix Cholesky decomposition
#> rho:  -0.2364499     function value:  -192.9523 
#> rho:  0.2354499  function value:  -183.542 
#> rho:  0.5271001  function value:  -182.7039 
#> rho:  0.4455543  function value:  -182.3974 
#> rho:  0.4267907  function value:  -182.391 
#> rho:  0.4311986  function value:  -182.3904 
#> rho:  0.4310114  function value:  -182.3904 
#> rho:  0.4310231  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310902  function value:  -182.3904 
#> rho:  0.4310233  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> rho:  0.4310232  function value:  -182.3904 
#> Computing eigenvalues ...
#> 
#>    user  system elapsed 
#>   0.436   0.000   0.438 
summary(COL.lag.sp)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "spam", quiet = FALSE)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079249   7.177346  6.2808 3.369e-10
#> INC         -1.031616   0.305143 -3.3808 0.0007229
#> HOVAL       -0.265926   0.088499 -3.0049 0.0026570
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Asymptotic standard error: 0.11768
#>     z-value: 3.6626, p-value: 0.00024962
#> Wald statistic: 13.415, p-value: 0.00024962
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: 0.31954, p-value: 0.57188
#> 
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"), control=list(pre_eig=ev))
summary(COL.lag.B)
#> 
#> Call:
#> lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = spdep::nb2listw(COL.nb, 
#>     style = "B"), control = list(pre_eig = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -33.6620  -4.8615  -1.3576   5.1567  25.7563 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 51.604815   6.075285  8.4942 < 2.2e-16
#> INC         -1.154463   0.301808 -3.8252 0.0001307
#> HOVAL       -0.251633   0.087612 -2.8721 0.0040773
#> 
#> Rho: 0.054543, LR test value: 13.453, p-value: 0.00024461
#> Asymptotic standard error: 0.014836
#>     z-value: 3.6763, p-value: 0.00023662
#> Wald statistic: 13.515, p-value: 0.00023662
#> 
#> Log likelihood: -180.6507 for lag model
#> ML residual variance (sigma squared): 93.22, (sigma: 9.655)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 371.3, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: 0.006827, p-value: 0.93415
#> 
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9,
 control=list(pre_eig=ev))
summary(COL.mixed.B)
#> 
#> Call:
#> lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = spdep::nb2listw(COL.nb, 
#>     style = "B"), type = "mixed", tol.solve = 1e-09, control = list(pre_eig = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -34.8460  -4.2057  -0.1195   4.6525  21.6112 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                    Estimate  Std. Error z value  Pr(>|z|)
#> (Intercept)      5.4215e+01  5.6639e+00  9.5719 < 2.2e-16
#> INC             -8.2386e-01  3.1643e-01 -2.6036 0.0092262
#> HOVAL           -3.0085e-01  8.5629e-02 -3.5134 0.0004424
#> lag.(Intercept) -7.5493e+00  1.7307e+00 -4.3620 1.289e-05
#> lag.INC          2.1531e-05  1.2216e-01  0.0002 0.9998594
#> lag.HOVAL        7.2458e-02  3.9007e-02  1.8576 0.0632281
#> 
#> Rho: 0.15212, LR test value: 7.435, p-value: 0.0063967
#> Asymptotic standard error: 0.015565
#>     z-value: 9.7735, p-value: < 2.22e-16
#> Wald statistic: 95.522, p-value: < 2.22e-16
#> 
#> Log likelihood: -177.7722 for mixed model
#> ML residual variance (sigma squared): 82.502, (sigma: 9.0831)
#> Number of observations: 49 
#> Number of parameters estimated: 8 
#> AIC: 371.54, (AIC for lm: 376.98)
#> LM test for residual autocorrelation
#> test value: -26.79, p-value: 1
#> 
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, type="mixed", control=list(pre_eig=ev))
summary(COL.mixed.W)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     type = "mixed", control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822413  12.667204  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520283   0.565129 -0.9206 0.3572355
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
COL.mixed.D00 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=TRUE, control=list(pre_eig=ev))
summary(COL.mixed.D00)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = TRUE, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822413  12.667204  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520283   0.565129 -0.9206 0.3572355
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
COL.mixed.D01 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=FALSE, control=list(pre_eig=ev))
summary(COL.mixed.D01)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = FALSE, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079251   7.177347  6.2808 3.369e-10
#> INC         -1.031616   0.305143 -3.3808 0.0007229
#> HOVAL       -0.265926   0.088499 -3.0049 0.0026570
#> 
#> Rho: 0.43102, LR test value: 9.9736, p-value: 0.001588
#> Asymptotic standard error: 0.11768
#>     z-value: 3.6626, p-value: 0.00024962
#> Wald statistic: 13.415, p-value: 0.00024962
#> 
#> Log likelihood: -182.3904 for lag model
#> ML residual variance (sigma squared): 95.494, (sigma: 9.7721)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.78, (AIC for lm: 382.75)
#> LM test for residual autocorrelation
#> test value: 0.31955, p-value: 0.57188
#> 
COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC + HOVAL, control=list(pre_eig=ev))
summary(COL.mixed.D1)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = ~INC + HOVAL, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822413  12.667204  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520283   0.565129 -0.9206 0.3572355
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
f <- CRIME ~ INC + HOVAL
COL.mixed.D2 <- lagsarlm(f, data=COL.OLD, listw,
 Durbin=as.formula(delete.response(terms(f))),
 control=list(pre_eig=ev))
summary(COL.mixed.D2)
#> 
#> Call:
#> lagsarlm(formula = f, data = COL.OLD, listw = listw, Durbin = as.formula(delete.response(terms(f))), 
#>     control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822413  12.667204  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520283   0.565129 -0.9206 0.3572355
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
COL.mixed.D1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC, control=list(pre_eig=ev))
summary(COL.mixed.D1a)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = ~INC, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.82800  -5.85207   0.12047   6.00137  23.19963 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 48.814687  12.198232  4.0018 6.287e-05
#> INC         -1.006620   0.330641 -3.0445  0.002331
#> HOVAL       -0.265514   0.088768 -2.9911  0.002780
#> lag.INC     -0.186684   0.530550 -0.3519  0.724936
#> 
#> Rho: 0.39229, LR test value: 4.5007, p-value: 0.033881
#> Asymptotic standard error: 0.1561
#>     z-value: 2.513, p-value: 0.011971
#> Wald statistic: 6.3151, p-value: 0.011971
#> 
#> Log likelihood: -182.3328 for mixed model
#> ML residual variance (sigma squared): 96.122, (sigma: 9.8042)
#> Number of observations: 49 
#> Number of parameters estimated: 6 
#> AIC: 376.67, (AIC for lm: 379.17)
#> LM test for residual autocorrelation
#> test value: 2.6134, p-value: 0.10596
#> 
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ inc + HOVAL, control=list(pre_eig=ev)))
#> Error in eval(predvars, data, env) : object 'inc' not found
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ DISCBD + HOVAL, control=list(pre_eig=ev)))
#> Error in create_Durbin(Durbin = Durbin, have_factor_preds = have_factor_preds,  : 
#>   WX variables not in X: DISCBD
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, na.action=na.exclude)
COL.lag.NA$na.action
#> 1020 1021 1022 1023 1024 1025 
#>   20   21   22   23   24   25 
#> attr(,"class")
#> [1] "exclude"
COL.lag.NA
#> 
#> Call:
#> lagsarlm(formula = CRIME ~ INC + HOVAL, data = NA.COL.OLD, listw = listw, 
#>     na.action = na.exclude)
#> Type: lag 
#> 
#> Coefficients:
#>         rho (Intercept)         INC       HOVAL 
#>   0.4537820  43.1054982  -0.9267352  -0.2715541 
#> 
#> Log likelihood: -160.8867 
resid(COL.lag.NA)
#>        1001        1002        1003        1004        1005        1006 
#>  -4.4352945 -13.2750950  -2.9233554 -37.3067541   1.3568011  -3.8828028 
#>        1007        1008        1009        1010        1011        1012 
#>   5.9980436 -12.8113317  -4.0482557  18.1813322   5.9714200   1.0035598 
#>        1013        1014        1015        1016        1017        1018 
#>  -1.7490668  -1.7490651   5.8456326  10.1074159  -4.3706893   3.3713771 
#>        1019        1020        1021        1022        1023        1024 
#>  -4.0823023          NA          NA          NA          NA          NA 
#>        1025        1026        1027        1028        1029        1030 
#>          NA  -6.0553296 -11.5813312  -8.0011389  -1.9047479   3.9744245 
#>        1031        1032        1033        1034        1035        1036 
#>   4.8148614   6.0673484  10.2059848  22.7419914  -2.0830996   0.1422778 
#>        1037        1038        1039        1040        1041        1042 
#>   8.6436388   8.8150865   6.4974684  15.4121790   9.4182148   5.6427605 
#>        1043        1044        1045        1046        1047        1048 
#>  -7.5727122  -7.5983734  -9.7792621 -10.0870765  -3.1242394   3.4921795 
#>        1049 
#>   0.7173251 
COL.lag.NA1 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, Durbin=~INC) # https://github.com/r-spatial/spatialreg/issues/10
COL.lag.NA1$na.action
#> 1020 1021 1022 1023 1024 1025 
#>   20   21   22   23   24   25 
#> attr(,"class")
#> [1] "omit"
COL.lag.NA2 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, Durbin=~INC, na.action=na.exclude)
COL.lag.NA2$na.action
#> 1020 1021 1022 1023 1024 1025 
#>   20   21   22   23   24   25 
#> attr(,"class")
#> [1] "exclude"
# https://github.com/r-spatial/spatialreg/issues/11
COL.lag.NA3 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, control=list(pre_eig=ev))
#> Warning: NAs found, precomputed eigenvalues ignored
COL.lag.NA3$na.action
#> 1020 1021 1022 1023 1024 1025 
#>   20   21   22   23   24   25 
#> attr(,"class")
#> [1] "omit"
# }

# \dontrun{
data(boston, package="spData")
gp2mM <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, spdep::nb2listw(boston.soi), type="mixed", method="Matrix")
#> Warning: use of spatially lagged factor (categorical variable)
#> CHAS
#> is not well-understood
summary(gp2mM)
#> 
#> Call:lagsarlm(formula = log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = spdep::nb2listw(boston.soi), 
#>     type = "mixed", method = "Matrix")
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.6316833 -0.0629790 -0.0090776  0.0682421  0.6991072 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                   Estimate  Std. Error  z value  Pr(>|z|)
#> (Intercept)     1.89816225  0.22759605   8.3400 < 2.2e-16
#> CRIM           -0.00571021  0.00093504  -6.1069 1.016e-09
#> ZN              0.00069091  0.00051869   1.3320  0.182851
#> INDUS          -0.00111343  0.00307354  -0.3623  0.717155
#> CHAS1          -0.04163225  0.02739364  -1.5198  0.128567
#> I(NOX^2)       -0.01034950  0.19360214  -0.0535  0.957367
#> I(RM^2)         0.00794979  0.00102063   7.7891 6.661e-15
#> AGE            -0.00128789  0.00048920  -2.6326  0.008473
#> log(DIS)       -0.12404108  0.09510940  -1.3042  0.192168
#> log(RAD)        0.05863502  0.02257078   2.5978  0.009382
#> TAX            -0.00049084  0.00012145  -4.0416 5.308e-05
#> PTRATIO        -0.01319853  0.00595352  -2.2169  0.026628
#> B               0.00056383  0.00011089   5.0847 3.682e-07
#> log(LSTAT)     -0.24724454  0.02262033 -10.9302 < 2.2e-16
#> lag.CRIM       -0.00464215  0.00172935  -2.6843  0.007267
#> lag.ZN         -0.00037937  0.00070584  -0.5375  0.590940
#> lag.INDUS       0.00025064  0.00385901   0.0649  0.948215
#> lag.CHAS1       0.12518252  0.04071559   3.0746  0.002108
#> lag.I(NOX^2)   -0.38640403  0.22157523  -1.7439  0.081177
#> lag.I(RM^2)    -0.00451252  0.00153180  -2.9459  0.003220
#> lag.AGE         0.00149678  0.00068418   2.1877  0.028693
#> lag.log(DIS)   -0.00453785  0.10046478  -0.0452  0.963973
#> lag.log(RAD)   -0.00940702  0.03104930  -0.3030  0.761912
#> lag.TAX         0.00041083  0.00017867   2.2994  0.021481
#> lag.PTRATIO     0.00060355  0.00788837   0.0765  0.939012
#> lag.B          -0.00050781  0.00014155  -3.5874  0.000334
#> lag.log(LSTAT)  0.09846780  0.03399423   2.8966  0.003772
#> 
#> Rho: 0.59578, LR test value: 181.68, p-value: < 2.22e-16
#> Asymptotic standard error: 0.038445
#>     z-value: 15.497, p-value: < 2.22e-16
#> Wald statistic: 240.16, p-value: < 2.22e-16
#> 
#> Log likelihood: 300.6131 for mixed model
#> ML residual variance (sigma squared): 0.016011, (sigma: 0.12654)
#> Number of observations: 506 
#> Number of parameters estimated: 29 
#> AIC: -543.23, (AIC for lm: -363.55)
#> LM test for residual autocorrelation
#> test value: 29.772, p-value: 4.8604e-08
#> 
W <- as(spdep::nb2listw(boston.soi), "CsparseMatrix")
trMatb <- trW(W, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, spdep::nb2listw(boston.soi), type="mixed", method="Matrix", 
trs=trMatb)
#> Warning: use of spatially lagged factor (categorical variable)
#> CHAS
#> is not well-understood
summary(gp2mMi)
#> 
#> Call:lagsarlm(formula = log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = spdep::nb2listw(boston.soi), 
#>     type = "mixed", method = "Matrix", trs = trMatb)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.6316833 -0.0629790 -0.0090776  0.0682421  0.6991072 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                   Estimate  Std. Error  z value  Pr(>|z|)
#> (Intercept)     1.89816225  0.22759605   8.3400 < 2.2e-16
#> CRIM           -0.00571021  0.00093504  -6.1069 1.016e-09
#> ZN              0.00069091  0.00051869   1.3320  0.182851
#> INDUS          -0.00111343  0.00307354  -0.3623  0.717155
#> CHAS1          -0.04163225  0.02739364  -1.5198  0.128567
#> I(NOX^2)       -0.01034950  0.19360214  -0.0535  0.957367
#> I(RM^2)         0.00794979  0.00102063   7.7891 6.661e-15
#> AGE            -0.00128789  0.00048920  -2.6326  0.008473
#> log(DIS)       -0.12404108  0.09510940  -1.3042  0.192168
#> log(RAD)        0.05863502  0.02257078   2.5978  0.009382
#> TAX            -0.00049084  0.00012145  -4.0416 5.308e-05
#> PTRATIO        -0.01319853  0.00595352  -2.2169  0.026628
#> B               0.00056383  0.00011089   5.0847 3.682e-07
#> log(LSTAT)     -0.24724454  0.02262033 -10.9302 < 2.2e-16
#> lag.CRIM       -0.00464215  0.00172935  -2.6843  0.007267
#> lag.ZN         -0.00037937  0.00070584  -0.5375  0.590940
#> lag.INDUS       0.00025064  0.00385901   0.0649  0.948215
#> lag.CHAS1       0.12518252  0.04071559   3.0746  0.002108
#> lag.I(NOX^2)   -0.38640403  0.22157523  -1.7439  0.081177
#> lag.I(RM^2)    -0.00451252  0.00153180  -2.9459  0.003220
#> lag.AGE         0.00149678  0.00068418   2.1877  0.028693
#> lag.log(DIS)   -0.00453785  0.10046478  -0.0452  0.963973
#> lag.log(RAD)   -0.00940702  0.03104930  -0.3030  0.761912
#> lag.TAX         0.00041083  0.00017867   2.2994  0.021481
#> lag.PTRATIO     0.00060355  0.00788837   0.0765  0.939012
#> lag.B          -0.00050781  0.00014155  -3.5874  0.000334
#> lag.log(LSTAT)  0.09846780  0.03399423   2.8966  0.003772
#> 
#> Rho: 0.59578, LR test value: 181.68, p-value: < 2.22e-16
#> Asymptotic standard error: 0.038445
#>     z-value: 15.497, p-value: < 2.22e-16
#> Wald statistic: 240.16, p-value: < 2.22e-16
#> 
#> Log likelihood: 300.6131 for mixed model
#> ML residual variance (sigma squared): 0.016011, (sigma: 0.12654)
#> Number of observations: 506 
#> Number of parameters estimated: 29 
#> AIC: -543.23, (AIC for lm: -363.55)
#> LM test for residual autocorrelation
#> test value: 29.772, p-value: 4.8604e-08
#> 
# }
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, quiet=FALSE, control=list(pre_eig=ev))
#> 
#> Spatial autoregressive error model
#> 
#> Jacobian calculated using neighbourhood matrix eigenvalues
#> 
#> lambda: -0.5674437  function: -195.8051  Jacobian: -1.636549  SSE: 7936.201 
#> lambda: 0.03126655  function: -187.0219  Jacobian: -0.005318373  SSE: 5927.009 
#> lambda: 0.4012898  function: -183.8422  Jacobian: -0.9953987  SSE: 4999.419 
#> lambda: 0.6299767  function: -183.4895  Jacobian: -2.818134  SSE: 4574.641 
#> lambda: 0.5811116  function: -183.3887  Jacobian: -2.314073  SSE: 4650.566 
#> lambda: 0.554104  function: -183.3817  Jacobian: -2.066354  SSE: 4696.49 
#> lambda: 0.5621834  function: -183.3805  Jacobian: -2.138326  SSE: 4682.474 
#> lambda: 0.5617028  function: -183.3805  Jacobian: -2.133995  SSE: 4683.301 
#> lambda: 0.5617888  function: -183.3805  Jacobian: -2.134769  SSE: 4683.153 
#> lambda: 0.5617902  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
#> lambda: 0.5617902  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
#> lambda: 0.5617902  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
summary(COL.errW.eig)
#> 
#> Call:errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     quiet = FALSE, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -34.81174  -6.44031  -0.72142   7.61476  23.33626 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 59.893220   5.366162 11.1613 < 2.2e-16
#> INC         -0.941312   0.330569 -2.8476 0.0044057
#> HOVAL       -0.302250   0.090476 -3.3407 0.0008358
#> 
#> Lambda: 0.56179, LR test value: 7.9935, p-value: 0.0046945
#> Asymptotic standard error: 0.13387
#>     z-value: 4.1966, p-value: 2.7098e-05
#> Wald statistic: 17.611, p-value: 2.7098e-05
#> 
#> Log likelihood: -183.3805 for error model
#> ML residual variance (sigma squared): 95.575, (sigma: 9.7762)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 376.76, (AIC for lm: 382.75)
#> 
COL.errW.eig_ev <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, control=list(pre_eig=ev))
all.equal(coefficients(COL.errW.eig), coefficients(COL.errW.eig_ev))
#> [1] TRUE
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"))
summary(COL.errB.eig)
#> 
#> Call:
#> errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = spdep::nb2listw(COL.nb, 
#>     style = "B"))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -32.19010  -5.22646  -0.69952   7.92588  24.23511 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 55.383118   5.449775 10.1625 < 2.2e-16
#> INC         -0.936595   0.319355 -2.9328 0.0033596
#> HOVAL       -0.299857   0.088678 -3.3814 0.0007212
#> 
#> Lambda: 0.12686, LR test value: 10.654, p-value: 0.0010983
#> Asymptotic standard error: 0.021745
#>     z-value: 5.8342, p-value: 5.4044e-09
#> Wald statistic: 34.038, p-value: 5.4044e-09
#> 
#> Log likelihood: -182.0502 for error model
#> ML residual variance (sigma squared): 88.744, (sigma: 9.4204)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 374.1, (AIC for lm: 382.75)
#> 
COL.errW.M <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix", quiet=FALSE, trs=trMatc)
#> 
#> Spatial autoregressive error model
#> 
#> Jacobian calculated using sparse matrix Cholesky decomposition
#> lambda: -0.2364499  function: -190.4287  Jacobian: -0.2899343  SSE: 6732.556 
#> lambda: 0.2354499  function: -185.0024  Jacobian: -0.3201148  SSE: 5388.346 
#> lambda: 0.5271001  function: -183.4053  Jacobian: -1.838306  SSE: 4744.975 
#> lambda: 0.728364  function: -184.1244  Jacobian: -4.106761  SSE: 4454.203 
#> lambda: 0.5304163  function: -183.4009  Jacobian: -1.865308  SSE: 4738.889 
#> lambda: 0.5557478  function: -183.3812  Jacobian: -2.080853  SSE: 4693.62 
#> lambda: 0.6216813  function: -183.4637  Jacobian: -2.727084  SSE: 4586.84 
#> lambda: 0.5627129  function: -183.3805  Jacobian: -2.143105  SSE: 4681.563 
#> lambda: 0.5618852  function: -183.3805  Jacobian: -2.135638  SSE: 4682.987 
#> lambda: 0.5617866  function: -183.3805  Jacobian: -2.134749  SSE: 4683.157 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134783  SSE: 4683.15 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134783  SSE: 4683.151 
#> lambda: 0.5617903  function: -183.3805  Jacobian: -2.134782  SSE: 4683.151 
summary(COL.errW.M)
#> 
#> Call:errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     method = "Matrix", quiet = FALSE, trs = trMatc)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -34.81174  -6.44031  -0.72142   7.61476  23.33626 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 59.893219   5.366163 11.1613 < 2.2e-16
#> INC         -0.941312   0.330569 -2.8476 0.0044057
#> HOVAL       -0.302250   0.090476 -3.3407 0.0008358
#> 
#> Lambda: 0.56179, LR test value: 7.9935, p-value: 0.0046945
#> Asymptotic standard error: 0.13387
#>     z-value: 4.1966, p-value: 2.7098e-05
#> Wald statistic: 17.611, p-value: 2.7098e-05
#> 
#> Log likelihood: -183.3805 for error model
#> ML residual variance (sigma squared): 95.575, (sigma: 9.7762)
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 376.76, (AIC for lm: 382.75)
#> 
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, etype="emixed", control=list(pre_eig=ev))
summary(COL.SDEM.eig)
#> 
#> Call:errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     etype = "emixed", control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.31635  -6.54376  -0.22212   6.44591  23.15801 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 73.545133   8.783543  8.3731 < 2.2e-16
#> INC         -1.051673   0.319514 -3.2915 0.0009966
#> HOVAL       -0.275608   0.091151 -3.0236 0.0024976
#> lag.INC     -1.156711   0.578629 -1.9991 0.0456024
#> lag.HOVAL    0.111691   0.198993  0.5613 0.5746048
#> 
#> Lambda: 0.4254, LR test value: 4.9871, p-value: 0.025537
#> Asymptotic standard error: 0.15842
#>     z-value: 2.6852, p-value: 0.0072485
#> Wald statistic: 7.2103, p-value: 0.0072485
#> 
#> Log likelihood: -181.5846 for error model
#> ML residual variance (sigma squared): 92.531, (sigma: 9.6193)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 377.17, (AIC for lm: 380.16)
#> 
# \dontrun{
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=TRUE, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
#> 
#> Call:errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = TRUE, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.31635  -6.54376  -0.22212   6.44591  23.15801 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 73.545133   8.783543  8.3731 < 2.2e-16
#> INC         -1.051673   0.319514 -3.2915 0.0009966
#> HOVAL       -0.275608   0.091151 -3.0236 0.0024976
#> lag.INC     -1.156711   0.578629 -1.9991 0.0456024
#> lag.HOVAL    0.111691   0.198993  0.5613 0.5746048
#> 
#> Lambda: 0.4254, LR test value: 4.9871, p-value: 0.025537
#> Asymptotic standard error: 0.15842
#>     z-value: 2.6852, p-value: 0.0072485
#> Wald statistic: 7.2103, p-value: 0.0072485
#> 
#> Log likelihood: -181.5846 for error model
#> ML residual variance (sigma squared): 92.531, (sigma: 9.6193)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 377.17, (AIC for lm: 380.16)
#> 
COL.SDEM.eig <- errorsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 listw, Durbin=~INC, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
#> 
#> Call:errorsarlm(formula = CRIME ~ DISCBD + INC + HOVAL, data = COL.OLD, 
#>     listw = listw, Durbin = ~INC, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -34.61867  -7.31993   0.82879   5.92877  17.82211 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 68.961912   6.985784  9.8717 < 2.2e-16
#> DISCBD      -5.412936   2.009281 -2.6940  0.007061
#> INC         -0.899425   0.315174 -2.8537  0.004321
#> HOVAL       -0.202846   0.090125 -2.2507  0.024403
#> lag.INC      0.161858   0.603859  0.2680  0.788669
#> 
#> Lambda: 0.24524, LR test value: 1.2719, p-value: 0.25942
#> Asymptotic standard error: 0.18393
#>     z-value: 1.3334, p-value: 0.18241
#> Wald statistic: 1.7778, p-value: 0.18241
#> 
#> Log likelihood: -178.9895 for error model
#> ML residual variance (sigma squared): 85.935, (sigma: 9.2701)
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 371.98, (AIC for lm: 371.25)
#> 
summary(impacts(COL.SDEM.eig))
#> Impact measures (SDEM, glht, n):
#>                  Direct  Indirect      Total
#> DISCBD dy/dx -5.4129365        NA -5.4129365
#> INC dy/dx    -0.8994251 0.1618581 -0.7375670
#> HOVAL dy/dx  -0.2028457        NA -0.2028457
#> ========================================================
#> Standard errors:
#>                  Direct  Indirect      Total
#> DISCBD dy/dx 2.00928140        NA 2.00928140
#> INC dy/dx    0.31517363 0.6038588 0.67559401
#> HOVAL dy/dx  0.09012493        NA 0.09012493
#> ========================================================
#> Z-values:
#>                 Direct  Indirect     Total
#> DISCBD dy/dx -2.693966        NA -2.693966
#> INC dy/dx    -2.853745 0.2680396 -1.091731
#> HOVAL dy/dx  -2.250717        NA -2.250717
#> 
#> p-values:
#>              Direct    Indirect Total    
#> DISCBD dy/dx 0.0070607 NA       0.0070607
#> INC dy/dx    0.0043207 0.78867  0.2749513
#> HOVAL dy/dx  0.0244035 NA       0.0244035
#> 
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, na.action=na.exclude)
COL.err.NA$na.action
#> 1020 1021 1022 1023 1024 1025 
#>   20   21   22   23   24   25 
#> attr(,"class")
#> [1] "exclude"
COL.err.NA
#> 
#> Call:
#> errorsarlm(formula = CRIME ~ INC + HOVAL, data = NA.COL.OLD, 
#>     listw = listw, na.action = na.exclude)
#> Type: error 
#> 
#> Coefficients:
#>      lambda (Intercept)         INC       HOVAL 
#>   0.5748430  58.2460531  -0.8473028  -0.3024909 
#> 
#> Log likelihood: -161.8763 
resid(COL.err.NA)
#>         1001         1002         1003         1004         1005         1006 
#>  -4.18270827 -11.44133870   0.31874930 -34.47163064   2.42244749  -4.32095076 
#>         1007         1008         1009         1010         1011         1012 
#>   8.66744161 -13.38669921  -1.92276572  17.85753938  -1.11484616  -2.30434797 
#>         1013         1014         1015         1016         1017         1018 
#>  -8.16935131  -5.80500232   0.14973709   5.93191455  -7.03028239   2.39112835 
#>         1019         1020         1021         1022         1023         1024 
#>  -8.95099917           NA           NA           NA           NA           NA 
#>         1025         1026         1027         1028         1029         1030 
#>           NA  -2.52940722  -9.60025406  -6.95635595  -0.43630590   5.98493691 
#>         1031         1032         1033         1034         1035         1036 
#>   6.25882676   7.75527047  10.83413252  23.23927270  -0.05594927   1.43808585 
#>         1037         1038         1039         1040         1041         1042 
#>   9.51995272  12.18295375   8.31031720  17.06834507   7.04418414   7.50088949 
#>         1043         1044         1045         1046         1047         1048 
#>  -7.78485971  -6.79207472  -7.94977561 -11.25362133  -5.68994646   5.04837373 
#>         1049 
#>   2.22497374 
print(system.time(ev <- eigenw(similar.listw(listw))))
#>    user  system elapsed 
#>   0.001   0.000   0.001 
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="eigen", control=list(pre_eig=ev))))
#>    user  system elapsed 
#>   0.161   0.000   0.163 
ocoef <- coefficients(COL.errW.eig)
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="eigen", control=list(pre_eig=ev, LAPACK=FALSE))))
#>    user  system elapsed 
#>   0.161   0.000   0.161 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="eigen", control=list(pre_eig=ev, compiled_sse=TRUE))))
#>    user  system elapsed 
#>   0.159   0.000   0.160 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix_J", control=list(super=TRUE))))
#> Warning: the default value of argument 'sqrt' of method 'determinant(<CHMfactor>, <logical>)' may change from TRUE to FALSE as soon as the next release of Matrix; set 'sqrt' when programming
#>    user  system elapsed 
#>   0.181   0.000   0.182 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix_J", control=list(super=FALSE))))
#>    user  system elapsed 
#>    0.18    0.00    0.18 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix_J", control=list(super=as.logical(NA)))))
#>    user  system elapsed 
#>   0.178   0.000   0.179 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix", control=list(super=TRUE))))
#>    user  system elapsed 
#>   0.163   0.000   0.165 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix", control=list(super=FALSE))))
#>    user  system elapsed 
#>   0.164   0.000   0.165 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="Matrix", control=list(super=as.logical(NA)))))
#>    user  system elapsed 
#>   0.164   0.000   0.165 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="spam", control=list(spamPivot="MMD"))))
#>    user  system elapsed 
#>   0.171   0.000   0.173 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="spam", control=list(spamPivot="RCM"))))
#>    user  system elapsed 
#>   0.173   0.000   0.174 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="spam_update", control=list(spamPivot="MMD"))))
#>    user  system elapsed 
#>   0.169   0.000   0.170 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, method="spam_update", control=list(spamPivot="RCM"))))
#>    user  system elapsed 
#>   0.169   0.000   0.170 
print(all.equal(ocoef, coefficients(COL.errW.eig)))
#> [1] TRUE
# }
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.sacW.eig)
#> 
#> Call:sacsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     control = list(pre_eig1 = ev, pre_eig2 = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.32081  -5.33662  -0.20219   6.59672  23.25604 
#> 
#> Type: sac 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 47.783764   9.902659  4.8253 1.398e-06
#> INC         -1.025894   0.326326 -3.1438  0.001668
#> HOVAL       -0.281651   0.090033 -3.1283  0.001758
#> 
#> Rho: 0.36807
#> Asymptotic standard error: 0.19668
#>     z-value: 1.8714, p-value: 0.061285
#> Lambda: 0.16668
#> Asymptotic standard error: 0.29661
#>     z-value: 0.56196, p-value: 0.57415
#> 
#> LR test value: 10.285, p-value: 0.0058432
#> 
#> Log likelihood: -182.2348 for sac model
#> ML residual variance (sigma squared): 95.604, (sigma: 9.7777)
#> Number of observations: 49 
#> Number of parameters estimated: 6 
#> AIC: 376.47, (AIC for lm: 382.75)
#> 
set.seed(1)
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
#> Impact measures (sac, trace):
#>                 Direct   Indirect      Total
#> INC dy/dx   -1.0632722 -0.5601502 -1.6234224
#> HOVAL dy/dx -0.2919129 -0.1537848 -0.4456977
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                 Direct  Indirect     Total
#> INC dy/dx   0.32656989 0.7312432 0.8605110
#> HOVAL dy/dx 0.09484958 0.2229292 0.2761647
#> 
#> Simulated z-values:
#>                Direct   Indirect     Total
#> INC dy/dx   -3.320010 -0.9264309 -2.047227
#> HOVAL dy/dx -3.126395 -0.8927070 -1.794392
#> 
#> Simulated p-values:
#>             Direct     Indirect Total   
#> INC dy/dx   0.00090014 0.35422  0.040636
#> HOVAL dy/dx 0.00176964 0.37201  0.072751
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 type="sacmixed", control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW.eig)
#> 
#> Call:sacsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     type = "sacmixed", control = list(pre_eig1 = ev, pre_eig2 = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -37.8045  -6.5244  -0.2207   5.9944  22.8691 
#> 
#> Type: sacmixed 
#> Coefficients: (asymptotic standard errors) 
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 50.92026   68.25722  0.7460 0.455664
#> INC         -0.95072    0.44033 -2.1591 0.030841
#> HOVAL       -0.28650    0.09994 -2.8667 0.004148
#> lag.INC     -0.69261    1.69113 -0.4096 0.682132
#> lag.HOVAL    0.20852    0.28702  0.7265 0.467546
#> 
#> Rho: 0.31557
#> Asymptotic standard error: 0.94581
#>     z-value: 0.33365, p-value: 0.73864
#> Lambda: 0.15415
#> Asymptotic standard error: 1.0643
#>     z-value: 0.14484, p-value: 0.88484
#> 
#> LR test value: 12.07, p-value: 0.016837
#> 
#> Log likelihood: -181.3422 for sacmixed model
#> ML residual variance (sigma squared): 93.149, (sigma: 9.6514)
#> Number of observations: 49 
#> Number of parameters estimated: 8 
#> AIC: 378.68, (AIC for lm: 382.75)
#> 
set.seed(1)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
#> Impact measures (sacmixed, trace):
#>                 Direct   Indirect      Total
#> INC dy/dx   -1.0317003 -1.3693141 -2.4010144
#> HOVAL dy/dx -0.2768608  0.1629265 -0.1139344
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct Indirect     Total
#> INC dy/dx   0.3758397 2.103274 2.2342655
#> HOVAL dy/dx 0.1146669 0.871566 0.9375493
#> 
#> Simulated z-values:
#>                Direct   Indirect       Total
#> INC dy/dx   -2.803315 -0.8032408 -1.22771113
#> HOVAL dy/dx -2.317372  0.2677616 -0.03450901
#> 
#> Simulated p-values:
#>             Direct   Indirect Total  
#> INC dy/dx   0.005058 0.42184  0.21956
#> HOVAL dy/dx 0.020483 0.78888  0.97247
COL.msacW1.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=TRUE, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW1.eig)
#> 
#> Call:sacsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     Durbin = TRUE, control = list(pre_eig1 = ev, pre_eig2 = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -37.8045  -6.5244  -0.2207   5.9944  22.8691 
#> 
#> Type: sacmixed 
#> Coefficients: (asymptotic standard errors) 
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 50.92026   68.25722  0.7460 0.455664
#> INC         -0.95072    0.44033 -2.1591 0.030841
#> HOVAL       -0.28650    0.09994 -2.8667 0.004148
#> lag.INC     -0.69261    1.69113 -0.4096 0.682132
#> lag.HOVAL    0.20852    0.28702  0.7265 0.467546
#> 
#> Rho: 0.31557
#> Asymptotic standard error: 0.94581
#>     z-value: 0.33365, p-value: 0.73864
#> Lambda: 0.15415
#> Asymptotic standard error: 1.0643
#>     z-value: 0.14484, p-value: 0.88484
#> 
#> LR test value: 12.07, p-value: 0.016837
#> 
#> Log likelihood: -181.3422 for sacmixed model
#> ML residual variance (sigma squared): 93.149, (sigma: 9.6514)
#> Number of observations: 49 
#> Number of parameters estimated: 8 
#> AIC: 378.68, (AIC for lm: 382.75)
#> 
set.seed(1)
summary(impacts(COL.msacW1.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
#> Impact measures (sacmixed, trace):
#>                 Direct   Indirect      Total
#> INC dy/dx   -1.0317003 -1.3693141 -2.4010144
#> HOVAL dy/dx -0.2768608  0.1629265 -0.1139344
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct Indirect     Total
#> INC dy/dx   0.3758397 2.103274 2.2342655
#> HOVAL dy/dx 0.1146669 0.871566 0.9375493
#> 
#> Simulated z-values:
#>                Direct   Indirect       Total
#> INC dy/dx   -2.803315 -0.8032408 -1.22771113
#> HOVAL dy/dx -2.317372  0.2677616 -0.03450901
#> 
#> Simulated p-values:
#>             Direct   Indirect Total  
#> INC dy/dx   0.005058 0.42184  0.21956
#> HOVAL dy/dx 0.020483 0.78888  0.97247
COL.msacW2.eig <- sacsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD, 
 listw, Durbin= ~ INC, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW2.eig)
#> 
#> Call:sacsarlm(formula = CRIME ~ DISCBD + INC + HOVAL, data = COL.OLD, 
#>     listw = listw, Durbin = ~INC, control = list(pre_eig1 = ev, 
#>         pre_eig2 = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -34.2794  -7.0786   1.0543   6.0019  17.8891 
#> 
#> Type: sacmixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 74.064502  37.738940  1.9625 0.049699
#> DISCBD      -5.707678   3.248629 -1.7569 0.078926
#> INC         -0.900975   0.314996 -2.8603 0.004233
#> HOVAL       -0.203399   0.092982 -2.1875 0.028705
#> lag.INC      0.061568   0.926947  0.0664 0.947043
#> 
#> Rho: -0.077599
#> Asymptotic standard error: 0.57955
#>     z-value: -0.1339, p-value: 0.89348
#> Lambda: 0.29646
#> Asymptotic standard error: 0.51021
#>     z-value: 0.58104, p-value: 0.56121
#> 
#> LR test value: 1.4982, p-value: 0.68269
#> 
#> Log likelihood: -178.963 for sacmixed model
#> ML residual variance (sigma squared): 85.135, (sigma: 9.2269)
#> Number of observations: 49 
#> Number of parameters estimated: 8 
#> AIC: 373.93, (AIC for lm: 369.42)
#> 
summary(impacts(COL.msacW2.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
#> Impact measures (sacmixed, trace):
#>                  Direct   Indirect      Total
#> DISCBD dy/dx -5.7150799 0.41841676 -5.2966632
#> INC dy/dx    -0.9031729 0.12421197 -0.7789609
#> HOVAL dy/dx  -0.2036631 0.01491074 -0.1887524
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                 Direct  Indirect    Total
#> DISCBD dy/dx 3.2143995 5.1439184 5.240818
#> INC dy/dx    0.3744544 1.9189292 2.100533
#> HOVAL dy/dx  0.1102170 0.6188162 0.683367
#> 
#> Simulated z-values:
#>                 Direct    Indirect      Total
#> DISCBD dy/dx -1.831506  0.05485563 -1.0694932
#> INC dy/dx    -2.385739  0.16719446 -0.2725575
#> HOVAL dy/dx  -1.984790 -0.19745864 -0.4989242
#> 
#> Simulated p-values:
#>              Direct   Indirect Total  
#> DISCBD dy/dx 0.067025 0.95625  0.28485
#> INC dy/dx    0.017045 0.86722  0.78519
#> HOVAL dy/dx  0.047168 0.84347  0.61783
# \dontrun{
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, type="mixed", method="eigen")
summary(COL.mix.eig, correlation=TRUE, Nagelkerke=TRUE)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     type = "mixed", method = "eigen")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822415  12.667205  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520284   0.565129 -0.9206 0.3572355
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Nagelkerke pseudo-R-squared: 0.6494 
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
#>  Correlation of coefficients 
#>             sigma rho   (Intercept) INC   HOVAL lag.INC
#> rho         -0.18                                      
#> (Intercept)  0.16 -0.89                                
#> INC         -0.03  0.14 -0.19                          
#> HOVAL        0.02 -0.09  0.03       -0.45              
#> lag.INC     -0.09  0.49 -0.53       -0.36  0.05        
#> lag.HOVAL   -0.04  0.19 -0.36        0.19 -0.24 -0.41  
#> 
COL.mix.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, type="mixed", method="Matrix")
summary(COL.mix.M, correlation=TRUE, Nagelkerke=TRUE)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw, 
#>     type = "mixed", method = "Matrix")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.47829  -6.46731  -0.33835   6.05200  22.62969 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 42.822418  12.667205  3.3806 0.0007233
#> INC         -0.914223   0.331094 -2.7612 0.0057586
#> HOVAL       -0.293738   0.089212 -3.2926 0.0009927
#> lag.INC     -0.520284   0.565129 -0.9206 0.3572354
#> lag.HOVAL    0.245640   0.178917  1.3729 0.1697756
#> 
#> Rho: 0.42634, LR test value: 5.3693, p-value: 0.020494
#> Asymptotic standard error: 0.15623
#>     z-value: 2.7288, p-value: 0.0063561
#> Wald statistic: 7.4465, p-value: 0.0063561
#> 
#> Log likelihood: -181.3935 for mixed model
#> ML residual variance (sigma squared): 91.791, (sigma: 9.5808)
#> Nagelkerke pseudo-R-squared: 0.6494 
#> Number of observations: 49 
#> Number of parameters estimated: 7 
#> AIC: 376.79, (AIC for lm: 380.16)
#> LM test for residual autocorrelation
#> test value: 0.28919, p-value: 0.59074
#> 
#>  Correlation of coefficients 
#>             sigma rho   (Intercept) INC   HOVAL lag.INC
#> rho         -0.18                                      
#> (Intercept)  0.16 -0.89                                
#> INC         -0.03  0.14 -0.19                          
#> HOVAL        0.02 -0.09  0.03       -0.45              
#> lag.INC     -0.09  0.49 -0.53       -0.36  0.05        
#> lag.HOVAL   -0.04  0.19 -0.36        0.19 -0.24 -0.41  
#> 
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
  spdep::nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, correlation=TRUE, Nagelkerke=TRUE, Hausman=TRUE)
#> 
#> Call:
#> errorsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = spdep::nb2listw(COL.nb, 
#>     style = "W"), method = "eigen")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -34.81174  -6.44031  -0.72142   7.61476  23.33626 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 59.893219   5.366163 11.1613 < 2.2e-16
#> INC         -0.941312   0.330569 -2.8476 0.0044057
#> HOVAL       -0.302250   0.090476 -3.3407 0.0008358
#> 
#> Lambda: 0.56179, LR test value: 7.9935, p-value: 0.0046945
#> Asymptotic standard error: 0.13387
#>     z-value: 4.1966, p-value: 2.7098e-05
#> Wald statistic: 17.611, p-value: 2.7098e-05
#> 
#> Log likelihood: -183.3805 for error model
#> ML residual variance (sigma squared): 95.575, (sigma: 9.7762)
#> Nagelkerke pseudo-R-squared: 0.61978 
#> Number of observations: 49 
#> Number of parameters estimated: 5 
#> AIC: 376.76, (AIC for lm: 382.75)
#> Hausman test: 4.902, df: 3, p-value: 0.17911
#> 
#>  Correlation of coefficients 
#>             sigma lambda (Intercept) INC  
#> lambda      -0.24                         
#> (Intercept)  0.00  0.00                   
#> INC          0.00  0.00  -0.56            
#> HOVAL        0.00  0.00  -0.26       -0.45
#> 
# }
```
