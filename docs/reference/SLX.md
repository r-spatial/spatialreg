# Spatial Durbin linear (SLX, spatially lagged X) model

`lmSLX` fits an `lm` model augmented with the spatially lagged RHS
variables, including the lagged intercept when the spatial weights are
not row-standardised. `create_WX` creates spatially lagged RHS
variables, and is exposed for use in model fitting functions.

## Usage

``` r
lmSLX(formula, data = list(), listw, na.action, weights=NULL, Durbin=TRUE,
 zero.policy=NULL, return_impacts=TRUE)
# S3 method for class 'SlX'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
# S3 method for class 'SlX'
summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)
# S3 method for class 'summary.SlX'
print(x, digits = max(3L, getOption("digits") - 3L),
 symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)
# S3 method for class 'SlX'
impacts(obj, ...)
# S3 method for class 'WXimpact'
print(x, ...)
# S3 method for class 'WXimpact'
summary(object, ..., adjust_k=(attr(object, "type") == "SDEM"))
# S3 method for class 'SlX'
predict(object, newdata, listw, zero.policy=NULL, ...)
create_WX(x, listw, zero.policy=NULL, prefix="")
```

## Arguments

- formula:

  a symbolic description of the model to be fit. The details of model
  specification are given for [`lm()`](https://rdrr.io/r/stats/lm.html)

- data:

  an optional data frame containing the variables in the model. By
  default the variables are taken from the environment which the
  function is called.

- listw:

  a `listw` object created for example by `nb2listw`

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the spatial weights list will be subsetted to remove NAs
  in the data. It may be necessary to set zero.policy to TRUE because
  this subsetting may create no-neighbour observations. Note that only
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

  default TRUE for `lmSLX` (Durbin model including WX); if TRUE, full
  spatial Durbin model; if a formula object, the subset of explanatory
  variables to lag. From version 1.3-7, the presence of factors
  (categorical variables) in the Durbin term will give a warning, as it
  is as yet unknown how spatial lags of categorical variables should be
  interpreted.

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE assign NA

- return_impacts:

  default TRUE; may be set FALSE to avoid problems calculating impacts
  with aliased variables

- digits:

  the number of significant digits to use when printing

- correlation:

  logical; if `TRUE`, the correlation matrix of the estimated parameters
  is returned and printed

- symbolic.cor:

  logical. If `TRUE`, print the correlations in a symbolic form (see
  'symnum') rather than as numbers

- signif.stars:

  logical. If `TRUE`, 'significance stars' are printed for each
  coefficient

- obj:

  A spatial regression object created by `lmSLX`

- ...:

  Arguments passed through

- prefix:

  default empty string, may be “lag” in some cases

- x, object:

  model matrix to be lagged; lagImpact objects created by `impacts`
  methods

- adjust_k:

  default TRUE if SDEM else FALSE, adjust internal OLS SDEM standard
  errors by dividing by n rather than (n-k) (default changed and bug
  fixed after 0.7-8; standard errors now ML in SDEM summary and impacts
  summary and identical - for SLX use FALSE)

- newdata:

  data frame in which to predict — if NULL, predictions are for the data
  on which the model was fitted. Should have row names corresponding to
  region.id. If row names are exactly the same than the ones used for
  training, it uses in-sample predictors for forecast.

## Value

The `lmSLX` function returns an “lm” object with a “mixedImps” list of
three impact matrixes (impacts and standard errors) for direct, indirect
and total impacts; total impacts and their standard errors calculated
using [`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html).

## Details

From version 1.4.1, functions for models including spatially lagged
independent variables warn on fitting if any of the right-hand side
variables are factors. This is because the interpretation of
coefficients that are not slopes is unclear when the variable is not
interpretable on an unbounded line, such as factors. Factor variable
names are shown with the suffix “(F)”, others “dy/dx” in output from
impact methods. A discussion can be found at
<https://github.com/rsbivand/eqc25_talk>.

## See also

[`lm`](https://rdrr.io/r/stats/lm.html)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Examples

``` r
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb, style="W")
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(COL.SLX)
#> 
#> Call:
#> lm(formula = formula(paste("y ~ ", paste(colnames(x)[-1], collapse = "+"))), 
#>     data = as.data.frame(x), weights = weights)
#> 
#> Coefficients:
#>              Estimate    Std. Error  t value     Pr(>|t|)  
#> (Intercept)   7.503e+01   6.626e+00   1.132e+01   1.261e-14
#> INC          -1.109e+00   3.738e-01  -2.967e+00   4.854e-03
#> HOVAL        -2.897e-01   1.014e-01  -2.858e+00   6.486e-03
#> lag.INC      -1.371e+00   5.613e-01  -2.443e+00   1.867e-02
#> lag.HOVAL     1.918e-01   2.003e-01   9.572e-01   3.437e-01
#> 
summary(impacts(COL.SLX))
#> Impact measures (SlX, glht, n-k):
#>                 Direct   Indirect       Total
#> INC dy/dx   -1.1089293 -1.3709725 -2.47990173
#> HOVAL dy/dx -0.2897283  0.1917608 -0.09796753
#> ========================================================
#> Standard errors:
#>                Direct  Indirect     Total
#> INC dy/dx   0.3738129 0.5612771 0.4965456
#> HOVAL dy/dx 0.1013673 0.2003335 0.2028016
#> ========================================================
#> Z-values:
#>                Direct   Indirect      Total
#> INC dy/dx   -2.966535 -2.4425945 -4.9943086
#> HOVAL dy/dx -2.858202  0.9572079 -0.4830709
#> 
#> p-values:
#>             Direct    Indirect Total     
#> INC dy/dx   0.0030118 0.014582 5.9047e-07
#> HOVAL dy/dx 0.0042605 0.338462 0.62905   
#> 
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=TRUE)
summary(impacts(COL.SLX))
#> Impact measures (SlX, glht, n-k):
#>                        Direct     Indirect       Total
#> INC dy/dx        -0.947594274 -1.275338647 -2.22293292
#> HOVAL dy/dx      -0.777427839 -0.355048446 -1.13247628
#> I(HOVAL^2) dy/dx  0.004639919  0.005608104  0.01024802
#> ========================================================
#> Standard errors:
#>                       Direct   Indirect      Total
#> INC dy/dx        0.398832844 0.59687399 0.56325440
#> HOVAL dy/dx      0.464456540 0.98213200 1.07348981
#> I(HOVAL^2) dy/dx 0.004226259 0.00908385 0.01025212
#> ========================================================
#> Z-values:
#>                     Direct   Indirect      Total
#> INC dy/dx        -2.375918 -2.1366966 -3.9465877
#> HOVAL dy/dx      -1.673844 -0.3615079 -1.0549483
#> I(HOVAL^2) dy/dx  1.097879  0.6173709  0.9996008
#> 
#> p-values:
#>                  Direct   Indirect Total     
#> INC dy/dx        0.017505 0.032623 7.9273e-05
#> HOVAL dy/dx      0.094161 0.717720 0.29145   
#> I(HOVAL^2) dy/dx 0.272258 0.536990 0.31750   
#> 
summary(COL.SLX)
#> 
#> Call:
#> lm(formula = formula(paste("y ~ ", paste(colnames(x)[-1], collapse = "+"))), 
#>     data = as.data.frame(x), weights = weights)
#> 
#> Coefficients:
#>                 Estimate    Std. Error  t value     Pr(>|t|)  
#> (Intercept)      9.246e+01   1.928e+01   4.796e+00   2.058e-05
#> INC             -9.476e-01   3.988e-01  -2.376e+00   2.214e-02
#> HOVAL           -7.774e-01   4.645e-01  -1.674e+00   1.016e-01
#> I.HOVAL.2.       4.640e-03   4.226e-03   1.098e+00   2.785e-01
#> lag.INC         -1.275e+00   5.969e-01  -2.137e+00   3.849e-02
#> lag.HOVAL       -3.550e-01   9.821e-01  -3.615e-01   7.195e-01
#> lag.I.HOVAL.2.   5.608e-03   9.084e-03   6.174e-01   5.403e-01
#> 
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=~INC)
summary(impacts(COL.SLX))
#> Impact measures (SlX, glht, n-k):
#>                        Direct  Indirect        Total
#> INC dy/dx        -1.079064628 -1.010896 -2.089960575
#> HOVAL dy/dx      -0.634518755        NA -0.634518755
#> I(HOVAL^2) dy/dx  0.003455273        NA  0.003455273
#> ========================================================
#> Standard errors:
#>                      Direct  Indirect      Total
#> INC dy/dx        0.38471071 0.4552667 0.44650193
#> HOVAL dy/dx      0.44760078        NA 0.44760078
#> I(HOVAL^2) dy/dx 0.00411036        NA 0.00411036
#> ========================================================
#> Z-values:
#>                      Direct  Indirect      Total
#> INC dy/dx        -2.8048728 -2.220448 -4.6807425
#> HOVAL dy/dx      -1.4175997        NA -1.4175997
#> I(HOVAL^2) dy/dx  0.8406254        NA  0.8406254
#> 
#> p-values:
#>                  Direct    Indirect Total     
#> INC dy/dx        0.0050336 0.026388 2.8584e-06
#> HOVAL dy/dx      0.1563077 NA       0.15631   
#> I(HOVAL^2) dy/dx 0.4005579 NA       0.40056   
#> 
summary(COL.SLX)
#> 
#> Call:
#> lm(formula = formula(paste("y ~ ", paste(colnames(x)[-1], collapse = "+"))), 
#>     data = as.data.frame(x), weights = weights)
#> 
#> Coefficients:
#>              Estimate    Std. Error  t value     Pr(>|t|)  
#> (Intercept)   8.368e+01   9.265e+00   9.032e+00   1.401e-11
#> INC          -1.079e+00   3.847e-01  -2.805e+00   7.466e-03
#> HOVAL        -6.345e-01   4.476e-01  -1.418e+00   1.634e-01
#> I.HOVAL.2.    3.455e-03   4.110e-03   8.406e-01   4.051e-01
#> lag.INC      -1.011e+00   4.553e-01  -2.220e+00   3.159e-02
#> 
COL.SLX <- lmSLX(CRIME ~ INC, data=COL.OLD, listw=lw)
summary(COL.SLX)
#> 
#> Call:
#> lm(formula = formula(paste("y ~ ", paste(colnames(x)[-1], collapse = "+"))), 
#>     data = as.data.frame(x), weights = weights)
#> 
#> Coefficients:
#>              Estimate    Std. Error  t value     Pr(>|t|)  
#> (Intercept)   7.398e+01   6.208e+00   1.192e+01   1.155e-15
#> INC          -1.589e+00   3.564e-01  -4.458e+00   5.276e-05
#> lag.INC      -1.086e+00   4.812e-01  -2.257e+00   2.882e-02
#> 
summary(impacts(COL.SLX))
#> Impact measures (SlX, glht, n-k):
#>              Direct  Indirect     Total
#> INC dy/dx -1.588901 -1.085867 -2.674768
#> ========================================================
#> Standard errors:
#>              Direct  Indirect    Total
#> INC dy/dx 0.3564039 0.4811809 0.407313
#> ========================================================
#> Z-values:
#>              Direct  Indirect     Total
#> INC dy/dx -4.458147 -2.256671 -6.566861
#> 
#> p-values:
#>           Direct     Indirect Total     
#> INC dy/dx 8.2671e-06 0.024029 5.1387e-11
#> 
# \dontrun{
crds <- cbind(COL.OLD$X, COL.OLD$Y)
mdist <- sqrt(sum(diff(apply(crds, 2, range))^2))
dnb <- spdep::dnearneigh(crds, 0, mdist)
dists <- spdep::nbdists(dnb, crds)
f <- function(x, form, data, dnb, dists, verbose) {
  glst <- lapply(dists, function(d) 1/(d^x))
  lw <- spdep::nb2listw(dnb, glist=glst, style="B")
  res <- logLik(lmSLX(form=form, data=data, listw=lw))
  if (verbose) cat("power:", x, "logLik:", res, "\n")
  res
}
opt <- optimize(f, interval=c(0.1, 4), form=CRIME ~ INC + HOVAL,
 data=COL.OLD, dnb=dnb, dists=dists, verbose=TRUE, maximum=TRUE)
#> power: 1.589667 logLik: -172.6864 
#> power: 2.510333 logLik: -177.741 
#> power: 1.020665 logLik: -171.5379 
#> power: 0.8721475 logLik: -171.7979 
#> power: 1.11302 logLik: -171.4973 
#> power: 1.107329 logLik: -171.497 
#> power: 1.105705 logLik: -171.497 
#> power: 1.105746 logLik: -171.497 
#> power: 1.105664 logLik: -171.497 
#> power: 1.105705 logLik: -171.497 
glst <- lapply(dists, function(d) 1/(d^opt$maximum))
lwa <- spdep::nb2listw(dnb, glist=glst, style="B")
SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lwa)
summary(SLX)
#> 
#> Call:
#> lm(formula = formula(paste("y ~ ", paste(colnames(x)[-1], collapse = "+"))), 
#>     data = as.data.frame(x), weights = weights)
#> 
#> Coefficients:
#>                  Estimate  Std. Error  t value   Pr(>|t|)
#> (Intercept)      18.78498  11.23795     1.67157   0.10187
#> INC              -0.65428   0.29463    -2.22064   0.03169
#> HOVAL            -0.18244   0.08138    -2.24174   0.03019
#> lag..Intercept.   6.27790   4.44284     1.41304   0.16484
#> lag.INC          -0.12958   0.28794    -0.45002   0.65496
#> lag.HOVAL         0.02668   0.11273     0.23664   0.81406
#> 
summary(impacts(SLX))
#> Impact measures (SlX, glht, n-k):
#>                 Direct    Indirect      Total
#> INC dy/dx   -0.6542760 -0.12957739 -0.7838534
#> HOVAL dy/dx -0.1824383  0.02667561 -0.1557627
#> ========================================================
#> Standard errors:
#>                 Direct  Indirect     Total
#> INC dy/dx   0.29463342 0.2879383 0.3742251
#> HOVAL dy/dx 0.08138257 0.1127279 0.1404591
#> ========================================================
#> Z-values:
#>                Direct   Indirect     Total
#> INC dy/dx   -2.220644 -0.4500179 -2.094604
#> HOVAL dy/dx -2.241736  0.2366371 -1.108954
#> 
#> p-values:
#>             Direct   Indirect Total   
#> INC dy/dx   0.026375 0.65270  0.036206
#> HOVAL dy/dx 0.024978 0.81294  0.267450
#> 
# }
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
pslx0 <- predict(COL.SLX)
pslx1 <- predict(COL.SLX, newdata=COL.OLD, listw=lw)
all.equal(pslx0, pslx1)
#> [1] TRUE
COL.OLD1 <- COL.OLD
COL.OLD1$INC <- COL.OLD1$INC + 1
pslx2 <- predict(COL.SLX, newdata=COL.OLD1, listw=lw)
sum(coef(COL.SLX)[c(2,4)])
#> [1] -2.479902
mean(pslx2-pslx1)
#> [1] -2.479902
```
