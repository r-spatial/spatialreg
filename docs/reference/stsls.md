# Generalized spatial two stage least squares

The function fits a spatial lag model by two stage least squares, with
the option of adjusting the results for heteroskedasticity.

## Usage

``` r
stsls(formula, data = list(), listw, zero.policy = NULL,
 na.action = na.fail, robust = FALSE, HC=NULL, legacy=FALSE, W2X = TRUE,
 sig2n_k=TRUE, adjust.n=FALSE)
# S3 method for class 'Stsls'
impacts(obj, ..., tr, R = NULL, listw = NULL, evalues=NULL, Q=NULL)
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

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE (default) assign
  NA - causing
  [`lagsarlm()`](https://r-spatial.github.io/spatialreg/reference/ML_models.md)
  to terminate with an error

- na.action:

  a function (default `na.fail`), can also be `na.omit` or `na.exclude`
  with consequences for residuals and fitted values - in these cases the
  weights list will be subsetted to remove NAs in the data. It may be
  necessary to set zero.policy to TRUE because this subsetting may
  create no-neighbour observations. Note that only weights lists created
  without using the glist argument to `nb2listw` may be subsetted.

- robust:

  default FALSE, if TRUE, apply a heteroskedasticity correction to the
  coefficients covariances

- HC:

  default NULL, if `robust` is TRUE, assigned “HC0”, may take values
  “HC0” or “HC1” for White estimates or MacKinnon-White estimates
  respectively

- legacy:

  the argument chooses between two implementations of the robustness
  correction: default FALSE - use the estimate of Omega only in the
  White consistent estimator of the variance-covariance matrix, if TRUE,
  use the original implementation which runs a GLS using the estimate of
  Omega, overrides sig2n_k, and yields different coefficient estimates
  as well - see example below

- W2X:

  default TRUE, if FALSE only WX are used as instruments in the spatial
  two stage least squares; until release 0.4-60, only WX were used - see
  example below; Python `spreg::GM_Lag` is default FALSE

- sig2n_k:

  default TRUE - use n-k to calculate sigma^2, if FALSE use n; Python
  `spreg::GM_Lag` is default FALSE

- adjust.n:

  default FALSE, used in creating spatial weights constants for the
  Anselin-Kelejian (1997) test

- obj:

  A spatial regression object created by `lagsarlm`, `lagmess` or by
  `lmSLX`; in `HPDinterval.LagImpact`, a LagImpact object

- ...:

  Arguments passed through to methods in the coda package

- tr:

  A vector of traces of powers of the spatial weights matrix created
  using `trW`, for approximate impact measures; if not given, `listw`
  must be given for exact measures (for small to moderate spatial
  weights matrices); the traces must be for the same spatial weights as
  were used in fitting the spatial regression, and must be
  row-standardised

- evalues:

  vector of eigenvalues of spatial weights matrix for impacts
  calculations

- R:

  If given, simulations are used to compute distributions for the impact
  measures, returned as `mcmc` objects; the objects are used for
  convenience but are not output by an MCMC process

- Q:

  default NULL, else an integer number of cumulative power series
  impacts to calculate if `tr` is given

## Details

The fitting implementation fits a spatial lag model:

\$\$y = \rho W y + X \beta + \varepsilon\$\$

by using spatially lagged X variables as instruments for the spatially
lagged dependent variable.

From version 1.3-6, the general Anselin-Kelejian (1997) test for
residual spatial autocorrelation is added.

## Value

an object of class "Stsls" containing:

- coefficients:

  coefficient estimates

- var:

  coefficient covariance matrix

- sse:

  sum of squared errors

- residuals:

  model residuals

- df:

  degrees of freedom

## References

Kelejian, H.H. and I.R. Prucha (1998). A generalized spatial two stage
least squares procedure for estimating a spatial autoregressive model
with autoregressive disturbances. *Journal of Real Estate Finance and
Economics* 17, 99-121.
[doi:10.1023/A:1007707430416](https://doi.org/10.1023/A%3A1007707430416)
.

Anselin, L., & Kelejian, H. H. (1997). Testing for Spatial Error
Autocorrelation in the Presence of Endogenous Regressors. International
Regional Science Review, 20(1-2), 153-182.
[doi:10.1177/016001769702000109](https://doi.org/10.1177/016001769702000109)
.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of
Estimation Methods for Spatial Econometrics. *Journal of Statistical
Software*, 63(18), 1-36.
[doi:10.18637/jss.v063.i18](https://doi.org/10.18637/jss.v063.i18) .

## Author

Luc Anselin, Gianfranco Piras and Roger Bivand

## See also

[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md)

## Examples

``` r
data(oldcol, package="spdep")
#require(spdep, quietly=TRUE)
lw <- spdep::nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
summary(COL.lag.eig, correlation=TRUE)
#> 
#> Call:lagsarlm(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.68585  -5.35636   0.05421   6.02013  23.20555 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 45.079250   7.177347  6.2808 3.369e-10
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
#>  Correlation of coefficients 
#>             sigma rho   (Intercept) INC  
#> rho         -0.14                        
#> (Intercept)  0.12 -0.83                  
#> INC         -0.05  0.35 -0.61            
#> HOVAL       -0.01  0.08 -0.25       -0.44
#> 
COL.lag.stsls <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
(x <- summary(COL.lag.stsls, correlation=TRUE))
#> 
#> Call:stsls(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -37.86437  -5.65096  -0.13669   6.23315  22.90823 
#> 
#> Coefficients: 
#>              Estimate Std. Error t value  Pr(>|t|)
#> Rho          0.454567   0.185118  2.4555  0.014067
#> (Intercept) 43.793442  10.952229  3.9986 6.372e-05
#> INC         -1.000716   0.383858 -2.6070  0.009134
#> HOVAL       -0.265489   0.091852 -2.8904  0.003847
#> 
#> Residual variance (sigma squared): 103.44, (sigma: 10.171)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 0.043401, p-value: 0.83497
#> 
#> Correlation of Coefficients:
#>             Rho   (Intercept) INC  
#> (Intercept) -0.92                  
#> INC          0.63 -0.76            
#> HOVAL        0.04 -0.16       -0.36
#> 
coef(x)
#>               Estimate  Std. Error   t value     Pr(>|t|)
#> Rho          0.4545669  0.18511845  2.455546 1.406706e-02
#> (Intercept) 43.7934425 10.95222944  3.998587 6.372175e-05
#> INC         -1.0007158  0.38385778 -2.606996 9.134037e-03
#> HOVAL       -0.2654890  0.09185167 -2.890410 3.847398e-03
W <- as(lw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
loobj1 <- impacts(COL.lag.stsls, R=200, tr=trMatc)
summary(loobj1, zstats=TRUE, short=TRUE)
#> Impact measures (lag, trace):
#>                 Direct   Indirect     Total
#> INC dy/dx   -1.0607411 -0.7739768 -1.834718
#> HOVAL dy/dx -0.2814136 -0.2053354 -0.486749
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct Indirect    Total
#> INC dy/dx   0.4037287 1.635349 1.685099
#> HOVAL dy/dx 0.1150209 1.323526 1.387249
#> 
#> Simulated z-values:
#>                Direct   Indirect      Total
#> INC dy/dx   -2.673468 -0.6208495 -1.2430491
#> HOVAL dy/dx -2.528055 -0.2991168 -0.4949854
#> 
#> Simulated p-values:
#>             Direct    Indirect Total  
#> INC dy/dx   0.0075071 0.53470  0.21385
#> HOVAL dy/dx 0.0114696 0.76485  0.62061
ev <- eigenw(lw)
loobj2 <- impacts(COL.lag.stsls, R=200, evalues=ev)
summary(loobj2, zstats=TRUE, short=TRUE)
#> Impact measures (lag, evalues):
#>                 Direct   Indirect     Total
#> INC dy/dx   -1.0607411 -0.7739768 -1.834718
#> HOVAL dy/dx -0.2814136 -0.2053354 -0.486749
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct  Indirect     Total
#> INC dy/dx   0.3827268 0.5916403 0.7619622
#> HOVAL dy/dx 0.1022968 0.2215918 0.2889025
#> 
#> Simulated z-values:
#>                Direct  Indirect     Total
#> INC dy/dx   -2.875069 -1.339200 -2.483969
#> HOVAL dy/dx -2.713089 -1.069861 -1.781268
#> 
#> Simulated p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0040394 0.18051  0.012993
#> HOVAL dy/dx 0.0066659 0.28468  0.074869
require(coda)
HPDinterval(loobj1)
#>                  lower      upper
#> INC dy/dx   -1.7997982 -0.2313338
#> HOVAL dy/dx -0.5510265 -0.1110313
#> attr(,"Probability")
#> [1] 0.95
COL.lag.stslsW <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw, W2X=FALSE)
summary(COL.lag.stslsW, correlation=TRUE)
#> 
#> Call:stsls(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw, 
#>     W2X = FALSE)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -37.785778  -5.442414  -0.052649   6.170104  23.039123 
#> 
#> Coefficients: 
#>              Estimate Std. Error t value  Pr(>|t|)
#> Rho          0.444202   0.189141  2.3485  0.018848
#> (Intercept) 44.359512  11.157079  3.9759 7.011e-05
#> INC         -1.014319   0.387469 -2.6178  0.008850
#> HOVAL       -0.265681   0.091954 -2.8893  0.003861
#> 
#> Residual variance (sigma squared): 103.67, (sigma: 10.182)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 0.058156, p-value: 0.80943
#> 
#> Correlation of Coefficients:
#>             Rho   (Intercept) INC  
#> (Intercept) -0.93                  
#> INC          0.64 -0.77            
#> HOVAL        0.04 -0.16       -0.36
#> 
COL.lag.stslsWn <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw, W2X=FALSE, sig2n_k=FALSE)
summary(COL.lag.stslsWn, correlation=TRUE)
#> 
#> Call:stsls(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw, 
#>     W2X = FALSE, sig2n_k = FALSE)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -37.785778  -5.442414  -0.052649   6.170104  23.039123 
#> 
#> Coefficients: 
#>              Estimate Std. Error t value  Pr(>|t|)
#> Rho          0.444202   0.181257  2.4507  0.014259
#> (Intercept) 44.359512  10.691995  4.1489 3.341e-05
#> INC         -1.014319   0.371318 -2.7317  0.006301
#> HOVAL       -0.265681   0.088121 -3.0150  0.002570
#> 
#> Residual variance (sigma squared): 95.203, (sigma: 9.7572)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 0.058156, p-value: 0.80943
#> 
#> Correlation of Coefficients:
#>             Rho   (Intercept) INC  
#> (Intercept) -0.93                  
#> INC          0.64 -0.77            
#> HOVAL        0.04 -0.16       -0.36
#> 
COL.lag.stslsR <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
robust=TRUE, W2X=FALSE)
summary(COL.lag.stslsR, correlation=TRUE)
#> 
#> Call:stsls(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw, 
#>     robust = TRUE, W2X = FALSE)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -37.785778  -5.442414  -0.052649   6.170104  23.039123 
#> 
#> Coefficients: 
#>             Estimate HC0 std. Error z value  Pr(>|z|)
#> Rho          0.44420        0.13748  3.2310  0.001234
#> (Intercept) 44.35951        7.67306  5.7812 7.417e-09
#> INC         -1.01432        0.44113 -2.2993  0.021486
#> HOVAL       -0.26568        0.17353 -1.5311  0.125752
#> 
#> Residual variance (sigma squared): 103.67, (sigma: 10.182)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 0.056852, p-value: 0.81154
#> 
#> Correlation of Coefficients:
#>             Rho   (Intercept) INC  
#> (Intercept) -0.90                  
#> INC          0.15 -0.28            
#> HOVAL        0.24 -0.24       -0.83
#> 
COL.lag.stslsRl <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
robust=TRUE, legacy=TRUE, W2X=FALSE)
summary(COL.lag.stslsRl, correlation=TRUE)
#> 
#> Call:stsls(formula = CRIME ~ INC + HOVAL, data = COL.OLD, listw = lw, 
#>     robust = TRUE, legacy = TRUE, W2X = FALSE)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -38.654607  -5.141303  -0.065221   5.864384  23.671589 
#> 
#> Coefficients: 
#>             Estimate HC0 std. Error z value  Pr(>|z|)
#> Rho          0.40138        0.13554  2.9613  0.003064
#> (Intercept) 47.37696        7.49975  6.3171 2.664e-10
#> INC         -1.15183        0.43490 -2.6485  0.008085
#> HOVAL       -0.25047        0.17333 -1.4450  0.148461
#> 
#> Asymptotic robust residual variance: 96.446, (sigma: 9.8207)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 0.10254, p-value: 0.7488
#> 
#> Correlation of Coefficients:
#>             Rho   (Intercept) INC  
#> (Intercept) -0.89                  
#> INC          0.12 -0.26            
#> HOVAL        0.25 -0.26       -0.83
#> 
data(boston, package="spData")
gp2a <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, spdep::nb2listw(boston.soi))
summary(gp2a)
#> 
#> Call:stsls(formula = log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = spdep::nb2listw(boston.soi))
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.5356002 -0.0758562 -0.0045074  0.0719613  0.7128012 
#> 
#> Coefficients: 
#>                Estimate  Std. Error  t value  Pr(>|t|)
#> Rho          4.5925e-01  3.8485e-02  11.9330 < 2.2e-16
#> (Intercept)  2.4025e+00  2.1710e-01  11.0661 < 2.2e-16
#> CRIM        -7.3557e-03  1.0345e-03  -7.1100 1.160e-12
#> ZN           3.6435e-04  3.9311e-04   0.9268 0.3540112
#> INDUS        1.1992e-03  1.8365e-03   0.6530 0.5137794
#> CHAS1        1.1929e-02  2.6632e-02   0.4479 0.6542202
#> I(NOX^2)    -2.8874e-01  9.2546e-02  -3.1199 0.0018091
#> I(RM^2)      6.6991e-03  1.0192e-03   6.5728 4.938e-11
#> AGE         -2.5810e-04  4.0940e-04  -0.6304 0.5284073
#> log(DIS)    -1.6043e-01  2.6107e-02  -6.1451 7.993e-10
#> log(RAD)     7.1704e-02  1.4926e-02   4.8038 1.557e-06
#> TAX         -3.6857e-04  9.5315e-05  -3.8668 0.0001103
#> PTRATIO     -1.2957e-02  4.1334e-03  -3.1347 0.0017203
#> B            2.8845e-04  8.0266e-05   3.5937 0.0003261
#> log(LSTAT)  -2.3984e-01  2.2470e-02 -10.6740 < 2.2e-16
#> 
#> Residual variance (sigma squared): 0.020054, (sigma: 0.14161)
#> Anselin-Kelejian (1997) test for residual autocorrelation
#> test value: 4.3942, p-value: 0.036062
#> 
```
