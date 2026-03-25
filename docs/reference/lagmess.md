# Matrix exponential spatial lag model

The function fits a matrix exponential spatial lag model, using `optim`
to find the value of `alpha`, the spatial coefficient.

## Usage

``` r
lagmess(formula, data = list(), listw, zero.policy = NULL, na.action = na.fail,
 q = 10, start = -2.5, control=list(), method="BFGS", verbose=NULL,
 use_expm=FALSE)
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

  a `listw` object created for example by
  [`spdep::nb2listw()`](https://r-spatial.github.io/spdep/reference/nb2listw.html)

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE assign NA - causing
  `lagmess()` to terminate with an error

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the weights list will be subsetted to remove NAs in the
  data. It may be necessary to set zero.policy to TRUE because this
  subsetting may create no-neighbour observations. Note that only
  weights lists created without using the glist argument to `nb2listw`
  may be subsetted.

- q:

  default 10; number of powers of the spatial weights to use

- start:

  starting value for numerical optimization, should be a small negative
  number

- control:

  control parameters passed to `optim`

- method:

  default `BFGS`, method passed to `optim`

- verbose:

  default NULL, use global option value; if TRUE report function values
  during optimization

- use_expm:

  default FALSE; if TRUE use
  [`expm::expAtv`](https://rdrr.io/pkg/expm/man/expAtv.html) instead of
  a truncated power series of W

## Details

The underlying spatial lag model:

\$\$y = \rho W y + X \beta + \varepsilon\$\$

where \\\rho\\ is the spatial parameter may be fitted by maximum
likelihood. In that case, the log likelihood function includes the
logarithm of cumbersome Jacobian term \\\|I - \rho W\|\\. If we rewrite
the model as:

\$\$S y = X \beta + \varepsilon\$\$

we see that in the ML case \\S y = (I - \rho W) y\\. If W is
row-stochastic, S may be expressed as a linear combination of
row-stochastic matrices. By pre-computing the matrix \\\[y, Wy, W^2y,
..., W^{q-1}y\]\\, the term \\S y (\alpha)\\ can readily be found by
numerical optimization using the matrix exponential approach. \\\alpha\\
and \\\rho\\ are related as \\\rho = 1 - \exp{\alpha}\\, conditional on
the number of matrix power terms taken `q`.

## Value

The function returns an object of class `Lagmess` with components:

- lmobj:

  the `lm` object returned after fitting `alpha`

- alpha:

  the spatial coefficient

- alphase:

  the standard error of the spatial coefficient using the numerical
  Hessian

- rho:

  the value of `rho` implied by `alpha`

- bestmess:

  the object returned by `optim`

- q:

  the number of powers of the spatial weights used

- start:

  the starting value for numerical optimization used

- na.action:

  (possibly) named vector of excluded or omitted observations if
  non-default na.action argument used

- nullLL:

  the log likelihood of the aspatial model for the same data

## References

J. P. LeSage and R. K. Pace (2007) A matrix exponential specification.
Journal of Econometrics, 140, 190-214; J. P. LeSage and R. K. Pace
(2009) Introduction to Spatial Econometrics. CRC Press, Chapter 9.

## Author

Roger Bivand <Roger.Bivand@nhh.no> and Eric Blankmeyer

## See also

[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`optim`](https://rdrr.io/r/stats/optim.html)

## Examples

``` r
#require(spdep, quietly=TRUE)
data(baltimore, package="spData")
baltimore$AGE <- ifelse(baltimore$AGE < 1, 1, baltimore$AGE)
lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(cbind(baltimore$X, baltimore$Y), k=7)))
obj1 <- lm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT),
 data=baltimore)
spdep::lm.morantest(obj1, lw)
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> weights: lw
#> 
#> Moran I statistic standard deviate = 7.4648, p-value = 4.171e-14
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.245149959     -0.007853660      0.001148722 
#> 
spdep::lm.LMtests(obj1, lw, test="all")
#> Please update scripts to use lm.RStests in place of lm.LMtests
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> test weights: listw
#> 
#> RSerr = 48.648, df = 1, p-value = 3.063e-12
#> 
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> test weights: listw
#> 
#> RSlag = 83.091, df = 1, p-value < 2.2e-16
#> 
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> test weights: listw
#> 
#> adjRSerr = 1.2535, df = 1, p-value = 0.2629
#> 
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> test weights: listw
#> 
#> adjRSlag = 35.696, df = 1, p-value = 2.306e-09
#> 
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data =
#> baltimore)
#> test weights: listw
#> 
#> SARMA = 84.344, df = 2, p-value < 2.2e-16
#> 
system.time(obj2 <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw))
#>    user  system elapsed 
#>    0.03    0.00    0.03 
(x <- summary(obj2))
#> Matrix exponential spatial lag model:
#> 
#> Call:
#> lagmess(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), 
#>     data = baltimore, listw = lw)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -2.026722 -0.141191  0.050831  0.223830  1.073114 
#> 
#> Coefficients:
#>              Estimate Std. Error t value  Pr(>|t|)
#> (Intercept)  1.546376   0.214513  7.2088 1.039e-11
#> PATIO        0.258287   0.086891  2.9726  0.003303
#> log(AGE)    -0.148174   0.035252 -4.2033 3.912e-05
#> log(SQFT)    0.300966   0.071598  4.2036 3.908e-05
#> 
#> Residual standard error: 0.41658 on 207 degrees of freedom
#> Multiple R-squared:  0.22373,    Adjusted R-squared:  0.21248 
#> F-statistic: 19.887 on 3 and 207 DF,  p-value: 2.2881e-11
#> 
#> Alpha: -0.64302, standard error: 0.1043
#>     z-value: -6.1649, p-value: 7.0511e-10
#> LR test value: 48.296, p-value: 3.6644e-12
#> Implied rho: 0.4742995 
#> 
coef(x)
#>               Estimate Std. Error   t value     Pr(>|t|)
#> (Intercept)  1.5463761 0.21451319  7.208769 1.038762e-11
#> PATIO        0.2582874 0.08689079  2.972552 3.303312e-03
#> log(AGE)    -0.1481738 0.03525174 -4.203305 3.912250e-05
#> log(SQFT)    0.3009659 0.07159771  4.203569 3.908038e-05
has_expm <- require("expm", quietly=TRUE)
#> 
#> Attaching package: ‘expm’
#> The following object is masked from ‘package:Matrix’:
#> 
#>     expm
if (has_expm) {
system.time(
obj2a <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw, use_expm=TRUE)
)
summary(obj2a)
}
#> Matrix exponential spatial lag model:
#> (calculated with expm)
#> 
#> Call:
#> lagmess(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), 
#>     data = baltimore, listw = lw, use_expm = TRUE)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -2.026722 -0.141191  0.050831  0.223830  1.073114 
#> 
#> Coefficients:
#>              Estimate Std. Error t value  Pr(>|t|)
#> (Intercept)  1.546376   0.214513  7.2088 1.039e-11
#> PATIO        0.258287   0.086891  2.9726  0.003303
#> log(AGE)    -0.148174   0.035252 -4.2033 3.912e-05
#> log(SQFT)    0.300966   0.071598  4.2036 3.908e-05
#> 
#> Residual standard error: 0.41658 on 207 degrees of freedom
#> Multiple R-squared:  0.22373,    Adjusted R-squared:  0.21248 
#> F-statistic: 19.887 on 3 and 207 DF,  p-value: 2.2881e-11
#> 
#> Alpha: -0.64302, standard error: 0.1043
#>     z-value: -6.1649, p-value: 7.0511e-10
#> LR test value: 48.296, p-value: 3.6644e-12
#> Implied rho: 0.4742995 
#> 
obj3 <- lagsarlm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw)
summary(obj3)
#> 
#> Call:lagsarlm(formula = log(PRICE) ~ PATIO + log(AGE) + log(SQFT), 
#>     data = baltimore, listw = lw)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.969554 -0.147514  0.042581  0.199567  1.076181 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)  1.255885   0.320679  3.9163 8.991e-05
#> PATIO        0.244225   0.083448  2.9267 0.0034262
#> log(AGE)    -0.131947   0.034510 -3.8235 0.0001316
#> log(SQFT)    0.278888   0.070259  3.9694 7.205e-05
#> 
#> Rho: 0.55765, LR test value: 52.661, p-value: 3.9635e-13
#> Asymptotic standard error: 0.072749
#>     z-value: 7.6653, p-value: 1.7764e-14
#> Wald statistic: 58.757, p-value: 1.7875e-14
#> 
#> Log likelihood: -110.4248 for lag model
#> ML residual variance (sigma squared): 0.1589, (sigma: 0.39862)
#> Number of observations: 211 
#> Number of parameters estimated: 6 
#> AIC: 232.85, (AIC for lm: 283.51)
#> LM test for residual autocorrelation
#> test value: 8.7942, p-value: 0.0030219
#> 
# \donttest{
data(boston, package="spData")
lw <- spdep::nb2listw(boston.soi)
gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw, method="Matrix")
summary(gp2)
#> 
#> Call:lagsarlm(formula = log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = lw, method = "Matrix")
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.5262308 -0.0749699 -0.0044237  0.0713409  0.7122121 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>                Estimate  Std. Error  z value  Pr(>|z|)
#> (Intercept)  2.2796e+00  1.7495e-01  13.0302 < 2.2e-16
#> CRIM        -7.1045e-03  9.6236e-04  -7.3824 1.554e-13
#> ZN           3.7985e-04  3.8510e-04   0.9864 0.3239507
#> INDUS        1.2572e-03  1.7986e-03   0.6990 0.4845472
#> CHAS1        7.3677e-03  2.5416e-02   0.2899 0.7719057
#> I(NOX^2)    -2.6892e-01  8.8026e-02  -3.0550 0.0022508
#> I(RM^2)      6.7243e-03  1.0039e-03   6.6985 2.106e-11
#> AGE         -2.7682e-04  4.0062e-04  -0.6910 0.4895829
#> log(DIS)    -1.5830e-01  2.5554e-02  -6.1947 5.841e-10
#> log(RAD)     7.0689e-02  1.4616e-02   4.8363 1.323e-06
#> TAX         -3.6569e-04  9.3744e-05  -3.9009 9.582e-05
#> PTRATIO     -1.2011e-02  3.9599e-03  -3.0330 0.0024211
#> B            2.8432e-04  7.9402e-05   3.5807 0.0003427
#> log(LSTAT)  -2.3216e-01  2.0425e-02 -11.3663 < 2.2e-16
#> 
#> Rho: 0.48537, LR test value: 214.06, p-value: < 2.22e-16
#> Asymptotic standard error: 0.029426
#>     z-value: 16.494, p-value: < 2.22e-16
#> Wald statistic: 272.06, p-value: < 2.22e-16
#> 
#> Log likelihood: 264.0089 for lag model
#> ML residual variance (sigma squared): 0.019276, (sigma: 0.13884)
#> Number of observations: 506 
#> Number of parameters estimated: 16 
#> AIC: -496.02, (AIC for lm: -283.96)
#> LM test for residual autocorrelation
#> test value: 10.74, p-value: 0.0010486
#> 
gp2a <- lagmess(CMEDV ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw)
summary(gp2a)
#> Matrix exponential spatial lag model:
#> 
#> Call:
#> lagmess(formula = CMEDV ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = lw)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -17.03755  -2.05386  -0.30295   1.67710  21.82120 
#> 
#> Coefficients:
#>               Estimate Std. Error  t value  Pr(>|t|)
#> (Intercept) 35.3206850  3.0128759  11.7232 < 2.2e-16
#> CRIM        -0.0986056  0.0242819  -4.0609 5.688e-05
#> ZN           0.0198500  0.0098638   2.0124 0.0447222
#> INDUS        0.0071211  0.0461104   0.1544 0.8773301
#> CHAS1        0.8158059  0.6477224   1.2595 0.2084473
#> I(NOX^2)    -9.2592911  2.2072427  -4.1950 3.238e-05
#> I(RM^2)      0.2434745  0.0256001   9.5107 < 2.2e-16
#> AGE         -0.0040683  0.0102667  -0.3963 0.6920816
#> log(DIS)    -5.3974116  0.6514323  -8.2855 1.125e-15
#> log(RAD)     1.7142905  0.3732772   4.5925 5.569e-06
#> TAX         -0.0087053  0.0023933  -3.6373 0.0003046
#> PTRATIO     -0.4118524  0.0977997  -4.2112 3.021e-05
#> B            0.0056141  0.0020116   2.7908 0.0054614
#> log(LSTAT)  -6.1484203  0.4878957 -12.6019 < 2.2e-16
#> 
#> Residual standard error: 3.5594 on 492 degrees of freedom
#> Multiple R-squared:  0.76221,    Adjusted R-squared:  0.75593 
#> F-statistic: 121.31 on 13 and 492 DF,  p-value: < 2.22e-16
#> 
#> Alpha: -0.41361, standard error: 0.038521
#>     z-value: -10.737, p-value: < 2.22e-16
#> LR test value: 121.4, p-value: < 2.22e-16
#> Implied rho: 0.3387434 
#> 
# }
```
