# Likelihood ratio test

The `LR.Sarlm()` function provides a likelihood ratio test for objects
for which a [`logLik()`](https://rdrr.io/r/stats/logLik.html) function
exists for their class, or for objects of class `logLik`. `LR1.Sarlm()`
and `Wald1.Sarlm()` are used internally in
[`summary.Sarlm()`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
but may be accessed directly; they report the values respectively of LR
and Wald tests for the absence of spatial dependence in spatial lag or
error models. The spatial Hausman test is available for models fitted
with `errorSarlm` and `GMerrorsar`.

## Usage

``` r
LR.Sarlm(x, y)
# S3 method for class 'Sarlm'
logLik(object, ...)
LR1.Sarlm(object)
Wald1.Sarlm(object)
# S3 method for class 'Sarlm'
Hausman.test(object, ..., tol=NULL)
# S3 method for class 'Sarlm'
anova(object, ...)
bptest.Sarlm(object, varformula=NULL, studentize = TRUE, data=list())
# S3 method for class 'Sarlm'
impacts(obj, ..., tr, R = NULL, listw = NULL, evalues=NULL,
 useHESS = NULL, Q=NULL)
```

## Arguments

- x:

  a `logLik` object or an object for which a
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) function exists

- y:

  a `logLik` object or an object for which a
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) function exists

- object, obj:

  a `Sarlm` object

- ...:

  further arguments passed to or from other methods

- tol:

  `tol` argument passed to `solve`, default NULL

- varformula:

  a formula describing only the potential explanatory variables for the
  variance (no dependent variable needed). By default the same
  explanatory variables are taken as in the main regression model

- studentize:

  logical. If set to `TRUE` Koenker's studentized version of the test
  statistic will be used.

- data:

  an optional data frame containing the variables in the varformula

- tr:

  A vector of traces of powers of the spatial weights matrix created
  using `trW`, for approximate impact measures; if not given, `listw`
  must be given for exact measures (for small to moderate spatial
  weights matrices); the traces must be for the same spatial weights as
  were used in fitting the spatial regression, and must be
  row-standardised

- listw:

  If `tr` is not given, a spatial weights object as created by
  `nb2listw`; they must be the same spatial weights as were used in
  fitting the spatial regression, but do not have to be row-standardised

- evalues:

  vector of eigenvalues of spatial weights matrix for impacts
  calculations

- R:

  If given, simulations are used to compute distributions for the impact
  measures, returned as `mcmc` objects; the objects are used for
  convenience but are not output by an MCMC process

- useHESS:

  Use the Hessian approximation (if available) even if the asymptotic
  coefficient covariance matrix is available; used for comparing methods

- Q:

  default NULL, else an integer number of cumulative power series
  impacts to calculate if `tr` is given

## Value

The tests return objects of class `htest` with:

- statistic:

  value of statistic

- parameter:

  degrees of freedom

- p.value:

  Probability value

- estimate:

  varies with test

- method:

  description of test method

`logLik.Sarlm()` returns an object of class `logLik` `LR1.Sarlm`,
`Hausman.Sarlm` and `Wald1.Sarlm` returm objects of class `htest`

## Note

The numbers of degrees of freedom returned by `logLik.Sarlm()` include
nuisance parameters, that is the number of regression coefficients, plus
sigma, plus spatial parameter esitmate(s).

## References

LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC
Press, Boca Raton, pp. 61–63; Pace RK and LeSage J (2008) A spatial
Hausman test. *Economics Letters* 101, 282–284. T.S. Breusch & A.R.
Pagan (1979), A Simple Test for Heteroscedasticity and Random
Coefficient Variation. *Econometrica* **47**, 1287–1294

W. Krämer & H. Sonnberger (1986), *The Linear Regression Model under
Test*. Heidelberg: Physica.

L. Anselin (1988) *Spatial econometrics: methods and models.* Dordrecht:
Kluwer, pp. 121–122.

## Author

Roger Bivand <Roger.Bivand@nhh.no>, `bptest`: Torsten Hothorn and Achim
Zeileis, modified by Roger Bivand

## See also

[`logLik.lm`](https://rdrr.io/r/stats/logLik.html), `anova.Sarlm`,
[`impacts`](https://r-spatial.github.io/spatialreg/reference/impacts.md)

## Examples

``` r
require("sf", quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
lm.mod <- lm(CRIME ~ HOVAL + INC, data=columbus)
lag <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb))
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb), type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb))
Hausman.test(error)
#> 
#>  Spatial Hausman test (asymptotic)
#> 
#> data:  NULL
#> Hausman test = 6.4729, df = 3, p-value = 0.09074
#> 
LR.Sarlm(mixed, error)
#> 
#>  Likelihood ratio for spatial linear models
#> 
#> data:  
#> Likelihood ratio = 4.2782, df = 2, p-value = 0.1178
#> sample estimates:
#> Log likelihood of mixed Log likelihood of error 
#>               -182.0161               -184.1552 
#> 
anova(lag, lm.mod)
#>        Model df    AIC  logLik Test L.Ratio   p-value
#> lag        1  5 376.34 -183.17    1                  
#> lm.mod     2  4 382.75 -187.38    2  8.4179 0.0037154
anova(lag, error, mixed)
#>       Model df    AIC  logLik Test L.Ratio p-value
#> lag       1  5 376.34 -183.17    1                
#> error     2  5 378.31 -184.16    1                
#> mixed     3  7 378.03 -182.02    2  4.2782 0.11776
AIC(lag, error, mixed)
#>       df      AIC
#> lag    5 376.3366
#> error  5 378.3104
#> mixed  7 378.0322
bptest.Sarlm(error)
#> 
#>  studentized Breusch-Pagan test
#> 
#> data:  
#> BP = 9.3694, df = 2, p-value = 0.009235
#> 
bptest.Sarlm(error, studentize=FALSE)
#> 
#>  Breusch-Pagan test
#> 
#> data:  
#> BP = 16.285, df = 2, p-value = 0.0002908
#> 
```
