# Semi-parametric spatial filtering

The function selects eigenvectors in a semi-parametric spatial filtering
approach to removing spatial dependence from linear models. Selection is
by brute force by finding the single eigenvector reducing the standard
variate of Moran's I for regression residuals most, and continuing until
no candidate eigenvector reduces the value by more than `tol`. It
returns a summary table from the selection process and a matrix of
selected eigenvectors for the specified model.

## Usage

``` r
SpatialFiltering(formula, lagformula=NULL, data=list(), na.action=na.fail,
 nb=NULL, glist = NULL,
 style = "C", zero.policy = NULL, tol = 0.1, zerovalue = 1e-04,
 ExactEV = FALSE, symmetric = TRUE, alpha=NULL, alternative="two.sided",
 verbose=NULL)
```

## Arguments

- formula:

  a symbolic description of the model to be fit, assuming a spatial
  error representation; when lagformula is given, it should include only
  the response and the intercept term

- lagformula:

  An extra one-sided formula to be used when a spatial lag
  representation is desired; the intercept is excluded within the
  function if present because it is part of the formula argument, but
  excluding it explicitly in the lagformula argument in the presence of
  factors generates a collinear model matrix

- data:

  an optional data frame containing the variables in the model

- nb:

  an object of class `nb`

- glist:

  list of general weights corresponding to neighbours

- style:

  `style` can take values W, B, C, U, and S

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the spatial weights list will be subsetted to remove NAs
  in the data. It may be necessary to set zero.policy to TRUE because
  this subsetting may create no-neighbour observations. Note that only
  weights lists created without using the glist argument to `nb2listw`
  may be subsetted.

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- tol:

  tolerance value for convergence of spatial filtering

- zerovalue:

  eigenvectors with eigenvalues of an absolute value smaller than
  zerovalue will be excluded in eigenvector search

- ExactEV:

  Set ExactEV=TRUE to use exact expectations and variances rather than
  the expectation and variance of Moran's I from the previous iteration,
  default FALSE

- symmetric:

  Should the spatial weights matrix be forced to symmetry, default TRUE

- alpha:

  if not NULL, used instead of the tol= argument as a stopping rule to
  choose all eigenvectors up to and including the one with a probability
  value exceeding alpha.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater, less or two.sided (default).

- verbose:

  default NULL, use global option value; if TRUE report eigenvectors
  selected

## Value

An `SfResult` object, with:

- selection:

  a matrix summarising the selection of eigenvectors for inclusion, with
  columns:

  Step

  :   Step counter of the selection procedure

  SelEvec

  :   number of selected eigenvector (sorted descending)

  Eval

  :   its associated eigenvalue

  MinMi

  :   value Moran's I for residual autocorrelation

  ZMinMi

  :   standardized value of Moran's I assuming a normal approximation

  pr(ZI)

  :   probability value of the permutation-based standardized deviate
      for the given value of the alternative argument

  R2

  :   R^2 of the model including exogenous variables and eigenvectors

  gamma

  :   regression coefficient of selected eigenvector in fit

  The first row is the value at the start of the search

- dataset:

  a matrix of the selected eigenvectors in order of selection

## References

Tiefelsdorf M, Griffith DA. (2007) Semiparametric Filtering of Spatial
Autocorrelation: The Eigenvector Approach. Environment and Planning A,
39 (5) 1193 - 1221.

## Author

Yongwan Chun, Michael Tiefelsdorf, Roger Bivand

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`eigen`](https://rdrr.io/r/base/eigen.html),
[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html),
[`listw2U`](https://r-spatial.github.io/spdep/reference/nb2listw.html)

## Examples

``` r
require("sf", quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
sarcol <- SpatialFiltering(CRIME ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", ExactEV=TRUE)
sarcol
#>   Step SelEvec      Eval        MinMi      ZMinMi      Pr(ZI)        R2
#> 0    0       0 0.0000000  0.212374153  2.68100025 0.007340246 0.5524040
#> 1    1       5 0.7148326  0.121528166  1.89037770 0.058707464 0.6209393
#> 2    2       3 0.8408661  0.065848648  1.54064108 0.123404165 0.6481722
#> 3    3       1 1.0206316 -0.005424824  1.08514557 0.277857187 0.6726114
#> 4    4      10 0.3658588 -0.039356232  0.80357070 0.421644951 0.7000258
#> 5    5      14 0.1831325 -0.072949543  0.47790213 0.632719864 0.7393770
#> 6    6      11 0.3144120 -0.108332631  0.18566599 0.852706701 0.7611907
#> 7    7       2 0.9157325 -0.153675621 -0.03464097 0.972366030 0.7713163
#>       gamma
#> 0   0.00000
#> 1  30.34786
#> 2  19.13010
#> 3 -18.12234
#> 4 -19.19379
#> 5  22.99586
#> 6  17.12127
#> 7  11.66487
lmsar <- lm(CRIME ~ INC + HOVAL + fitted(sarcol), data=columbus)
(x <- summary(lmsar))
#> 
#> Call:
#> lm(formula = CRIME ~ INC + HOVAL + fitted(sarcol), data = columbus)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -27.6527  -5.3084   0.0804   5.6844  15.6912 
#> 
#> Coefficients:
#>                      Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)          68.61896    3.67609  18.666  < 2e-16 ***
#> INC                  -1.59731    0.25938  -6.158 3.12e-07 ***
#> HOVAL                -0.27393    0.08011  -3.419  0.00148 ** 
#> fitted(sarcol)vec5   30.34786    8.87679   3.419  0.00149 ** 
#> fitted(sarcol)vec3   19.13010    8.87679   2.155  0.03739 *  
#> fitted(sarcol)vec1  -18.12234    8.87679  -2.042  0.04800 *  
#> fitted(sarcol)vec10 -19.19379    8.87679  -2.162  0.03679 *  
#> fitted(sarcol)vec14  22.99586    8.87679   2.591  0.01341 *  
#> fitted(sarcol)vec11  17.12127    8.87679   1.929  0.06106 .  
#> fitted(sarcol)vec2   11.66487    8.87679   1.314  0.19649    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 8.877 on 39 degrees of freedom
#> Multiple R-squared:  0.7713, Adjusted R-squared:  0.7185 
#> F-statistic: 14.62 on 9 and 39 DF,  p-value: 5.579e-10
#> 
coef(x)
#>                        Estimate Std. Error   t value     Pr(>|t|)
#> (Intercept)          68.6189611 3.67608656 18.666307 4.802352e-21
#> INC                  -1.5973108 0.25938068 -6.158172 3.122801e-07
#> HOVAL                -0.2739315 0.08011158 -3.419374 1.483698e-03
#> fitted(sarcol)vec5   30.3478552 8.87679493  3.418785 1.486164e-03
#> fitted(sarcol)vec3   19.1300996 8.87679493  2.155068 3.738943e-02
#> fitted(sarcol)vec1  -18.1223409 8.87679493 -2.041541 4.800339e-02
#> fitted(sarcol)vec10 -19.1937947 8.87679493 -2.162244 3.679422e-02
#> fitted(sarcol)vec14  22.9958588 8.87679493  2.590559 1.340783e-02
#> fitted(sarcol)vec11  17.1212741 8.87679493  1.928768 6.106079e-02
#> fitted(sarcol)vec2   11.6648669 8.87679493  1.314085 1.964945e-01
anova(lmbase, lmsar)
#> Analysis of Variance Table
#> 
#> Model 1: CRIME ~ INC + HOVAL
#> Model 2: CRIME ~ INC + HOVAL + fitted(sarcol)
#>   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
#> 1     46 6014.9                                  
#> 2     39 3073.1  7    2941.8 5.3334 0.0002445 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
spdep::lm.morantest(lmsar, spdep::nb2listw(col.gal.nb))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ INC + HOVAL + fitted(sarcol), data =
#> columbus)
#> weights: spdep::nb2listw(col.gal.nb)
#> 
#> Moran I statistic standard deviate = -0.034641, p-value = 0.5138
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>     -0.153675621     -0.150918131      0.006336477 
#> 
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL - 1, data=columbus,
 nb=col.gal.nb, style="W")
lagcol
#>    Step SelEvec      Eval       MinMi      ZMinMi      Pr(ZI)        R2
#> 0     0       0 0.0000000  0.21237415  2.68100025 0.007340246 0.5524040
#> 1     1       6 0.7161123  0.11782248  1.84511963 0.065020139 0.6038801
#> 2     2       4 0.8682938  0.06242664  1.49482111 0.134961136 0.6531288
#> 3     3       1 1.0310063 -0.02066604  0.88134183 0.378132834 0.6924845
#> 4     4       5 0.7905397 -0.04619973  0.84746904 0.396733736 0.7136578
#> 5     5      15 0.1753342 -0.07609524  0.55233191 0.580720971 0.7558543
#> 6     6       9 0.5501433 -0.10190889  0.43919419 0.660520837 0.7626784
#> 7     7       8 0.5721041 -0.12232942  0.41846803 0.675604953 0.7757314
#> 8     8       3 0.9026222 -0.14991822  0.38315383 0.701605709 0.7908693
#> 9     9       2 0.9649166 -0.21756342 -0.28556733 0.775209527 0.8078727
#> 10   10       7 0.6219404 -0.22017920 -0.04856547 0.961265592 0.8082842
#>         gamma
#> 0    0.000000
#> 1   19.848854
#> 2  -35.542595
#> 3  -30.697851
#> 4  -24.540372
#> 5   25.227798
#> 6    7.590082
#> 7  -16.933168
#> 8  -20.556931
#> 9  -18.434534
#> 10  -2.597572
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
lmlag
#> 
#> Call:
#> lm(formula = CRIME ~ INC + HOVAL + fitted(lagcol), data = columbus)
#> 
#> Coefficients:
#>         (Intercept)                  INC                HOVAL  
#>             56.7977              -0.4857              -0.3821  
#>  fitted(lagcol)vec6   fitted(lagcol)vec4   fitted(lagcol)vec1  
#>             19.8489             -35.5426             -30.6979  
#>  fitted(lagcol)vec5  fitted(lagcol)vec15   fitted(lagcol)vec9  
#>            -24.5404              25.2278               7.5901  
#>  fitted(lagcol)vec8   fitted(lagcol)vec3   fitted(lagcol)vec2  
#>            -16.9332             -20.5569             -18.4345  
#>  fitted(lagcol)vec7  
#>             -2.5976  
#> 
anova(lmbase, lmlag)
#> Analysis of Variance Table
#> 
#> Model 1: CRIME ~ INC + HOVAL
#> Model 2: CRIME ~ INC + HOVAL + fitted(lagcol)
#>   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
#> 1     46 6014.9                                  
#> 2     36 2576.3 10    3438.6 4.8049 0.0002165 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
spdep::lm.morantest(lmlag, spdep::nb2listw(col.gal.nb))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ INC + HOVAL + fitted(lagcol), data =
#> columbus)
#> weights: spdep::nb2listw(col.gal.nb)
#> 
#> Moran I statistic standard deviate = -0.048565, p-value = 0.5194
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>     -0.220179195     -0.217083975      0.004061888 
#> 
NA.columbus <- columbus
NA.columbus$CRIME[20:25] <- NA
COL.SF.NA <- SpatialFiltering(CRIME ~ INC + HOVAL, data=NA.columbus,
 nb=col.gal.nb, style="W", na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
COL.SF.NA$na.action
#> 20 21 22 23 24 25 
#> 20 21 22 23 24 25 
#> attr(,"class")
#> [1] "exclude"
summary(lm(CRIME ~ INC + HOVAL + fitted(COL.SF.NA), data=NA.columbus,
 na.action=na.exclude))
#> 
#> Call:
#> lm(formula = CRIME ~ INC + HOVAL + fitted(COL.SF.NA), data = NA.columbus, 
#>     na.action = na.exclude)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -23.6712  -4.7984   0.1761   6.7460  11.3353 
#> 
#> Coefficients:
#>                     Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)         69.04674    3.91643  17.630  < 2e-16 ***
#> INC                 -1.60115    0.26017  -6.154 3.88e-07 ***
#> HOVAL               -0.28742    0.07716  -3.725 0.000649 ***
#> fitted(COL.SF.NA)1 -39.91305    8.23597  -4.846 2.27e-05 ***
#> fitted(COL.SF.NA)2  19.81805    8.23597   2.406 0.021226 *  
#> fitted(COL.SF.NA)3 -35.07336    8.23597  -4.259 0.000135 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 8.236 on 37 degrees of freedom
#>   (6 observations deleted due to missingness)
#> Multiple R-squared:  0.7772, Adjusted R-squared:  0.7471 
#> F-statistic: 25.81 on 5 and 37 DF,  p-value: 3.996e-11
#> 
```
