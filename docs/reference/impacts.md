# Impacts in spatial lag models

The calculation of impacts for spatial lag and spatial Durbin models is
needed in order to interpret the regression coefficients correctly,
because of the spillovers between the terms in these data generation
processes (unlike the spatial error model). Methods for “SLX” and
Bayesian fitted models are also provided, the former do not need MC
simulations, while the latter pass through MCMC draws.

## Usage

``` r
#\method{impacts}{sarlm}(obj, \dots, tr, R = NULL, listw = NULL, evalues=NULL,
# useHESS = NULL, Q=NULL)
#\method{impacts}{lagmess}(obj, ..., R=NULL, listw=NULL)
#\method{impacts}{SLX}(obj, ...)
#\method{impacts}{MCMC_sar_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
#\method{impacts}{MCMC_sem_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
#\method{impacts}{MCMC_sac_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
# S3 method for class 'LagImpact'
plot(x, ..., choice="direct", trace=FALSE, density=TRUE)
# S3 method for class 'LagImpact'
print(x, ..., reportQ=NULL)
# S3 method for class 'LagImpact'
summary(object, ..., zstats=FALSE, short=FALSE, reportQ=NULL)
#\method{print}{WXImpact}(x, ...)
#\method{summary}{WXImpact}(object, ..., adjust_k=(attr(object, "type") == "SDEM"))
# S3 method for class 'LagImpact'
HPDinterval(obj, prob = 0.95, ..., choice="direct")
intImpacts(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames, interval,
 type, tr, R, listw, evalues, tol, empirical, Q, icept, iicept, p, mess=FALSE,
 samples=NULL, zero_fill = NULL, dvars = NULL)
```

## Arguments

- obj:

  A spatial regression object created by `lagsarlm` or by `lmSLX`; in
  `HPDinterval.LagImpact`, a LagImpact object

- ...:

  Arguments passed through to methods in the coda package

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

- n:

  defaults to `length(obj$residuals)`; in the method for `gmsar` objects
  it may be used in panel settings to compute the impacts for
  cross-sectional weights only, suggested by Angela Parenti

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

- reportQ:

  default NULL; if TRUE and `Q` given as an argument to `impacts`,
  report impact components

- x, object:

  LagImpact objects created by `impacts` methods

- zstats:

  default FALSE, if TRUE, also return z-values and p-values for the
  impacts based on the simulations

- short:

  default FALSE, if TRUE passed to the print summary method to omit
  printing of the mcmc summaries

- choice:

  One of three impacts: direct, indirect, or total

- trace:

  Argument passed to `plot.mcmc`: plot trace plots

- density:

  Argument passed to `plot.mcmc`: plot density plots

- prob:

  Argument passed to `HPDinterval.mcmc`: a numeric scalar in the
  interval (0,1) giving the target probability content of the intervals

- adjust_k:

  default TRUE if SDEM else FALSE, adjust internal OLS SDEM standard
  errors by dividing by n rather than (n-k) (default changed and bug
  fixed after 0.7-8; standard errors now ML in SDEM summary and impacts
  summary and identical - for SLX use FALSE)

- rho, beta, P, mu, Sigma, irho, drop2beta, bnames, interval, type,
  icept, iicept, p, mess, samples, zero_fill, dvars:

  internal arguments shared inside impacts methods

- tol, empirical:

  deprecated arguments

## Details

If called without `R` being set, the method returns the direct, indirect
and total impacts for the variables in the model, for the variables
themselves in tha spatial lag model case, for the variables and their
spatial lags in the spatial Durbin (mixed) model case. The spatial lag
impact measures are computed using eq. 2.46 (LeSage and Pace, 2009, p.
38), either using the exact dense matrix (when `listw` is given), or
traces of powers of the weights matrix (when `tr` is given). When the
traces are created by powering sparse matrices, the exact and the trace
methods should give very similar results, unless the number of powers
used is very small, or the spatial coefficient is close to its bounds.

If `R` is given, simulations will be used to create distributions for
the impact measures, provided that the fitted model object contains a
coefficient covariance matrix. The simulations are made using
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) with the
coefficients and their covariance matrix from the fitted model.

The simulations are stored as `mcmc` objects as defined in the coda
package; the objects are used for convenience but are not output by an
MCMC process. The simulated values of the coefficients are checked to
see that the spatial coefficient remains within its valid interval —
draws outside the interval are discarded.

If a model is fitted with the “Durbin=” set to a formula subsetting the
explanatory variables, the impacts object returned reports Durbin
impacts for variables included in the formula and lag impacts for the
other variables.

When `Q` and `tr` are given, addition impact component results are
provided for each step in the traces of powers of the weights matrix up
to and including the `Q`'th power. This increases computing time because
the output object is substantially increased in size in proportion to
the size of `Q`.

The method for `gmsar` objects is only for those of `type` `SARAR`
output by `gstsls`, and assume that the spatial error coefficient is
fixed, and thus omitted from the coefficients and covariance matrix used
for simulation.

From version 1.4.1, functions for models including spatially lagged
independent variables warn on fitting if any of the right-hand side
variables are factors. This is because the interpretation of
coefficients that are not slopes is unclear when the variable is not
interpretable on an unbounded line, such as factors. Factor variable
names are shown with the suffix “(F)”, others “dy/dx” in output from
impact methods. A discussion can be found at
<https://github.com/rsbivand/eqc25_talk>.

## Value

An object of class LagImpact.

If no simulation is carried out, the object returned is a list with:

- direct:

  numeric vector

- indirect:

  numeric vector

- total:

  numeric vector

and a matching `Qres` list attribute if `Q` was given.

If simulation is carried out, the object returned is a list with:

- res:

  a list with three components as for the non-simulation case, with a
  matching `Qres` list attribute if `Q` was given

- sres:

  a list with three `mcmc` matrices, for the direct, indirect and total
  impacts with a matching `Qmcmc` list attribute if `Q` was given

## References

LeSage J and RK Pace (2009) *Introduction to Spatial Econometrics*. CRC
Press, Boca Raton, pp. 33–42, 114–115; LeSage J and MM Fischer (2008)
Spatial growth regressions: model specification, estimation and
interpretation. *Spatial Economic Analysis* 3 (3), pp. 275–304.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of
Estimation Methods for Spatial Econometrics. *Journal of Statistical
Software*, 63(18), 1-36.
[doi:10.18637/jss.v063.i18](https://doi.org/10.18637/jss.v063.i18) .

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`trW`](https://r-spatial.github.io/spatialreg/reference/trW.md),
[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html),
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html),
[`plot.mcmc`](https://rdrr.io/pkg/coda/man/plot.mcmc.html),
[`summary.mcmc`](https://rdrr.io/pkg/coda/man/summary.mcmc.html),
[`HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html)

## Examples

``` r
require("sf", quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
columbus$fEW <- factor(columbus$EW)
columbus$fDISCBD <- ordered(cut(columbus$DISCBD, c(0, 1.5, 3, 4.5, 6)))
run <- require("codingMatrices", quietly=TRUE)
f <- formula(log(CRIME) ~ INC + HOVAL + fDISCBD + fEW)
listw <- spdep::nb2listw(col.gal.nb)
ev <- eigenw(listw)
lobj <- lagsarlm(f, columbus, listw, control=list(pre_eig=ev))
summary(lobj)
#> 
#> Call:
#> lagsarlm(formula = f, data = columbus, listw = listw, control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -3.887021 -0.265705  0.075038  0.351142  1.682544 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>               Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)  6.1392333  0.7750510  7.9211 2.442e-15
#> INC         -0.0520746  0.0297726 -1.7491   0.08028
#> HOVAL       -0.0290994  0.0077396 -3.7598   0.00017
#> fDISCBD.L   -0.1168456  0.3981785 -0.2935   0.76918
#> fDISCBD.Q    0.3670863  0.2549784  1.4397   0.14996
#> fDISCBD.C   -0.0370120  0.2296689 -0.1612   0.87197
#> fEW1         0.0837988  0.2460444  0.3406   0.73342
#> 
#> Rho: -0.30114, LR test value: 2.2103, p-value: 0.13709
#> Asymptotic standard error: 0.19387
#>     z-value: -1.5534, p-value: 0.12034
#> Wald statistic: 2.4129, p-value: 0.12034
#> 
#> Log likelihood: -60.69996 for lag model
#> ML residual variance (sigma squared): 0.68415, (sigma: 0.82714)
#> Number of observations: 49 
#> Number of parameters estimated: 9 
#> AIC: 139.4, (AIC for lm: 139.61)
#> LM test for residual autocorrelation
#> test value: 1.1851, p-value: 0.27632
#> 
if (run) {
contrasts(columbus$fDISCBD) <- "code_diff"
lobjd <- lagsarlm(f, columbus, listw, control=list(pre_eig=ev))
}
mobj <- lagsarlm(f, columbus, listw, Durbin=TRUE, control=list(pre_eig=ev))
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(mobj)
#> 
#> Call:lagsarlm(formula = f, data = columbus, listw = listw, Durbin = TRUE, 
#>     control = list(pre_eig = ev))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -3.485440 -0.303568  0.057651  0.412917  1.618395 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                    Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)       7.2505614  1.5705214  4.6167   3.9e-06
#> INC              -0.0511360  0.0317741 -1.6094 0.1075380
#> HOVAL            -0.0284100  0.0078672 -3.6112 0.0003048
#> fDISCBDm2-m1      0.2649701  0.4829325  0.5487 0.5832325
#> fDISCBDm3-m2      0.2832158  0.4781799  0.5923 0.5536639
#> fDISCBDm4-m3      0.5248924  0.5018269  1.0460 0.2955780
#> fEW1             -0.0391368  0.5485068 -0.0714 0.9431180
#> lag.INC          -0.0340158  0.0626300 -0.5431 0.5870447
#> lag.HOVAL        -0.0056745  0.0174504 -0.3252 0.7450467
#> lag.fDISCBDm2-m1 -1.7168592  0.7865735 -2.1827 0.0290574
#> lag.fDISCBDm3-m2  0.2492410  0.7088775  0.3516 0.7251386
#> lag.fDISCBDm4-m3 -0.3272604  1.1177652 -0.2928 0.7696896
#> lag.fEW1          0.2155133  0.6636093  0.3248 0.7453632
#> 
#> Rho: -0.43661, LR test value: 3.3537, p-value: 0.067054
#> Asymptotic standard error: 0.21519
#>     z-value: -2.0289, p-value: 0.042467
#> Wald statistic: 4.1165, p-value: 0.042467
#> 
#> Log likelihood: -58.21561 for mixed model
#> ML residual variance (sigma squared): 0.60545, (sigma: 0.7781)
#> Number of observations: 49 
#> Number of parameters estimated: 15 
#> AIC: 146.43, (AIC for lm: 147.78)
#> LM test for residual autocorrelation
#> test value: 5.1201, p-value: 0.02365
#> 
mobj1 <- lagsarlm(f, columbus, listw, Durbin= ~ INC + HOVAL, control=list(pre_eig=ev))
summary(mobj1)
#> 
#> Call:lagsarlm(formula = f, data = columbus, listw = listw, Durbin = ~INC + 
#>     HOVAL, control = list(pre_eig = ev))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.86709 -0.26174  0.01044  0.33923  1.70763 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)   6.6778949  1.3996025  4.7713 1.831e-06
#> INC          -0.0512618  0.0303927 -1.6866 0.0916708
#> HOVAL        -0.0294172  0.0079682 -3.6918 0.0002227
#> fDISCBDm2-m1 -0.4432283  0.3456691 -1.2822 0.1997610
#> fDISCBDm3-m2  0.0682508  0.4063270  0.1680 0.8666068
#> fDISCBDm4-m3  0.3660184  0.4599471  0.7958 0.4261577
#> fEW1          0.1093682  0.2634028  0.4152 0.6779863
#> lag.INC      -0.0163520  0.0587022 -0.2786 0.7805832
#> lag.HOVAL    -0.0036691  0.0165743 -0.2214 0.8248025
#> 
#> Rho: -0.35002, LR test value: 2.1961, p-value: 0.13836
#> Asymptotic standard error: 0.21643
#>     z-value: -1.6173, p-value: 0.10582
#> Wald statistic: 2.6155, p-value: 0.10582
#> 
#> Log likelihood: -60.61159 for mixed model
#> ML residual variance (sigma squared): 0.67718, (sigma: 0.82291)
#> Number of observations: 49 
#> Number of parameters estimated: 11 
#> AIC: 143.22, (AIC for lm: 143.42)
#> LM test for residual autocorrelation
#> test value: 4.3874, p-value: 0.036205
#> 
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
set.seed(1)
impacts(lobj, listw=listw)
#> Impact measures (lag, exact):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05306202  0.013039791 -0.04002223
#> HOVAL dy/dx   -0.02965119  0.007286669 -0.02236453
#> fDISCBD.L (F) -0.11906113  0.029258822 -0.08980231
#> fDISCBD.Q (F)  0.37404673 -0.091920565  0.28212616
#> fDISCBD.C (F) -0.03771380  0.009268022 -0.02844577
#> fEW1 (F)       0.08538769 -0.020983702  0.06440399
if (run) {
impacts(lobjd, listw=listw)
}
#> Impact measures (lag, exact):
#>                        Direct      Indirect        Total
#> INC dy/dx        -0.053062024  0.0130397910 -0.040022233
#> HOVAL dy/dx      -0.029651195  0.0072866685 -0.022364526
#> fDISCBDm2-m1 (F) -0.461024722  0.1132950753 -0.347729647
#> fDISCBDm3-m2 (F) -0.002647386  0.0006505851 -0.001996801
#> fDISCBDm4-m3 (F)  0.287068727 -0.0705460500  0.216522677
#> fEW1 (F)          0.085387695 -0.0209837019  0.064403993
impacts(lobj, tr=trMatc)
#> Impact measures (lag, trace):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05306202  0.013039791 -0.04002223
#> HOVAL dy/dx   -0.02965119  0.007286669 -0.02236453
#> fDISCBD.L (F) -0.11906113  0.029258822 -0.08980231
#> fDISCBD.Q (F)  0.37404673 -0.091920565  0.28212616
#> fDISCBD.C (F) -0.03771380  0.009268022 -0.02844577
#> fEW1 (F)       0.08538769 -0.020983702  0.06440399
impacts(lobj, tr=trMC)
#> Impact measures (lag, trace):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05303288  0.013010644 -0.04002223
#> HOVAL dy/dx   -0.02963491  0.007270381 -0.02236453
#> fDISCBD.L (F) -0.11899573  0.029193421 -0.08980231
#> fDISCBD.Q (F)  0.37384126 -0.091715099  0.28212616
#> fDISCBD.C (F) -0.03769308  0.009247306 -0.02844577
#> fEW1 (F)       0.08534079 -0.020936798  0.06440399
impacts(lobj, evalues=ev)
#> Impact measures (lag, evalues):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05306202  0.013039791 -0.04002223
#> HOVAL dy/dx   -0.02965119  0.007286669 -0.02236453
#> fDISCBD.L (F) -0.11906113  0.029258822 -0.08980231
#> fDISCBD.Q (F)  0.37404673 -0.091920565  0.28212616
#> fDISCBD.C (F) -0.03771380  0.009268022 -0.02844577
#> fEW1 (F)       0.08538769 -0.020983702  0.06440399
library(coda)
lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)
summary(lobjIQ5, zstats=TRUE, short=TRUE)
#> Impact measures (lag, trace):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05306202  0.013039791 -0.04002223
#> HOVAL dy/dx   -0.02965119  0.007286669 -0.02236453
#> fDISCBD.L (F) -0.11906113  0.029258822 -0.08980231
#> fDISCBD.Q (F)  0.37404673 -0.091920565  0.28212616
#> fDISCBD.C (F) -0.03771380  0.009268022 -0.02844577
#> fEW1 (F)       0.08538769 -0.020983702  0.06440399
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                    Direct    Indirect       Total
#> INC dy/dx     0.032198021 0.010586529 0.027343925
#> HOVAL dy/dx   0.008536286 0.004965342 0.007295719
#> fDISCBD.L (F) 0.427784454 0.124236614 0.327396547
#> fDISCBD.Q (F) 0.272884419 0.088311572 0.226247266
#> fDISCBD.C (F) 0.242343959 0.063144429 0.196654588
#> fEW1 (F)      0.282959370 0.070702863 0.227608861
#> 
#> Simulated z-values:
#>                   Direct   Indirect      Total
#> INC dy/dx     -1.6136470  1.0027299 -1.5118829
#> HOVAL dy/dx   -3.5534014  1.3533361 -3.2365657
#> fDISCBD.L (F) -0.3012284  0.4339303 -0.2289297
#> fDISCBD.Q (F)  1.3345289 -0.8796413  1.2662679
#> fDISCBD.C (F) -0.1840497  0.1315488 -0.1845712
#> fEW1 (F)       0.2853271 -0.1160991  0.3186494
#> 
#> Simulated p-values:
#>               Direct     Indirect Total    
#> INC dy/dx     0.10660401 0.31599  0.1305637
#> HOVAL dy/dx   0.00038028 0.17595  0.0012098
#> fDISCBD.L (F) 0.76324036 0.66434  0.8189236
#> fDISCBD.Q (F) 0.18203057 0.37905  0.2054172
#> fDISCBD.C (F) 0.85397444 0.89534  0.8535654
#> fEW1 (F)      0.77539357 0.90757  0.7499924
summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
#> Impact measures (lag, trace):
#>                    Direct     Indirect       Total
#> INC dy/dx     -0.05306202  0.013039791 -0.04002223
#> HOVAL dy/dx   -0.02965119  0.007286669 -0.02236453
#> fDISCBD.L (F) -0.11906113  0.029258822 -0.08980231
#> fDISCBD.Q (F)  0.37404673 -0.091920565  0.28212616
#> fDISCBD.C (F) -0.03771380  0.009268022 -0.02844577
#> fEW1 (F)       0.08538769 -0.020983702  0.06440399
#> =================================
#> Impact components
#> $direct
#>        INC dy/dx   HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)
#> Q1 -5.207463e-02 -2.909943e-02 -0.1168455984  0.3670863248 -3.701200e-02
#> Q2  0.000000e+00  0.000000e+00  0.0000000000  0.0000000000  0.000000e+00
#> Q3 -1.051311e-03 -5.874752e-04 -0.0023589426  0.0074109386 -7.472185e-04
#> Q4  1.059547e-04  5.920777e-05  0.0002377423 -0.0007468999  7.530725e-05
#> Q5 -4.911361e-05 -2.744481e-05 -0.0001102016  0.0003462134 -3.490746e-05
#>         fEW1 (F)
#> Q1  8.379877e-02
#> Q2  0.000000e+00
#> Q3  1.691775e-03
#> Q4 -1.705029e-04
#> Q5  7.903388e-05
#> 
#> $indirect
#>        INC dy/dx   HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)
#> Q1  0.0000000000  0.0000000000  0.0000000000   0.000000000  0.0000000000
#> Q2  0.0156818807  0.0087630750  0.0351871702  -0.110545276  0.0111458856
#> Q3 -0.0036711691 -0.0020514587 -0.0082374082   0.025878937 -0.0026092809
#> Q4  0.0013161845  0.0007354873  0.0029532687  -0.009278095  0.0009354772
#> Q5 -0.0003791529 -0.0002118716 -0.0008507473   0.002672738 -0.0002694826
#>         fEW1 (F)
#> Q1  0.0000000000
#> Q2 -0.0252353665
#> Q3  0.0059076650
#> Q4 -0.0021180111
#> Q5  0.0006101349
#> 
#> $total
#>        INC dy/dx   HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)
#> Q1 -0.0520746262 -0.0290994341 -0.1168455984   0.367086325 -0.0370120037
#> Q2  0.0156818807  0.0087630750  0.0351871702  -0.110545276  0.0111458856
#> Q3 -0.0047224800 -0.0026389339 -0.0105963508   0.033289876 -0.0033564993
#> Q4  0.0014221392  0.0007946950  0.0031910111  -0.010024995  0.0010107844
#> Q5 -0.0004282665 -0.0002393164 -0.0009609489   0.003018952 -0.0003043901
#>         fEW1 (F)
#> Q1  0.0837987677
#> Q2 -0.0252353665
#> Q3  0.0075994402
#> Q4 -0.0022885141
#> Q5  0.0006891687
#> 
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                    Direct    Indirect       Total
#> INC dy/dx     0.032198021 0.010586529 0.027343925
#> HOVAL dy/dx   0.008536286 0.004965342 0.007295719
#> fDISCBD.L (F) 0.427784454 0.124236614 0.327396547
#> fDISCBD.Q (F) 0.272884419 0.088311572 0.226247266
#> fDISCBD.C (F) 0.242343959 0.063144429 0.196654588
#> fEW1 (F)      0.282959370 0.070702863 0.227608861
#> 
#> Simulated z-values:
#>                   Direct   Indirect      Total
#> INC dy/dx     -1.6136470  1.0027299 -1.5118829
#> HOVAL dy/dx   -3.5534014  1.3533361 -3.2365657
#> fDISCBD.L (F) -0.3012284  0.4339303 -0.2289297
#> fDISCBD.Q (F)  1.3345289 -0.8796413  1.2662679
#> fDISCBD.C (F) -0.1840497  0.1315488 -0.1845712
#> fEW1 (F)       0.2853271 -0.1160991  0.3186494
#> 
#> Simulated p-values:
#>               Direct     Indirect Total    
#> INC dy/dx     0.10660401 0.31599  0.1305637
#> HOVAL dy/dx   0.00038028 0.17595  0.0012098
#> fDISCBD.L (F) 0.76324036 0.66434  0.8189236
#> fDISCBD.Q (F) 0.18203057 0.37905  0.2054172
#> fDISCBD.C (F) 0.85397444 0.89534  0.8535654
#> fEW1 (F)      0.77539357 0.90757  0.7499924
#> ========================================================
#> Simulated impact components z-values:
#> $Direct
#>     INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)    fEW1 (F)
#> Q1 -1.6068710  -3.5684081    -0.2920634     1.3343061   -0.18490469  0.28995270
#> Q2        NaN         NaN           NaN           NaN           NaN         NaN
#> Q3 -0.6654716  -0.8392850    -0.3689039     0.5729230   -0.09080301  0.05026793
#> Q4  0.4283678   0.5529226     0.2840155    -0.3627465    0.07612433  0.01148692
#> Q5 -0.3003730  -0.3915966    -0.2212287     0.2570945   -0.06566033 -0.04834433
#> 
#> $Indirect
#>     INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)    fEW1 (F)
#> Q1        NaN         NaN           NaN           NaN           NaN         NaN
#> Q2  0.9836815   1.3049540     0.4359765    -0.8644240    0.12534191 -0.11183678
#> Q3 -0.6654716  -0.8392850    -0.3689039     0.5729230   -0.09080301  0.05026793
#> Q4  0.4283678   0.5529226     0.2840155    -0.3627465    0.07612433  0.01148692
#> Q5 -0.3003730  -0.3915966    -0.2212287     0.2570945   -0.06566033 -0.04834433
#> 
#> $Total
#>     INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F)    fEW1 (F)
#> Q1 -1.6068710  -3.5684081    -0.2920634     1.3343061   -0.18490469  0.28995270
#> Q2  0.9836815   1.3049540     0.4359765    -0.8644240    0.12534191 -0.11183678
#> Q3 -0.6654716  -0.8392850    -0.3689039     0.5729230   -0.09080301  0.05026793
#> Q4  0.4283678   0.5529226     0.2840155    -0.3627465    0.07612433  0.01148692
#> Q5 -0.3003730  -0.3915966    -0.2212287     0.2570945   -0.06566033 -0.04834433
#> 
#> 
#> Simulated impact components p-values:
#> $Direct
#>    INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F) fEW1 (F)
#> Q1 0.10808   0.00035916  0.77024       0.18210       0.85330       0.77185 
#> Q2 NA        NA          NA            NA            NA            NA      
#> Q3 0.50575   0.40130938  0.71220       0.56670       0.92765       0.95991 
#> Q4 0.66838   0.58031642  0.77640       0.71679       0.93932       0.99083 
#> Q5 0.76389   0.69535631  0.82491       0.79711       0.94765       0.96144 
#> 
#> $Indirect
#>    INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F) fEW1 (F)
#> Q1 NA        NA          NA            NA            NA            NA      
#> Q2 0.32527   0.19191     0.66285       0.38736       0.90025       0.91095 
#> Q3 0.50575   0.40131     0.71220       0.56670       0.92765       0.95991 
#> Q4 0.66838   0.58032     0.77640       0.71679       0.93932       0.99083 
#> Q5 0.76389   0.69536     0.82491       0.79711       0.94765       0.96144 
#> 
#> $Total
#>    INC dy/dx HOVAL dy/dx fDISCBD.L (F) fDISCBD.Q (F) fDISCBD.C (F) fEW1 (F)
#> Q1 0.10808   0.00035916  0.77024       0.18210       0.85330       0.77185 
#> Q2 0.32527   0.19190852  0.66285       0.38736       0.90025       0.91095 
#> Q3 0.50575   0.40130938  0.71220       0.56670       0.92765       0.95991 
#> Q4 0.66838   0.58031642  0.77640       0.71679       0.93932       0.99083 
#> Q5 0.76389   0.69535631  0.82491       0.79711       0.94765       0.96144 
#> 
impacts(mobj, listw=listw)
#> Impact measures (mixed, exact):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136498 -1.442000098 -1.01063512
#> fDISCBDm3-m2 (F)  0.27181264  0.098821437  0.37063408
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848169  0.13756830
#> fEW1 (F)         -0.06025686  0.183029567  0.12277271
impacts(mobj, tr=trMatc)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136498 -1.442000098 -1.01063512
#> fDISCBDm3-m2 (F)  0.27181264  0.098821437  0.37063408
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848169  0.13756830
#> fEW1 (F)         -0.06025686  0.183029567  0.12277271
impacts(mobj, tr=trMC)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05011897 -0.009153760 -0.05927273
#> HOVAL dy/dx      -0.02899549  0.005269884 -0.02372561
#> fDISCBDm2-m1 (F)  0.42440866 -1.435043779 -1.01063512
#> fDISCBDm3-m2 (F)  0.27228937  0.098344715  0.37063408
#> fDISCBDm4-m3 (F)  0.57330426 -0.435735956  0.13756830
#> fEW1 (F)         -0.05937391  0.182146618  0.12277271
impacts(mobj1, tr=trMatc)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05137786  0.001294350 -0.05008351
#> HOVAL dy/dx      -0.02990055  0.005392585 -0.02450797
#> fDISCBDm2-m1 (F) -0.45454225  0.126230607 -0.32831165
#> fDISCBDm3-m2 (F)  0.06999298 -0.019437700  0.05055528
#> fDISCBDm4-m3 (F)  0.37536154 -0.104241388  0.27112015
#> fEW1 (F)          0.11215995 -0.031147860  0.08101209
impacts(mobj1, listw=listw)
#> Impact measures (mixed, exact):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05137786  0.001294350 -0.05008351
#> HOVAL dy/dx      -0.02990055  0.005392585 -0.02450797
#> fDISCBDm2-m1 (F) -0.45454225  0.126230607 -0.32831165
#> fDISCBDm3-m2 (F)  0.06999298 -0.019437700  0.05055528
#> fDISCBDm4-m3 (F)  0.37536154 -0.104241388  0.27112015
#> fEW1 (F)          0.11215995 -0.031147860  0.08101209
# \dontrun{
try(impacts(mobj, evalues=ev), silent=TRUE)
#> Impact measures (mixed, evalues):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136498 -1.442000098 -1.01063512
#> fDISCBDm3-m2 (F)  0.27181264  0.098821437  0.37063408
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848169  0.13756830
#> fEW1 (F)         -0.06025686  0.183029567  0.12277271
# }
summary(impacts(mobj, tr=trMatc, R=200), short=TRUE, zstats=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136498 -1.442000098 -1.01063512
#> fDISCBDm3-m2 (F)  0.27181264  0.098821437  0.37063408
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848169  0.13756830
#> fEW1 (F)         -0.06025686  0.183029567  0.12277271
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                      Direct  Indirect     Total
#> INC dy/dx        0.04007337 0.4246233 0.4395851
#> HOVAL dy/dx      0.01147091 0.1760879 0.1813305
#> fDISCBDm2-m1 (F) 0.58904879 6.9113761 7.0433550
#> fDISCBDm3-m2 (F) 0.59974853 3.3448213 3.3831617
#> fDISCBDm4-m3 (F) 0.60661743 1.2316474 1.0088726
#> fEW1 (F)         0.67401119 0.7694596 0.3121581
#> 
#> Simulated z-values:
#>                       Direct    Indirect       Total
#> INC dy/dx        -1.29161747  0.04169967 -0.07746578
#> HOVAL dy/dx      -2.51461026  0.09503796 -0.06678320
#> fDISCBDm2-m1 (F)  0.84297356 -0.14998260 -0.07667278
#> fDISCBDm3-m2 (F)  0.41106515 -0.02024585  0.05285498
#> fDISCBDm4-m3 (F)  1.06912710 -0.31239394  0.26147200
#> fEW1 (F)         -0.08051817  0.23472324  0.40473050
#> 
#> Simulated p-values:
#>                  Direct   Indirect Total  
#> INC dy/dx        0.196490 0.96674  0.93825
#> HOVAL dy/dx      0.011916 0.92428  0.94675
#> fDISCBDm2-m1 (F) 0.399243 0.88078  0.93888
#> fDISCBDm3-m2 (F) 0.681025 0.98385  0.95785
#> fDISCBDm4-m3 (F) 0.285012 0.75474  0.79373
#> fEW1 (F)         0.935825 0.81442  0.68568
summary(impacts(mobj1, tr=trMatc, R=200), short=TRUE, zstats=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05137786  0.001294350 -0.05008351
#> HOVAL dy/dx      -0.02990055  0.005392585 -0.02450797
#> fDISCBDm2-m1 (F) -0.45454225  0.126230607 -0.32831165
#> fDISCBDm3-m2 (F)  0.06999298 -0.019437700  0.05055528
#> fDISCBDm4-m3 (F)  0.37536154 -0.104241388  0.27112015
#> fEW1 (F)          0.11215995 -0.031147860  0.08101209
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                       Direct   Indirect      Total
#> INC dy/dx        0.030366547 0.05056312 0.04476444
#> HOVAL dy/dx      0.008890799 0.01437362 0.01275784
#> fDISCBDm2-m1 (F) 0.340527773 0.13712780 0.25110468
#> fDISCBDm3-m2 (F) 0.410061367 0.12934391 0.30174660
#> fDISCBDm4-m3 (F) 0.451945152 0.14404190 0.34013015
#> fEW1 (F)         0.260039301 0.07482156 0.20105110
#> 
#> Simulated z-values:
#>                       Direct    Indirect       Total
#> INC dy/dx        -1.66905591 -0.01740920 -1.15189036
#> HOVAL dy/dx      -3.36304376  0.51780714 -1.76028090
#> fDISCBDm2-m1 (F) -1.28019493  0.86393144 -1.26430502
#> fDISCBDm3-m2 (F)  0.07577798 -0.05188593  0.08073825
#> fDISCBDm4-m3 (F)  0.83601090 -0.71734836  0.80705238
#> fEW1 (F)          0.37971460 -0.38196374  0.34897394
#> 
#> Simulated p-values:
#>                  Direct     Indirect Total  
#> INC dy/dx        0.09510630 0.98611  0.24937
#> HOVAL dy/dx      0.00077088 0.60459  0.07836
#> fDISCBDm2-m1 (F) 0.20047659 0.38763  0.20612
#> fDISCBDm3-m2 (F) 0.93959573 0.95862  0.93565
#> fDISCBDm4-m3 (F) 0.40314877 0.47316  0.41964
#> fEW1 (F)         0.70415728 0.70249  0.72711
xobj <- lmSLX(f, columbus, listw)
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(impacts(xobj))
#> Impact measures (SlX, glht, n-k):
#>                       Direct    Indirect       Total
#> INC dy/dx        -0.05127030 -0.01300677 -0.06427708
#> HOVAL dy/dx      -0.02823116  0.00776633 -0.02046483
#> fDISCBDm2-m1 (F)  0.30550205 -1.48251556 -1.17701351
#> fDISCBDm3-m2 (F)  0.19545891  0.36450618  0.55996509
#> fDISCBDm4-m3 (F)  0.52090107 -0.56818138 -0.04728031
#> fEW1 (F)         -0.05759815  0.25527403  0.19767588
#> ========================================================
#> Standard errors:
#>                       Direct   Indirect      Total
#> INC dy/dx        0.039113827 0.07595938 0.08320131
#> HOVAL dy/dx      0.009684776 0.02041006 0.01996799
#> fDISCBDm2-m1 (F) 0.589593465 0.96642237 0.67596041
#> fDISCBDm3-m2 (F) 0.588701391 0.87056187 0.75333500
#> fDISCBDm4-m3 (F) 0.617427646 1.37625313 1.23485161
#> fEW1 (F)         0.675512118 0.81732233 0.38451537
#> ========================================================
#> Z-values:
#>                      Direct   Indirect       Total
#> INC dy/dx        -1.3107974 -0.1712333 -0.77254884
#> HOVAL dy/dx      -2.9150036  0.3805147 -1.02488147
#> fDISCBDm2-m1 (F)  0.5181571 -1.5340245 -1.74124623
#> fDISCBDm3-m2 (F)  0.3320171  0.4187022  0.74331484
#> fDISCBDm4-m3 (F)  0.8436633 -0.4128466 -0.03828825
#> fEW1 (F)         -0.0852659  0.3123297  0.51409098
#> 
#> p-values:
#>                  Direct    Indirect Total  
#> INC dy/dx        0.1899262 0.86404  0.43979
#> HOVAL dy/dx      0.0035568 0.70356  0.30542
#> fDISCBDm2-m1 (F) 0.6043487 0.12502  0.08164
#> fDISCBDm3-m2 (F) 0.7398764 0.67543  0.45729
#> fDISCBDm4-m3 (F) 0.3988576 0.67972  0.96946
#> fEW1 (F)         0.9320500 0.75479  0.60719
#> 
eobj <- errorsarlm(f, columbus, listw, etype="emixed")
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(impacts(eobj), adjust_k=TRUE)
#> Impact measures (SDEM, glht, n):
#>                         Direct     Indirect       Total
#> INC dy/dx          -0.05700667 -0.014643459 -0.07165013
#> HOVAL dy/dx        -0.02851468  0.002480578 -0.02603410
#> fDISCBDm2-m1 dy/dx  0.53615183 -1.821801604 -1.28564977
#> fDISCBDm3-m2 dy/dx  0.44837745 -0.008451605  0.43992584
#> fDISCBDm4-m3 dy/dx  0.57159014 -0.134565014  0.43702512
#> fEW1 dy/dx         -0.09418651  0.247893670  0.15370716
#> ========================================================
#> Standard errors:
#>                         Direct   Indirect      Total
#> INC dy/dx          0.033522267 0.05500753 0.05177221
#> HOVAL dy/dx        0.008460333 0.01570895 0.01266713
#> fDISCBDm2-m1 dy/dx 0.526900798 0.80321952 0.43061797
#> fDISCBDm3-m2 dy/dx 0.525906684 0.70085294 0.45630345
#> fDISCBDm4-m3 dy/dx 0.535329561 1.11949891 0.87187542
#> fEW1 dy/dx         0.579445970 0.67806872 0.22642612
#> ========================================================
#> Z-values:
#>                        Direct    Indirect      Total
#> INC dy/dx          -1.7005614 -0.26620825 -1.3839496
#> HOVAL dy/dx        -3.3703965  0.15790853 -2.0552477
#> fDISCBDm2-m1 dy/dx  1.0175575 -2.26812418 -2.9855925
#> fDISCBDm3-m2 dy/dx  0.8525799 -0.01205903  0.9641081
#> fDISCBDm4-m3 dy/dx  1.0677351 -0.12020111  0.5012472
#> fEW1 dy/dx         -0.1625458  0.36558783  0.6788402
#> 
#> p-values:
#>                    Direct    Indirect Total    
#> INC dy/dx          0.0890254 0.790079 0.1663739
#> HOVAL dy/dx        0.0007506 0.874529 0.0398551
#> fDISCBDm2-m1 dy/dx 0.3088883 0.023322 0.0028303
#> fDISCBDm3-m2 dy/dx 0.3938923 0.990379 0.3349918
#> fDISCBDm4-m3 dy/dx 0.2856400 0.904324 0.6161972
#> fEW1 dy/dx         0.8708761 0.714673 0.4972391
#> 
# \dontrun{
mobj1 <- lagsarlm(f, columbus, listw, type="mixed", 
method="Matrix", control=list(fdHess=TRUE))
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(mobj1)
#> 
#> Call:lagsarlm(formula = f, data = columbus, listw = listw, type = "mixed", 
#>     method = "Matrix", control = list(fdHess = TRUE))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -3.485440 -0.303568  0.057651  0.412917  1.618395 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                    Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)       7.2505613  1.6930845  4.2825 1.848e-05
#> INC              -0.0511360  0.0317573 -1.6102 0.1073513
#> HOVAL            -0.0284100  0.0078616 -3.6138 0.0003018
#> fDISCBDm2-m1      0.2649702  0.4784131  0.5539 0.5796800
#> fDISCBDm3-m2      0.2832158  0.4791218  0.5911 0.5544438
#> fDISCBDm4-m3      0.5248924  0.5011489  1.0474 0.2949252
#> fEW1             -0.0391368  0.5390224 -0.0726 0.9421189
#> lag.INC          -0.0340158  0.0626273 -0.5431 0.5870287
#> lag.HOVAL        -0.0056745  0.0180359 -0.3146 0.7530494
#> lag.fDISCBDm2-m1 -1.7168592  0.7939419 -2.1624 0.0305835
#> lag.fDISCBDm3-m2  0.2492410  0.7084643  0.3518 0.7249848
#> lag.fDISCBDm4-m3 -0.3272604  1.1247038 -0.2910 0.7710707
#> lag.fEW1          0.2155133  0.6538676  0.3296 0.7417039
#> 
#> Rho: -0.43661, LR test value: 3.3537, p-value: 0.067054
#> Asymptotic standard error: 0.23242
#>     z-value: -1.8785, p-value: 0.060309
#> Wald statistic: 3.5289, p-value: 0.060309
#> 
#> Log likelihood: -58.21561 for mixed model
#> ML residual variance (sigma squared): 0.60545, (sigma: 0.7781)
#> Number of observations: 49 
#> Number of parameters estimated: 15 
#> AIC: 146.43, (AIC for lm: 147.78)
#> LM test for residual autocorrelation
#> test value: 5.1201, p-value: 0.02365
#> 
set.seed(1)
summary(impacts(mobj1, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136497 -1.442000095 -1.01063513
#> fDISCBDm3-m2 (F)  0.27181264  0.098821456  0.37063409
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848175  0.13756829
#> fEW1 (F)         -0.06025686  0.183029571  0.12277271
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                       Direct   Indirect      Total
#> INC dy/dx        0.034914580 0.05078645 0.04877969
#> HOVAL dy/dx      0.009168988 0.01503549 0.01169550
#> fDISCBDm2-m1 (F) 0.561090264 0.77511618 0.41334943
#> fDISCBDm3-m2 (F) 0.554689479 0.67343464 0.44513630
#> fDISCBDm4-m3 (F) 0.589205278 0.96930226 0.70220843
#> fEW1 (F)         0.648157823 0.72875362 0.23065107
#> 
#> Simulated z-values:
#>                      Direct   Indirect      Total
#> INC dy/dx        -1.4004887 -0.1571855 -1.1660665
#> HOVAL dy/dx      -3.2020383  0.3485614 -2.0622166
#> fDISCBDm2-m1 (F)  0.7698209 -1.8637180 -2.4498859
#> fDISCBDm3-m2 (F)  0.4854517  0.1178708  0.7832505
#> fDISCBDm4-m3 (F)  0.9955892 -0.4559156  0.2060447
#> fEW1 (F)         -0.1171795  0.2584545  0.4873111
#> 
#> Simulated p-values:
#>                  Direct    Indirect Total   
#> INC dy/dx        0.1613670 0.875099 0.243588
#> HOVAL dy/dx      0.0013646 0.727419 0.039187
#> fDISCBDm2-m1 (F) 0.4414062 0.062361 0.014290
#> fDISCBDm3-m2 (F) 0.6273560 0.906170 0.433480
#> fDISCBDm4-m3 (F) 0.3194498 0.648451 0.836756
#> fEW1 (F)         0.9067178 0.796056 0.626038
summary(impacts(mobj, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136498 -1.442000098 -1.01063512
#> fDISCBDm3-m2 (F)  0.27181264  0.098821437  0.37063408
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848169  0.13756830
#> fEW1 (F)         -0.06025686  0.183029567  0.12277271
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                       Direct   Indirect      Total
#> INC dy/dx        0.034361118 0.05820791 0.05682386
#> HOVAL dy/dx      0.009110538 0.01758732 0.01560487
#> fDISCBDm2-m1 (F) 0.600232069 0.93223139 0.66010753
#> fDISCBDm3-m2 (F) 0.566538983 0.69081171 0.48306608
#> fDISCBDm4-m3 (F) 0.581495739 1.02565992 0.78503361
#> fEW1 (F)         0.620713378 0.69372797 0.22641489
#> 
#> Simulated z-values:
#>                       Direct   Indirect      Total
#> INC dy/dx        -1.47846015 -0.1659242 -1.0639835
#> HOVAL dy/dx      -3.19397024  0.3457892 -1.4750061
#> fDISCBDm2-m1 (F)  0.72934305 -1.5269344 -1.4932129
#> fDISCBDm3-m2 (F)  0.51095323  0.1233050  0.7755781
#> fDISCBDm4-m3 (F)  1.01000469 -0.4304817  0.1857062
#> fEW1 (F)         -0.09272081  0.2678059  0.5663560
#> 
#> Simulated p-values:
#>                  Direct    Indirect Total  
#> INC dy/dx        0.1392847 0.86822  0.28734
#> HOVAL dy/dx      0.0014033 0.72950  0.14021
#> fDISCBDm2-m1 (F) 0.4657918 0.12678  0.13538
#> fDISCBDm3-m2 (F) 0.6093838 0.90187  0.43800
#> fDISCBDm4-m3 (F) 0.3124930 0.66685  0.85268
#> fEW1 (F)         0.9261254 0.78885  0.57115
mobj2 <- lagsarlm(f, columbus, listw, type="mixed", 
method="Matrix", control=list(fdHess=TRUE, optimHess=TRUE))
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(impacts(mobj2, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136497 -1.442000095 -1.01063513
#> fDISCBDm3-m2 (F)  0.27181264  0.098821456  0.37063409
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848175  0.13756829
#> fEW1 (F)         -0.06025686  0.183029571  0.12277271
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                       Direct   Indirect      Total
#> INC dy/dx        0.034560681 0.05202389 0.04934873
#> HOVAL dy/dx      0.008849542 0.01436242 0.01182192
#> fDISCBDm2-m1 (F) 0.583984901 0.79453308 0.40775202
#> fDISCBDm3-m2 (F) 0.534746156 0.68537387 0.45273784
#> fDISCBDm4-m3 (F) 0.585406472 1.01779477 0.74927843
#> fEW1 (F)         0.644626049 0.70756107 0.21952770
#> 
#> Simulated z-values:
#>                       Direct   Indirect      Total
#> INC dy/dx        -1.42127203 -0.1953328 -1.2012893
#> HOVAL dy/dx      -3.26432418  0.3833840 -1.9778059
#> fDISCBDm2-m1 (F)  0.76537261 -1.8647255 -2.5373755
#> fDISCBDm3-m2 (F)  0.46477572  0.1976818  0.8482237
#> fDISCBDm4-m3 (F)  0.95459098 -0.4237352  0.1702281
#> fEW1 (F)         -0.06488861  0.2291142  0.5479190
#> 
#> Simulated p-values:
#>                  Direct    Indirect Total   
#> INC dy/dx        0.1552377 0.84513  0.229639
#> HOVAL dy/dx      0.0010973 0.70144  0.047951
#> fDISCBDm2-m1 (F) 0.4440497 0.06222  0.011169
#> fDISCBDm3-m2 (F) 0.6420921 0.84329  0.396313
#> fDISCBDm4-m3 (F) 0.3397846 0.67176  0.864831
#> fEW1 (F)         0.9482627 0.81878  0.583748
mobj3 <- lagsarlm(f, columbus, listw, type="mixed", 
method="spam", control=list(fdHess=TRUE))
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
summary(impacts(mobj3, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
#> Impact measures (mixed, trace):
#>                       Direct     Indirect       Total
#> INC dy/dx        -0.05007460 -0.009198133 -0.05927273
#> HOVAL dy/dx      -0.02902104  0.005295429 -0.02372561
#> fDISCBDm2-m1 (F)  0.43136496 -1.442000094 -1.01063513
#> fDISCBDm3-m2 (F)  0.27181263  0.098821462  0.37063409
#> fDISCBDm4-m3 (F)  0.57541647 -0.437848177  0.13756829
#> fEW1 (F)         -0.06025686  0.183029573  0.12277271
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                       Direct   Indirect      Total
#> INC dy/dx        0.033745257 0.05307187 0.04915650
#> HOVAL dy/dx      0.008896117 0.01395402 0.01099169
#> fDISCBDm2-m1 (F) 0.561279468 0.78344519 0.41791425
#> fDISCBDm3-m2 (F) 0.544197996 0.66423748 0.43240507
#> fDISCBDm4-m3 (F) 0.593325811 0.97933534 0.71602876
#> fEW1 (F)         0.644079849 0.71604484 0.21992686
#> 
#> Simulated z-values:
#>                       Direct   Indirect      Total
#> INC dy/dx        -1.41272587 -0.1906571 -1.1756598
#> HOVAL dy/dx      -3.35155566  0.4551439 -2.1347701
#> fDISCBDm2-m1 (F)  0.77887522 -1.8508292 -2.4235991
#> fDISCBDm3-m2 (F)  0.52183626  0.1101818  0.8260059
#> fDISCBDm4-m3 (F)  0.97978462 -0.4581765  0.1852203
#> fEW1 (F)         -0.07870018  0.2331363  0.5285704
#> 
#> Simulated p-values:
#>                  Direct     Indirect Total   
#> INC dy/dx        0.15773634 0.848794 0.239731
#> HOVAL dy/dx      0.00080359 0.649006 0.032780
#> fDISCBDm2-m1 (F) 0.43605322 0.064194 0.015368
#> fDISCBDm3-m2 (F) 0.60178434 0.912265 0.408801
#> fDISCBDm4-m3 (F) 0.32719245 0.646826 0.853056
#> fEW1 (F)         0.93727110 0.815656 0.597104
# }
# \dontrun{
data(boston, package="spData")
Wb <- as(spdep::nb2listw(boston.soi), "CsparseMatrix")
trMatb <- trW(Wb, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, spdep::nb2listw(boston.soi), Durbin=TRUE, method="Matrix", 
control=list(fdHess=TRUE), trs=trMatb)
#> Warning: use of spatially lagged factor (categorical variable)
#> CHAS
#> is not well-understood
summary(gp2mMi)
#> 
#> Call:lagsarlm(formula = log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
#>     I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + 
#>     log(LSTAT), data = boston.c, listw = spdep::nb2listw(boston.soi), 
#>     Durbin = TRUE, method = "Matrix", trs = trMatb, control = list(fdHess = TRUE))
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.6316833 -0.0629790 -0.0090776  0.0682421  0.6991072 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                   Estimate  Std. Error  z value  Pr(>|z|)
#> (Intercept)     1.89816225  0.24400294   7.7793 7.327e-15
#> CRIM           -0.00571021  0.00093857  -6.0839 1.173e-09
#> ZN              0.00069091  0.00051874   1.3319 0.1828931
#> INDUS          -0.00111343  0.00307380  -0.3622 0.7171778
#> CHAS1          -0.04163225  0.02738840  -1.5201 0.1284937
#> I(NOX^2)       -0.01034950  0.19358633  -0.0535 0.9573639
#> I(RM^2)         0.00794979  0.00102109   7.7856 6.883e-15
#> AGE            -0.00128789  0.00048929  -2.6322 0.0084838
#> log(DIS)       -0.12404108  0.09510145  -1.3043 0.1921304
#> log(RAD)        0.05863502  0.02257529   2.5973 0.0093957
#> TAX            -0.00049084  0.00012146  -4.0412 5.317e-05
#> PTRATIO        -0.01319853  0.00595331  -2.2170 0.0266227
#> B               0.00056383  0.00011084   5.0867 3.643e-07
#> log(LSTAT)     -0.24724454  0.02265150 -10.9152 < 2.2e-16
#> lag.CRIM       -0.00464215  0.00173900  -2.6694 0.0075978
#> lag.ZN         -0.00037937  0.00070703  -0.5366 0.5915659
#> lag.INDUS       0.00025064  0.00385911   0.0649 0.9482165
#> lag.CHAS1       0.12518252  0.04083949   3.0652 0.0021750
#> lag.I(NOX^2)   -0.38640403  0.22253432  -1.7364 0.0824968
#> lag.I(RM^2)    -0.00451252  0.00148919  -3.0302 0.0024440
#> lag.AGE         0.00149678  0.00068470   2.1860 0.0288128
#> lag.log(DIS)   -0.00453785  0.10056962  -0.0451 0.9640104
#> lag.log(RAD)   -0.00940702  0.03127788  -0.3008 0.7636002
#> lag.TAX         0.00041083  0.00017859   2.3004 0.0214237
#> lag.PTRATIO     0.00060355  0.00789995   0.0764 0.9391011
#> lag.B          -0.00050781  0.00014107  -3.5996 0.0003187
#> lag.log(LSTAT)  0.09846780  0.03574187   2.7550 0.0058697
#> 
#> Rho: 0.59578, LR test value: 181.68, p-value: < 2.22e-16
#> Asymptotic standard error: 0.037474
#>     z-value: 15.898, p-value: < 2.22e-16
#> Wald statistic: 252.76, p-value: < 2.22e-16
#> 
#> Log likelihood: 300.6131 for mixed model
#> ML residual variance (sigma squared): 0.016011, (sigma: 0.12654)
#> Number of observations: 506 
#> Number of parameters estimated: 29 
#> AIC: -543.23, (AIC for lm: -363.55)
#> LM test for residual autocorrelation
#> test value: 29.772, p-value: 4.8604e-08
#> 
summary(impacts(gp2mMi, tr=trMatb, R=1000), zstats=TRUE, short=TRUE)
#> Impact measures (mixed, trace):
#>                         Direct      Indirect         Total
#> CRIM dy/dx       -0.0074555753 -0.0181548470 -0.0256104223
#> ZN dy/dx          0.0006979073  0.0000727985  0.0007707058
#> INDUS dy/dx      -0.0012029822 -0.0009314672 -0.0021344494
#> CHAS1 (F)        -0.0198526431  0.2265453878  0.2066927447
#> I(NOX^2) dy/dx   -0.0955268251 -0.8859909283 -0.9815177534
#> I(RM^2) dy/dx     0.0079983430  0.0005050175  0.0085033605
#> AGE dy/dx        -0.0011296134  0.0016463650  0.0005167515
#> log(DIS) dy/dx   -0.1410601685 -0.1770277708 -0.3180879393
#> log(RAD) dy/dx    0.0641735546  0.0576102594  0.1217838140
#> TAX dy/dx        -0.0004651543  0.0002672119 -0.0001979424
#> PTRATIO dy/dx    -0.0147737151 -0.0163846522 -0.0311583673
#> B dy/dx           0.0005265343 -0.0003879424  0.0001385920
#> log(LSTAT) dy/dx -0.2578403220 -0.1102143882 -0.3680547102
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                        Direct     Indirect        Total
#> CRIM dy/dx       0.0010187475 0.0034731614 0.0039572691
#> ZN dy/dx         0.0004953559 0.0012372378 0.0012649951
#> INDUS dy/dx      0.0028029574 0.0057197966 0.0051920475
#> CHAS1 (F)        0.0277058122 0.0826448563 0.0915400508
#> I(NOX^2) dy/dx   0.1710489720 0.2781318249 0.2416975710
#> I(RM^2) dy/dx    0.0010767044 0.0030410036 0.0034987019
#> AGE dy/dx        0.0004786794 0.0012375171 0.0013167320
#> log(DIS) dy/dx   0.0852410189 0.1084881326 0.0760790272
#> log(RAD) dy/dx   0.0204736160 0.0506486494 0.0490276895
#> TAX dy/dx        0.0001195085 0.0003163806 0.0003249098
#> PTRATIO dy/dx    0.0058207887 0.0123533893 0.0119644745
#> B dy/dx          0.0001065122 0.0002368501 0.0002380004
#> log(LSTAT) dy/dx 0.0218703588 0.0580567220 0.0619685540
#> 
#> Simulated z-values:
#>                       Direct    Indirect      Total
#> CRIM dy/dx        -7.3203741 -5.22085306 -6.4667016
#> ZN dy/dx           1.4291232  0.05285644  0.6113229
#> INDUS dy/dx       -0.4039075 -0.14710985 -0.3801148
#> CHAS1 (F)         -0.7277550  2.75235020  2.2646322
#> I(NOX^2) dy/dx    -0.5618575 -3.17743923 -4.0540421
#> I(RM^2) dy/dx      7.4516491  0.17482235  2.4451523
#> AGE dy/dx         -2.3631765  1.33754263  0.3979763
#> log(DIS) dy/dx    -1.6500370 -1.64668260 -4.1969038
#> log(RAD) dy/dx     3.0978694  1.16856102  2.5008445
#> TAX dy/dx         -3.8848395  0.78890626 -0.6607271
#> PTRATIO dy/dx     -2.4869923 -1.33375639 -2.5870479
#> B dy/dx            4.9617418 -1.59801003  0.6302394
#> log(LSTAT) dy/dx -11.8101114 -1.87918924 -5.9286673
#> 
#> Simulated p-values:
#>                  Direct     Indirect  Total     
#> CRIM dy/dx       2.4736e-13 1.781e-07 1.0017e-10
#> ZN dy/dx         0.1529688  0.9578463 0.5409858 
#> INDUS dy/dx      0.6862807  0.8830453 0.7038602 
#> CHAS1 (F)        0.4667636  0.0059169 0.0235353 
#> I(NOX^2) dy/dx   0.5742131  0.0014858 5.0340e-05
#> I(RM^2) dy/dx    9.2149e-14 0.8612192 0.0144791 
#> AGE dy/dx        0.0181190  0.1810456 0.6906476 
#> log(DIS) dy/dx   0.0989354  0.0996233 2.7059e-05
#> log(RAD) dy/dx   0.0019492  0.2425805 0.0123898 
#> TAX dy/dx        0.0001024  0.4301668 0.5087873 
#> PTRATIO dy/dx    0.0128828  0.1822837 0.0096802 
#> B dy/dx          6.9864e-07 0.1100407 0.5285379 
#> log(LSTAT) dy/dx < 2.22e-16 0.0602187 3.0540e-09
#data(house, package="spData")
#lw <- spdep::nb2listw(LO_nb)
#form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
#   rooms + log(TLA) + beds + syear)
#lobj <- lagsarlm(form, house, lw, method="Matrix",
# control=list(fdHess=TRUE), trs=trMat)
#summary(lobj)
#loobj <- impacts(lobj, tr=trMat, R=1000)
#summary(loobj, zstats=TRUE, short=TRUE)
#lobj1 <- stsls(form, house, lw)
#loobj1 <- impacts(lobj1, tr=trMat, R=1000)
#summary(loobj1, zstats=TRUE, short=TRUE)
#mobj <- lagsarlm(form, house, lw, type="mixed",
# method="Matrix", control=list(fdHess=TRUE), trs=trMat)
#summary(mobj)
#moobj <- impacts(mobj, tr=trMat, R=1000)
#summary(moobj, zstats=TRUE, short=TRUE)
# }
```
