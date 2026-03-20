# Bayesian MCMC spatial simultaneous autoregressive model estimation

The `spBreg_lag` function is an early-release version of the Matlab
Spatial Econometrics Toolbox function `sar_g.m`, using drawing by
inversion, and not accommodating heteroskedastic disturbances.

## Usage

``` r
spBreg_lag(formula, data = list(), listw, na.action, Durbin, type,
    zero.policy=NULL, control=list())
spBreg_sac(formula, data = list(), listw, listw2=NULL, na.action, 
    Durbin, type, zero.policy=NULL, control=list())
spBreg_err(formula, data = list(), listw, na.action, Durbin, etype,
    zero.policy=NULL, control=list())
# S3 method for class 'MCMC_sar_G'
impacts(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
# S3 method for class 'MCMC_sem_G'
impacts(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
# S3 method for class 'MCMC_sac_G'
impacts(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
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

  a `listw` object created for example by `nb2listw`

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the weights list will be subsetted to remove NAs in the
  data. It may be necessary to set zero.policy to TRUE because this
  subsetting may create no-neighbour observations. Note that only
  weights lists created without using the glist argument to `nb2listw`
  may be subsetted.

- Durbin:

  default FALSE (spatial lag model); if TRUE, full spatial Durbin model;
  if a formula object, the subset of explanatory variables to lag. From
  version 1.3-7, the presence of factors (categorical variables) in the
  Durbin term will give a warning, as it is as yet unknown how spatial
  lags of categorical variables should be interpreted.

- type, etype:

  (use the ‘Durbin=’ argument - retained for backwards compatibility
  only) default "lag", may be set to "mixed"; when "mixed", the lagged
  intercept is dropped for spatial weights style "W", that is
  row-standardised weights, but otherwise included; “Durbin” may be used
  instead of “mixed”

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE (default) assign NA

- control:

  list of extra control arguments - see section below

- obj:

  A spatial regression object

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

- Q:

  default NULL, else an integer number of cumulative power series
  impacts to calculate if `tr` is given

## Control arguments

- tol.opt::

  the desired accuracy of the optimization - passed to
  [`optimize()`](https://rdrr.io/r/stats/optimize.html) (default=square
  root of double precision machine tolerance, a larger root may be used
  needed, see help(boston) for an example)

- fdHess::

  default NULL, then set to (method != "eigen") internally; use `fdHess`
  to compute an approximate Hessian using finite differences when using
  sparse matrix methods; used to make a coefficient covariance matrix
  when the number of observations is large; may be turned off to save
  resources if need be

- optimHess::

  default FALSE, use `fdHess` from nlme, if TRUE, use `optim` to
  calculate Hessian at optimum

- optimHessMethod::

  default “optimHess”, may be “nlm” or one of the `optim` methods

- compiled_sse::

  default FALSE; logical value used in the log likelihood function to
  choose compiled code for computing SSE

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

- OrdVsign:

  default 1; used to set the sign of the final component to negative if
  -1 (alpha times ((sigma squared) squared) in Ord (1975) equation B.1).

## Extra Bayesian control arguments

- ldet_method:

  default “SE_classic”; equivalent to the `method` argument in
  `lagsarlm`

- interval:

  default `c(-1, 1)`; used unmodified or set internally by
  `jacobianSetup`

- ndraw:

  default `2500L`; integer total number of draws

- nomit:

  default `500L`; integer total number of omitted burn-in draws

- thin:

  default `1L`; integer thinning proportion

- verbose:

  default `FALSE`; inverse of `quiet` argument in `lagsarlm`

- detval:

  default `NULL`; not yet in use, precomputed matrix of log determinants

- prior:

  a list with the following components:

  rhoMH, lambdaMH

  :   default FALSE; use Metropolis or griddy Gibbs

  Tbeta

  :   default `NULL`; values of the betas variance-covariance matrix,
      set to `diag(k)*1e+12` if `NULL`

  c_beta

  :   default `NULL`; values of the betas set to 0 if `NULL`

  rho

  :   default `0.5`; value of the autoregressive coefficient

  sige

  :   default `1`; value of the residual variance

  nu

  :   default `0`; informative Gamma(nu,d0) prior on sige

  d0

  :   default `0`; informative Gamma(nu,d0) prior on sige

  a1

  :   default `1.01`; parameter for beta(a1,a2) prior on rho

  a2

  :   default `1.01`; parameter for beta(a1,a2) prior on rho

  cc

  :   default `0.2`; initial tuning parameter for M-H sampling

  gG_sige

  :   default TRUE; include sige in lambda griddy Gibbs update

  cc1

  :   default `0.2`; initial tuning parameter for M-H sampling

  cc2

  :   default `0.2`; initial tuning parameter for M-H sampling

## Details

From version 1.4.1, functions for models including spatially lagged
independent variables warn on fitting if any of the right-hand side
variables are factors. This is because the interpretation of
coefficients that are not slopes is unclear when the variable is not
interpretable on an unbounded line, such as factors. Factor variable
names are shown with the suffix “(F)”, others “dy/dx” in output from
impact methods. A discussion can be found at
<https://github.com/rsbivand/eqc25_talk>.

## References

LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC
Press, Boca Raton.

## Author

Roger Bivand <Roger.Bivand@nhh.no>, with thanks to Abhirup Mallik and
Virgilio Gómez-Rubio for initial coding GSoC 2011

## Examples

``` r
#require("spdep", quietly=TRUE)
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb, style="W")
require("coda", quietly=TRUE)
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD Naive SE Time-series SE
#> (Intercept)  59.6679  6.9955 0.156423       0.156423
#> INC          -0.9460  0.4018 0.008985       0.008684
#> HOVAL        -0.2987  0.1001 0.002239       0.002239
#> lambda        0.5664  0.1536 0.003435       0.003435
#> sige        115.0536 26.1966 0.585775       0.636948
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%      50%      75%    97.5%
#> (Intercept) 45.6943 55.3584  60.0629  64.5402  72.2670
#> INC         -1.7514 -1.2092  -0.9559  -0.6698  -0.1607
#> HOVAL       -0.4976 -0.3668  -0.2982  -0.2320  -0.1008
#> lambda       0.2335  0.4717   0.5773   0.6758   0.8309
#> sige        74.7007 96.1409 111.5963 129.1220 177.0618
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 3        1052  937          1.120     
#>  INC         2        969   937          1.030     
#>  HOVAL       2        969   937          1.030     
#>  lambda      3        1010  937          1.080     
#>  sige        2        930   937          0.993     
#> 
# \dontrun{
ev <- eigenw(lw)
W <- as(lw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  59.4996  7.82235 0.174913       0.215238
#> INC          -0.9212  0.39143 0.008753       0.013191
#> HOVAL        -0.2982  0.09838 0.002200       0.002368
#> lambda        0.6005  0.14838 0.003318       0.007353
#> sige        116.6722 27.41993 0.613128       0.664469
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%      50%      75%    97.5%
#> (Intercept) 44.4995 55.3154  59.5282  64.3852  73.5281
#> INC         -1.6697 -1.1806  -0.9244  -0.6632  -0.1450
#> HOVAL       -0.4871 -0.3637  -0.2979  -0.2317  -0.1033
#> lambda       0.2925  0.5008   0.6080   0.7029   0.8732
#> sige        75.2439 96.7819 113.3399 131.7463 181.0479
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 4        1295  937          1.380     
#>  INC         2        969   937          1.030     
#>  HOVAL       3        1010  937          1.080     
#>  lambda      12       3372  937          3.600     
#>  sige        2        930   937          0.993     
#> 
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE)
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                  Mean      SD Naive SE Time-series SE
#> (Intercept)  72.65867 11.7955 0.263755       0.263755
#> INC          -1.02285  0.3806 0.008510       0.008510
#> HOVAL        -0.27805  0.1055 0.002358       0.002358
#> lag.INC      -1.07102  0.7277 0.016272       0.016272
#> lag.HOVAL     0.09426  0.2377 0.005315       0.005150
#> lambda        0.48878  0.1732 0.003873       0.003873
#> sige        114.74074 26.8443 0.600257       0.653655
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%      25%      50%      75%    97.5%
#> (Intercept) 48.2209 65.56843  73.0461  80.0125  95.7722
#> INC         -1.7289 -1.28403  -1.0319  -0.7678  -0.2625
#> HOVAL       -0.4880 -0.34810  -0.2771  -0.2068  -0.0765
#> lag.INC     -2.4400 -1.57056  -1.1057  -0.5923   0.3666
#> lag.HOVAL   -0.3749 -0.06271   0.1014   0.2545   0.5513
#> lambda       0.1165  0.37669   0.4977   0.6118   0.7969
#> sige        73.4577 95.82398 110.8881 129.0908 177.2002
#> 
print(summary(impacts(COL.err.Bayes)))
#> Impact measures (SDEM, MCMC, n):
#>                 Direct    Indirect      Total
#> INC dy/dx   -1.0228452 -1.07102436 -2.0938696
#> HOVAL dy/dx -0.2780522  0.09426141 -0.1837908
#> ========================================================
#> Standard errors:
#>                 Direct  Indirect     Total
#> INC dy/dx   0.36063747 0.6895763 0.8293833
#> HOVAL dy/dx 0.09993534 0.2252272 0.2714330
#> ========================================================
#> Z-values:
#>                Direct   Indirect      Total
#> INC dy/dx   -2.836214 -1.5531629 -2.5246103
#> HOVAL dy/dx -2.782321  0.4185171 -0.6771128
#> 
#> p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0045652 0.12038  0.011583
#> HOVAL dy/dx 0.0053972 0.67557  0.498334
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 2        930   937          0.993     
#>  INC         3        1010  937          1.080     
#>  HOVAL       2        969   937          1.030     
#>  lag.INC     2        969   937          1.030     
#>  lag.HOVAL   2        930   937          0.993     
#>  lambda      2        930   937          0.993     
#>  sige        2        930   937          0.993     
#> 
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                  Mean      SD Naive SE Time-series SE
#> (Intercept)  71.85960 14.3318 0.320468       0.320468
#> INC          -0.99458  0.3896 0.008711       0.009242
#> HOVAL        -0.28057  0.1071 0.002395       0.002395
#> lag.INC      -0.99171  0.8076 0.018058       0.018989
#> lag.HOVAL     0.07993  0.2534 0.005667       0.006244
#> lambda        0.56313  0.1956 0.004374       0.011987
#> sige        118.58119 28.8191 0.644415       0.904218
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%      25%       50%      75%     97.5%
#> (Intercept) 42.5325 64.10911  72.28760  79.9509  98.37754
#> INC         -1.7588 -1.25887  -0.99557  -0.7325  -0.20956
#> HOVAL       -0.4815 -0.35123  -0.28296  -0.2080  -0.06975
#> lag.INC     -2.4108 -1.54399  -1.03869  -0.5151   0.83601
#> lag.HOVAL   -0.4475 -0.09169   0.08381   0.2529   0.56757
#> lambda       0.1178  0.44211   0.57614   0.7027   0.90737
#> sige        75.6445 97.86837 113.74179 133.5623 187.50180
#> 
print(summary(impacts(COL.err.Bayes)))
#> Impact measures (SDEM, MCMC, n):
#>                 Direct    Indirect      Total
#> INC dy/dx   -0.9945809 -0.99171062 -1.9862916
#> HOVAL dy/dx -0.2805719  0.07993348 -0.2006384
#> ========================================================
#> Standard errors:
#>                Direct  Indirect     Total
#> INC dy/dx   0.3691657 0.7652838 0.9322001
#> HOVAL dy/dx 0.1014867 0.2401648 0.2943969
#> ========================================================
#> Z-values:
#>                Direct   Indirect      Total
#> INC dy/dx   -2.694131 -1.2958731 -2.1307566
#> HOVAL dy/dx -2.764617  0.3328276 -0.6815234
#> 
#> p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0070572 0.19502  0.033109
#> HOVAL dy/dx 0.0056990 0.73926  0.495540
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 3        1052  937          1.120     
#>  INC         3        1010  937          1.080     
#>  HOVAL       2        930   937          0.993     
#>  lag.INC     2        930   937          0.993     
#>  lag.HOVAL   2        892   937          0.952     
#>  lambda      24       6635  937          7.080     
#>  sige        2        930   937          0.993     
#> 
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC)
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  74.5363 10.99444 0.245843       0.245843
#> INC          -1.0046  0.36864 0.008243       0.008243
#> HOVAL        -0.2897  0.09787 0.002188       0.002188
#> lag.INC      -0.9385  0.58827 0.013154       0.012742
#> lambda        0.4879  0.17196 0.003845       0.003672
#> sige        111.8908 25.45327 0.569152       0.566145
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%      50%      75%    97.5%
#> (Intercept) 52.5957 67.7552  74.8200  81.6559  95.1525
#> INC         -1.7327 -1.2421  -1.0142  -0.7527  -0.2977
#> HOVAL       -0.4813 -0.3538  -0.2898  -0.2244  -0.1006
#> lag.INC     -2.1220 -1.3114  -0.9574  -0.5528   0.2343
#> lambda       0.1355  0.3737   0.4987   0.6088   0.8010
#> sige        73.2004 93.5798 108.3480 125.6844 174.2070
#> 
print(summary(impacts(COL.err.Bayes)))
#> Impact measures (SDEM, MCMC, n):
#>                 Direct   Indirect      Total
#> INC dy/dx   -1.0046355 -0.9385296 -1.9431651
#> HOVAL dy/dx -0.2897097         NA -0.2897097
#> ========================================================
#> Standard errors:
#>                 Direct  Indirect      Total
#> INC dy/dx   0.35327149 0.5637475 0.68482488
#> HOVAL dy/dx 0.09379138        NA 0.09379138
#> ========================================================
#> Z-values:
#>                Direct  Indirect     Total
#> INC dy/dx   -2.843806 -1.664805 -2.837463
#> HOVAL dy/dx -3.088873        NA -3.088873
#> 
#> p-values:
#>             Direct    Indirect Total    
#> INC dy/dx   0.0044578 0.095952 0.0045474
#> HOVAL dy/dx 0.0020092 NA       0.0020092
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 2        969   937          1.030     
#>  INC         2        892   937          0.952     
#>  HOVAL       2        892   937          0.952     
#>  lag.INC     2        930   937          0.993     
#>  lambda      2        892   937          0.952     
#>  sige        2        969   937          1.030     
#> 
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  73.0959 14.04252 0.314000       0.343190
#> INC          -0.9880  0.38386 0.008583       0.008583
#> HOVAL        -0.2876  0.09897 0.002213       0.002213
#> lag.INC      -0.8747  0.68850 0.015395       0.016417
#> lambda        0.5826  0.16471 0.003683       0.009205
#> sige        115.7170 27.07951 0.605516       0.767425
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%      50%      75%     97.5%
#> (Intercept) 44.3417 65.4362  73.5329  81.4403  99.71035
#> INC         -1.7299 -1.2493  -0.9925  -0.7241  -0.24351
#> HOVAL       -0.4812 -0.3536  -0.2873  -0.2222  -0.09316
#> lag.INC     -2.1791 -1.3231  -0.8918  -0.4329   0.51982
#> lambda       0.2286  0.4699   0.5922   0.7004   0.87677
#> sige        74.7196 96.6325 111.5361 130.5380 178.38293
#> 
print(summary(impacts(COL.err.Bayes)))
#> Impact measures (SDEM, MCMC, n):
#>                 Direct   Indirect      Total
#> INC dy/dx   -0.9880208 -0.8747392 -1.8627600
#> HOVAL dy/dx -0.2875558         NA -0.2875558
#> ========================================================
#> Standard errors:
#>                 Direct  Indirect      Total
#> INC dy/dx   0.36785864 0.6597952 0.82204787
#> HOVAL dy/dx 0.09484316        NA 0.09484316
#> ========================================================
#> Z-values:
#>                Direct  Indirect     Total
#> INC dy/dx   -2.685871 -1.325774 -2.265999
#> HOVAL dy/dx -3.031909        NA -3.031909
#> 
#> p-values:
#>             Direct    Indirect Total    
#> INC dy/dx   0.0072341 0.18491  0.0234514
#> HOVAL dy/dx 0.0024301 NA       0.0024301
#> 
print(raftery.diag(COL.err.Bayes, r=0.01))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.01
#> Probability (s) = 0.95 
#>                                                    
#>              Burn-in  Total Lower bound  Dependence
#>              (M)      (N)   (Nmin)       factor (I)
#>  (Intercept) 4        1192  937          1.270     
#>  INC         2        930   937          0.993     
#>  HOVAL       2        892   937          0.952     
#>  lag.INC     2        969   937          1.030     
#>  lambda      11       3046  937          3.250     
#>  sige        2        930   937          0.993     
#> 
set.seed(1)
COL.sacW.B0 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=FALSE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B0))
#> 
#> Iterations = 501:1500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 1000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD Naive SE Time-series SE
#> (Intercept)  49.3764  9.2344 0.292017       0.845878
#> INC          -1.0228  0.3490 0.011038       0.011716
#> HOVAL        -0.2765  0.1004 0.003176       0.003774
#> rho           0.3117  0.1963 0.006208       0.020392
#> lambda        0.1927  0.2713 0.008579       0.027162
#> sige        105.0462 23.9681 0.757936       0.757936
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%      25%      50%      75%    97.5%
#> (Intercept) 31.7652 43.16427  48.9447  55.3343  68.1554
#> INC         -1.7141 -1.26321  -1.0246  -0.7873  -0.3626
#> HOVAL       -0.4703 -0.34727  -0.2756  -0.2074  -0.0829
#> rho         -0.1548  0.21285   0.3406   0.4438   0.6176
#> lambda      -0.3788  0.01455   0.2171   0.3794   0.6988
#> sige        68.5025 88.13909 101.6888 118.2545 160.2792
#> 
print(summary(impacts(COL.sacW.B0, tr=trMatc), zstats=TRUE, short=TRUE))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': error in evaluating the argument 'object' in selecting a method for function 'summary': !is.null(have_factor_preds) is not TRUE
set.seed(1)
COL.sacW.B1 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B1))
#> 
#> Iterations = 501:1500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 1000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  60.1933 24.16199 0.764069       3.169235
#> INC          -0.9933  0.35907 0.011355       0.015218
#> HOVAL        -0.2853  0.09553 0.003021       0.003325
#> lag.INC      -0.8653  0.78593 0.024853       0.057898
#> lag.HOVAL     0.1677  0.22105 0.006990       0.013123
#> rho           0.1808  0.32007 0.010121       0.044327
#> lambda        0.2106  0.32071 0.010142       0.041657
#> sige        100.9463 24.32264 0.769149       0.946662
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%       25%     50%      75%    97.5%
#> (Intercept) 22.1665 42.233072 56.3130  76.2570 113.7785
#> INC         -1.6894 -1.233410 -1.0113  -0.7459  -0.2608
#> HOVAL       -0.4686 -0.352240 -0.2824  -0.2201  -0.1038
#> lag.INC     -2.5043 -1.361233 -0.8048  -0.3001   0.5206
#> lag.HOVAL   -0.2868  0.029384  0.1760   0.3192   0.5787
#> rho         -0.5076 -0.021227  0.2386   0.4320   0.6655
#> lambda      -0.5106 -0.007594  0.2129   0.4482   0.7389
#> sige        64.6867 84.395156 96.9599 114.4101 159.7252
#> 
print(summary(impacts(COL.sacW.B1, tr=trMatc), zstats=TRUE, short=TRUE))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': error in evaluating the argument 'object' in selecting a method for function 'summary': !is.null(have_factor_preds) is not TRUE
set.seed(1)
COL.lag.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=lw)
print(summary(COL.lag.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  45.8263  8.05132 0.180033       0.180033
#> INC          -1.0494  0.34943 0.007813       0.007559
#> HOVAL        -0.2640  0.09535 0.002132       0.002132
#> rho           0.4146  0.12434 0.002780       0.002780
#> sige        107.8827 23.75095 0.531087       0.561152
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%      50%      75%     97.5%
#> (Intercept) 30.2538 40.2712  45.9327  51.0968  61.45198
#> INC         -1.7554 -1.2823  -1.0501  -0.8146  -0.37317
#> HOVAL       -0.4546 -0.3259  -0.2642  -0.2007  -0.07478
#> rho          0.1645  0.3327   0.4157   0.5008   0.64982
#> sige        70.5778 90.6316 104.8556 120.4706 164.75421
#> 
print(summary(impacts(COL.lag.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
#> Impact measures (lag, trace):
#>                Direct   Indirect      Total
#> INC dy/dx   -1.099838 -0.6926466 -1.7924845
#> HOVAL dy/dx -0.276692 -0.1742527 -0.4509447
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct  Indirect     Total
#> INC dy/dx   0.3746613 0.5380947 0.8255332
#> HOVAL dy/dx 0.1008362 0.1265110 0.2037510
#> 
#> Simulated z-values:
#>                Direct  Indirect     Total
#> INC dy/dx   -2.963118 -1.451877 -2.291141
#> HOVAL dy/dx -2.766217 -1.525740 -2.316345
#> 
#> Simulated p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0030454 0.14654  0.021955
#> HOVAL dy/dx 0.0056711 0.12707  0.020539
print(summary(impacts(COL.lag.Bayes, evalues=ev), short=TRUE, zstats=TRUE))
#> Impact measures (lag, evalues):
#>                Direct   Indirect      Total
#> INC dy/dx   -1.099838 -0.6926466 -1.7924845
#> HOVAL dy/dx -0.276692 -0.1742527 -0.4509447
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct  Indirect     Total
#> INC dy/dx   0.3746642 0.5385717 0.8258978
#> HOVAL dy/dx 0.1008363 0.1265680 0.2037897
#> 
#> Simulated z-values:
#>                Direct  Indirect     Total
#> INC dy/dx   -2.963099 -1.450670 -2.290182
#> HOVAL dy/dx -2.766215 -1.525111 -2.315941
#> 
#> Simulated p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0030456 0.14687  0.022011
#> HOVAL dy/dx 0.0056711 0.12723  0.020561
set.seed(1)
COL.D0.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=lw, Durbin=TRUE)
print(summary(COL.D0.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  45.3083 13.75028 0.307466       0.307466
#> INC          -0.9310  0.36730 0.008213       0.008213
#> HOVAL        -0.2930  0.09694 0.002168       0.002168
#> lag.INC      -0.5808  0.62019 0.013868       0.013868
#> lag.HOVAL     0.2404  0.19496 0.004359       0.004359
#> rho           0.3932  0.16296 0.003644       0.003644
#> sige        109.7428 24.97328 0.558419       0.607251
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%     25%      50%      75%    97.5%
#> (Intercept) 19.29927 35.9473  44.9971  54.5571  72.7712
#> INC         -1.64191 -1.1710  -0.9344  -0.6855  -0.1947
#> HOVAL       -0.48457 -0.3590  -0.2923  -0.2268  -0.1030
#> lag.INC     -1.77614 -1.0145  -0.5787  -0.1582   0.6293
#> lag.HOVAL   -0.13457  0.1025   0.2446   0.3694   0.6273
#> rho          0.06053  0.2826   0.3977   0.5078   0.6949
#> sige        71.19421 91.4872 106.1690 123.7826 170.9503
#> 
print(summary(impacts(COL.D0.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
#> Impact measures (mixed, trace):
#>                Direct   Indirect       Total
#> INC dy/dx   -1.033240 -1.4584485 -2.49168807
#> HOVAL dy/dx -0.279485  0.1928629 -0.08662218
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                Direct  Indirect     Total
#> INC dy/dx   0.3830078 1.2651790 1.4115054
#> HOVAL dy/dx 0.1009060 0.3398633 0.3770575
#> 
#> Simulated z-values:
#>                Direct   Indirect      Total
#> INC dy/dx   -2.747755 -1.2977011 -1.9087676
#> HOVAL dy/dx -2.786365  0.5517863 -0.2483154
#> 
#> Simulated p-values:
#>             Direct    Indirect Total   
#> INC dy/dx   0.0060005 0.19439  0.056292
#> HOVAL dy/dx 0.0053303 0.58109  0.803890
set.seed(1)
COL.D1.Bayes <- spBreg_lag(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 listw=lw, Durbin= ~ INC)
print(summary(COL.D1.Bayes))
#> 
#> Iterations = 501:2500
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 2000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD Naive SE Time-series SE
#> (Intercept)  56.4264 13.44686 0.300681       0.300681
#> DISCBD       -4.7754  2.08410 0.046602       0.046602
#> INC          -0.9480  0.35344 0.007903       0.007903
#> HOVAL        -0.1810  0.09891 0.002212       0.002212
#> lag.INC       0.4314  0.64207 0.014357       0.014357
#> rho           0.1873  0.18174 0.004064       0.004064
#> sige        103.8701 23.29384 0.520866       0.566802
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%      25%      50%      75%     97.5%
#> (Intercept) 30.5876 47.91538  56.5688  65.0907  84.59242
#> DISCBD      -8.7372 -6.20095  -4.7088  -3.3578  -0.69914
#> INC         -1.6262 -1.18846  -0.9496  -0.7131  -0.24344
#> HOVAL       -0.3701 -0.24813  -0.1822  -0.1145   0.01959
#> lag.INC     -0.8156  0.01182   0.4262   0.8481   1.72379
#> rho         -0.2116  0.07254   0.1836   0.3067   0.53579
#> sige        67.8946 87.05494 100.9815 116.4528 161.25826
#> 
print(summary(impacts(COL.D1.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
#> Impact measures (mixed, trace):
#>                  Direct    Indirect      Total
#> DISCBD dy/dx -4.8154723 -1.06081543 -5.8762877
#> INC dy/dx    -0.9365992  0.30093299 -0.6356662
#> HOVAL dy/dx  -0.1825344 -0.04021108 -0.2227455
#> ========================================================
#> Simulation results ( variance matrix):
#> ========================================================
#> Simulated standard errors
#>                 Direct  Indirect     Total
#> DISCBD dy/dx 2.1266139 1.7089367 3.2153297
#> INC dy/dx    0.3577465 0.8401827 0.9226330
#> HOVAL dy/dx  0.1009231 0.0708645 0.1460887
#> 
#> Simulated z-values:
#>                 Direct   Indirect      Total
#> DISCBD dy/dx -2.287425 -0.7858770 -1.9305901
#> INC dy/dx    -2.636247  0.3294488 -0.7221841
#> HOVAL dy/dx  -1.826954 -0.7109182 -1.6069734
#> 
#> Simulated p-values:
#>              Direct    Indirect Total   
#> DISCBD dy/dx 0.0221710 0.43194  0.053534
#> INC dy/dx    0.0083829 0.74182  0.470181
#> HOVAL dy/dx  0.0677067 0.47713  0.108060
#data(elect80, package="spData")
#lw <- spdep::nb2listw(e80_queen, zero.policy=TRUE)
#el_ml <- lagsarlm(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
# + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE, method="LU")
#print(summary(el_ml))
#set.seed(1)
#el_B <- spBreg_lag(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
# + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE)
#print(summary(el_B))
#print(el_ml$timings)
#print(attr(el_B, "timings"))
# }
```
