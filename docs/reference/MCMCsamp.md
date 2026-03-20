# MCMC sample from fitted spatial regression

The `MCMCsamp` method uses
[`rwmetrop`](https://rdrr.io/pkg/LearnBayes/man/rwmetrop.html), a random
walk Metropolis algorithm, from LearnBayes to make MCMC samples from
fitted maximum likelihood spatial regression models.

## Usage

``` r
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...)
# S3 method for class 'Spautolm'
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...,
 burnin = 0L, scale=1, listw, control = list())
# S3 method for class 'Sarlm'
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, scale=1, listw, listw2=NULL, control=list())
```

## Arguments

- object:

  A spatial regression model object fitted by maximum likelihood with
  [`spautolm`](https://r-spatial.github.io/spatialreg/reference/spautolm.md)

- mcmc:

  The number of MCMC iterations after burnin

- verbose:

  default NULL, use global option value; if TRUE, reports progress

- ...:

  Arguments passed through

- burnin:

  The number of burn-in iterations for the sampler

- scale:

  a positive scale parameter

- listw, listw2:

  `listw` objects created for example by `nb2listw`; should be the same
  object(s) used for fitting the model

- control:

  list of extra control arguments - see
  [`spautolm`](https://r-spatial.github.io/spatialreg/reference/spautolm.md)

## Value

An object of class “mcmc” suited to coda, with attributes: “accept”
acceptance rate; “type” input ML fitted model type “SAR”, “CAR”, “SMA”,
“lag”, “mixed”, “error”, “sac”, “sacmixed”; “timings” run times

## Note

If the acceptance rate is below 0.05, a warning will be issued; consider
increasing mcmc.

## References

Jim Albert (2007) Bayesian Computation with R, Springer, New York, pp.
104-105.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`rwmetrop`](https://rdrr.io/pkg/LearnBayes/man/rwmetrop.html),
[`spautolm`](https://r-spatial.github.io/spatialreg/reference/spautolm.md),
[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`errorsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`sacsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md)

## Examples

``` r
require("sf", quietly=TRUE)
nydata <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData")[1], quiet=TRUE)
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
suppressMessages(ID <- as.character(names(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1]))[-1]))
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))
#> [1] TRUE
#require("spdep", quietly=TRUE)
listw_NY <- spdep::mat2listw(nyadjmat, as.character(nydata$AREAKEY), style="B")
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
#> 
#> Call: 
#> spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, family = "SAR", method = "eigen")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.56754 -0.38239 -0.02643  0.33109  4.01219 
#> 
#> Coefficients: 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.618193   0.176784 -3.4969 0.0004707
#> PEXPOSURE    0.071014   0.042051  1.6888 0.0912635
#> PCTAGE65P    3.754200   0.624722  6.0094 1.862e-09
#> PCTOWNHOME  -0.419890   0.191329 -2.1946 0.0281930
#> 
#> Lambda: 0.040487 LR test value: 5.2438 p-value: 0.022026 
#> Numerical Hessian standard error of lambda: 0.017194 
#> 
#> Log likelihood: -276.1069 
#> ML residual variance (sigma squared): 0.41388, (sigma: 0.64333)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: 564.21
#> 
res <- MCMCsamp(esar1f, mcmc=1000, burnin=200, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:1000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 1000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.04312 0.01737 0.0005491       0.002192
#> (Intercept) -0.64875 0.18226 0.0057637       0.023218
#> PEXPOSURE    0.08481 0.04907 0.0015518       0.006311
#> PCTAGE65P    3.78234 0.58076 0.0183653       0.066333
#> PCTOWNHOME  -0.39620 0.22091 0.0069858       0.030267
#> 
#> 2. Quantiles for each variable:
#> 
#>                  2.5%      25%      50%     75%     97.5%
#> lambda       0.008417  0.03208  0.04214  0.0550  0.076103
#> (Intercept) -1.039139 -0.77553 -0.64401 -0.5121 -0.326717
#> PEXPOSURE   -0.004897  0.05362  0.08142  0.1152  0.202214
#> PCTAGE65P    2.627489  3.40004  3.73461  4.1834  4.961677
#> PCTOWNHOME  -0.819693 -0.56249 -0.39268 -0.2268  0.003843
#> 
# \dontrun{
esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1fw)
#> 
#> Call: 
#> spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, weights = POP8, family = "SAR", method = "eigen")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.48488 -0.26823  0.09489  0.46552  4.28343 
#> 
#> Coefficients: 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.797063   0.144054 -5.5331 3.146e-08
#> PEXPOSURE    0.080545   0.028334  2.8428  0.004473
#> PCTAGE65P    3.816731   0.576037  6.6258 3.453e-11
#> PCTOWNHOME  -0.380778   0.156507 -2.4330  0.014975
#> 
#> Lambda: 0.0095636 LR test value: 0.32665 p-value: 0.56764 
#> Numerical Hessian standard error of lambda: 0.016576 
#> 
#> Log likelihood: -251.6017 
#> ML residual variance (sigma squared): 1104.1, (sigma: 33.229)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: NA (not available for weighted model)
#> 
res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.01300 0.01674 0.0002368      0.0009254
#> (Intercept) -0.80039 0.14839 0.0020986      0.0084178
#> PEXPOSURE    0.08015 0.02889 0.0004086      0.0015790
#> PCTAGE65P    3.78384 0.58436 0.0082641      0.0333591
#> PCTOWNHOME  -0.37235 0.16125 0.0022804      0.0089901
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%       25%      50%      75%    97.5%
#> lambda      -0.02006  0.001984  0.01314  0.02452  0.04696
#> (Intercept) -1.07354 -0.900021 -0.80332 -0.70216 -0.51025
#> PEXPOSURE    0.02362  0.060131  0.07984  0.09844  0.13595
#> PCTAGE65P    2.61573  3.396785  3.81384  4.18217  4.86523
#> PCTOWNHOME  -0.69549 -0.479514 -0.37158 -0.26255 -0.06704
#> 
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
#> 
#> Call: 
#> spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, family = "CAR", method = "eigen")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.539732 -0.384311 -0.030646  0.335126  3.808848 
#> 
#> Coefficients: 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.648362   0.181129 -3.5796 0.0003442
#> PEXPOSURE    0.077899   0.043692  1.7829 0.0745986
#> PCTAGE65P    3.703830   0.627185  5.9055 3.516e-09
#> PCTOWNHOME  -0.382789   0.195564 -1.9574 0.0503053
#> 
#> Lambda: 0.084123 LR test value: 5.8009 p-value: 0.016018 
#> Numerical Hessian standard error of lambda: 0.03086 
#> 
#> Log likelihood: -275.8283 
#> ML residual variance (sigma squared): 0.40758, (sigma: 0.63842)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: 563.66
#> 
res <- MCMCsamp(ecar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.08485 0.03006 0.0004251       0.001840
#> (Intercept) -0.66320 0.21524 0.0030439       0.014580
#> PEXPOSURE    0.08242 0.04989 0.0007056       0.003037
#> PCTAGE65P    3.65269 0.64311 0.0090950       0.039341
#> PCTOWNHOME  -0.35760 0.22514 0.0031840       0.014849
#> 
#> 2. Quantiles for each variable:
#> 
#>                  2.5%      25%      50%     75%   97.5%
#> lambda       0.021098  0.06566  0.08492  0.1070  0.1398
#> (Intercept) -1.172915 -0.78604 -0.63980 -0.5233 -0.2865
#> PEXPOSURE   -0.008489  0.04920  0.07832  0.1125  0.1965
#> PCTAGE65P    2.400190  3.17776  3.65485  4.0906  4.9433
#> PCTOWNHOME  -0.761080 -0.51641 -0.37279 -0.2168  0.1504
#> 
esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1fw)
#> 
#> Call: 
#> spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, weights = POP8, family = "SAR", method = "eigen")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.48488 -0.26823  0.09489  0.46552  4.28343 
#> 
#> Coefficients: 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.797063   0.144054 -5.5331 3.146e-08
#> PEXPOSURE    0.080545   0.028334  2.8428  0.004473
#> PCTAGE65P    3.816731   0.576037  6.6258 3.453e-11
#> PCTOWNHOME  -0.380778   0.156507 -2.4330  0.014975
#> 
#> Lambda: 0.0095636 LR test value: 0.32665 p-value: 0.56764 
#> Numerical Hessian standard error of lambda: 0.016576 
#> 
#> Log likelihood: -251.6017 
#> ML residual variance (sigma squared): 1104.1, (sigma: 33.229)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: NA (not available for weighted model)
#> 
res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.01285 0.01605 0.0002269      0.0008353
#> (Intercept) -0.78979 0.16154 0.0022845      0.0098840
#> PEXPOSURE    0.07897 0.03125 0.0004420      0.0018876
#> PCTAGE65P    3.80337 0.62549 0.0088458      0.0373827
#> PCTOWNHOME  -0.39118 0.17236 0.0024375      0.0101411
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%       25%      50%      75%    97.5%
#> lambda      -0.01870  0.002134  0.01270  0.02419  0.04464
#> (Intercept) -1.09016 -0.895820 -0.78818 -0.68302 -0.46802
#> PEXPOSURE    0.02011  0.057760  0.07845  0.10179  0.14063
#> PCTAGE65P    2.51194  3.411864  3.78334  4.21963  5.01253
#> PCTOWNHOME  -0.72202 -0.511085 -0.39158 -0.27288 -0.05472
#> 
ecar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="CAR", method="eigen")
summary(ecar1fw)
#> 
#> Call: 
#> spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, weights = POP8, family = "CAR", method = "eigen")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.491042 -0.270906  0.081435  0.451556  4.198134 
#> 
#> Coefficients: 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.790154   0.144862 -5.4545 4.910e-08
#> PEXPOSURE    0.081922   0.028593  2.8651  0.004169
#> PCTAGE65P    3.825858   0.577720  6.6223 3.536e-11
#> PCTOWNHOME  -0.386820   0.157436 -2.4570  0.014010
#> 
#> Lambda: 0.022419 LR test value: 0.38785 p-value: 0.53343 
#> Numerical Hessian standard error of lambda: 0.038391 
#> 
#> Log likelihood: -251.5711 
#> ML residual variance (sigma squared): 1102.9, (sigma: 33.21)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: NA (not available for weighted model)
#> 
res <- MCMCsamp(ecar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.03478 0.04079 0.0005769       0.002511
#> (Intercept) -0.83251 0.16481 0.0023308       0.011061
#> PEXPOSURE    0.08962 0.03570 0.0005049       0.002349
#> PCTAGE65P    3.74934 0.58194 0.0082299       0.034899
#> PCTOWNHOME  -0.34673 0.17171 0.0024283       0.010821
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%       25%      50%      75%     97.5%
#> lambda      -0.04817  0.006844  0.03556  0.06383  0.112814
#> (Intercept) -1.19744 -0.939036 -0.82513 -0.72555 -0.511357
#> PEXPOSURE    0.02580  0.065603  0.08504  0.11214  0.164547
#> PCTAGE65P    2.63663  3.329179  3.75917  4.15136  4.903759
#> PCTOWNHOME  -0.67001 -0.454367 -0.35147 -0.22914 -0.001242
#> 
# }
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
#> 
#> Call:errorsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, 
#>     data = nydata, listw = listw_NY)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.56754 -0.38239 -0.02643  0.33109  4.01219 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.618193   0.176784 -3.4969 0.0004707
#> PEXPOSURE    0.071014   0.042051  1.6888 0.0912635
#> PCTAGE65P    3.754200   0.624722  6.0094 1.862e-09
#> PCTOWNHOME  -0.419890   0.191329 -2.1946 0.0281930
#> 
#> Lambda: 0.040487, LR test value: 5.2438, p-value: 0.022026
#> Asymptotic standard error: 0.016214
#>     z-value: 2.4971, p-value: 0.01252
#> Wald statistic: 6.2356, p-value: 0.01252
#> 
#> Log likelihood: -276.1069 for error model
#> ML residual variance (sigma squared): 0.41388, (sigma: 0.64333)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: 564.21, (AIC for lm: 567.46)
#> 
res <- MCMCsamp(esar0, mcmc=1000, burnin=200, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:1000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 1000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.04682 0.01895 0.0005993       0.002739
#> (Intercept) -0.65775 0.21896 0.0069243       0.031633
#> PEXPOSURE    0.08543 0.05260 0.0016634       0.007626
#> PCTAGE65P    3.65044 0.57371 0.0181422       0.060983
#> PCTOWNHOME  -0.36116 0.23745 0.0075089       0.036472
#> 
#> 2. Quantiles for each variable:
#> 
#>                   2.5%      25%      50%      75%    97.5%
#> lambda       0.0109828  0.03340  0.04839  0.05958  0.08139
#> (Intercept) -1.0771095 -0.81427 -0.62781 -0.51408 -0.23029
#> PEXPOSURE   -0.0004135  0.04738  0.07927  0.11564  0.20607
#> PCTAGE65P    2.5911599  3.20769  3.63005  4.09097  4.74100
#> PCTOWNHOME  -0.8037868 -0.54136 -0.35322 -0.20277  0.17058
#> 
# \dontrun{
esar0w <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8)
summary(esar0)
#> 
#> Call:errorsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, 
#>     data = nydata, listw = listw_NY)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.56754 -0.38239 -0.02643  0.33109  4.01219 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.618193   0.176784 -3.4969 0.0004707
#> PEXPOSURE    0.071014   0.042051  1.6888 0.0912635
#> PCTAGE65P    3.754200   0.624722  6.0094 1.862e-09
#> PCTOWNHOME  -0.419890   0.191329 -2.1946 0.0281930
#> 
#> Lambda: 0.040487, LR test value: 5.2438, p-value: 0.022026
#> Asymptotic standard error: 0.016214
#>     z-value: 2.4971, p-value: 0.01252
#> Wald statistic: 6.2356, p-value: 0.01252
#> 
#> Log likelihood: -276.1069 for error model
#> ML residual variance (sigma squared): 0.41388, (sigma: 0.64333)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: 564.21, (AIC for lm: 567.46)
#> 
res <- MCMCsamp(esar0w, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> lambda       0.01177 0.01571 0.0002221      0.0008781
#> (Intercept) -0.79575 0.14274 0.0020186      0.0078442
#> PEXPOSURE    0.08036 0.02929 0.0004143      0.0017545
#> PCTAGE65P    3.81500 0.57129 0.0080793      0.0330363
#> PCTOWNHOME  -0.38426 0.15508 0.0021932      0.0087813
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%      25%      50%      75%    97.5%
#> lambda      -0.01899  0.00175  0.01220  0.02213  0.04328
#> (Intercept) -1.07330 -0.89063 -0.79366 -0.70020 -0.51399
#> PEXPOSURE    0.02200  0.06082  0.08123  0.10067  0.13675
#> PCTAGE65P    2.70738  3.41885  3.80976  4.20150  4.94589
#> PCTOWNHOME  -0.70906 -0.47931 -0.38839 -0.28703 -0.06345
#> 
esar1 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, etype="emixed")
summary(esar1)
#> 
#> Call:errorsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, 
#>     data = nydata, listw = listw_NY, etype = "emixed")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.81562 -0.37641 -0.02224  0.33638  4.00054 
#> 
#> Type: error 
#> Coefficients: (asymptotic standard errors) 
#>                  Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)     -1.118019   0.247425 -4.5186 6.225e-06
#> PEXPOSURE        0.218279   0.079245  2.7545  0.005879
#> PCTAGE65P        3.416477   0.645587  5.2920 1.210e-07
#> PCTOWNHOME       0.036593   0.249835  0.1465  0.883551
#> lag.(Intercept)  0.121515   0.057636  2.1083  0.035003
#> lag.PEXPOSURE   -0.035075   0.015943 -2.2000  0.027808
#> lag.PCTAGE65P    0.263096   0.220118  1.1953  0.231989
#> lag.PCTOWNHOME  -0.155680   0.059213 -2.6291  0.008560
#> 
#> Lambda: 0.022723, LR test value: 1.6846, p-value: 0.19432
#> Asymptotic standard error: 0.017169
#>     z-value: 1.3235, p-value: 0.18567
#> Wald statistic: 1.7516, p-value: 0.18567
#> 
#> Log likelihood: -269.5398 for error model
#> ML residual variance (sigma squared): 0.39759, (sigma: 0.63055)
#> Number of observations: 281 
#> Number of parameters estimated: 10 
#> AIC: 559.08, (AIC for lm: 558.76)
#> 
res <- MCMCsamp(esar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                     Mean      SD  Naive SE Time-series SE
#> lambda           0.02791 0.01696 0.0002398       0.001257
#> (Intercept)     -1.11089 0.24295 0.0034359       0.019066
#> PEXPOSURE        0.20998 0.08066 0.0011407       0.006616
#> PCTAGE65P        3.44550 0.61861 0.0087484       0.045120
#> PCTOWNHOME       0.03249 0.25058 0.0035437       0.022380
#> lag.(Intercept)  0.11788 0.05850 0.0008273       0.004548
#> lag.PEXPOSURE   -0.03381 0.01626 0.0002299       0.001358
#> lag.PCTAGE65P    0.25893 0.23503 0.0033238       0.019521
#> lag.PCTOWNHOME  -0.15337 0.05824 0.0008236       0.004496
#> 
#> 2. Quantiles for each variable:
#> 
#>                      2.5%      25%      50%      75%      97.5%
#> lambda          -0.004792  0.01604  0.02829  0.03980  0.0616941
#> (Intercept)     -1.549205 -1.28823 -1.12499 -0.91520 -0.6509996
#> PEXPOSURE        0.063165  0.15363  0.20577  0.26455  0.3711420
#> PCTAGE65P        2.335975  3.01533  3.44864  3.91005  4.6644758
#> PCTOWNHOME      -0.505045 -0.11695  0.05420  0.19038  0.5171339
#> lag.(Intercept)  0.006027  0.07616  0.11689  0.16220  0.2310606
#> lag.PEXPOSURE   -0.064845 -0.04488 -0.03385 -0.02339 -0.0009086
#> lag.PCTAGE65P   -0.195440  0.09430  0.24778  0.43584  0.6912104
#> lag.PCTOWNHOME  -0.263176 -0.19213 -0.15386 -0.11592 -0.0317401
#> 
lsar0 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(lsar0)
#> 
#> Call:lagsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.586752 -0.391580 -0.022469  0.338017  4.029430 
#> 
#> Type: lag 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.514495   0.156154 -3.2948  0.000985
#> PEXPOSURE    0.047627   0.034509  1.3801  0.167542
#> PCTAGE65P    3.648198   0.599046  6.0900 1.129e-09
#> PCTOWNHOME  -0.414601   0.169554 -2.4453  0.014475
#> 
#> Rho: 0.038893, LR test value: 6.9683, p-value: 0.0082967
#> Asymptotic standard error: 0.015053
#>     z-value: 2.5837, p-value: 0.0097755
#> Wald statistic: 6.6754, p-value: 0.0097755
#> 
#> Log likelihood: -275.2447 for lag model
#> ML residual variance (sigma squared): 0.41166, (sigma: 0.6416)
#> Number of observations: 281 
#> Number of parameters estimated: 6 
#> AIC: 562.49, (AIC for lm: 567.46)
#> LM test for residual autocorrelation
#> test value: 1.4633, p-value: 0.22641
#> 
res <- MCMCsamp(lsar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD  Naive SE Time-series SE
#> rho          0.03924 0.01531 0.0002166      0.0009369
#> (Intercept) -0.51482 0.16721 0.0023647      0.0103798
#> PEXPOSURE    0.05057 0.03374 0.0004771      0.0019411
#> PCTAGE65P    3.58543 0.66827 0.0094508      0.0442173
#> PCTOWNHOME  -0.41083 0.18638 0.0026358      0.0118147
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%      25%      50%      75%    97.5%
#> rho          0.01129  0.02805  0.03879  0.05012  0.07010
#> (Intercept) -0.84632 -0.62884 -0.50852 -0.39309 -0.20739
#> PEXPOSURE   -0.01768  0.02886  0.04981  0.07205  0.11808
#> PCTAGE65P    2.28982  3.12634  3.60888  4.03227  4.88549
#> PCTOWNHOME  -0.74664 -0.54787 -0.40497 -0.28735 -0.03815
#> 
lsar1 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="mixed")
summary(lsar1)
#> 
#> Call:lagsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, type = "mixed")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.799308 -0.390125 -0.021371  0.346128  3.965251 
#> 
#> Type: mixed 
#> Coefficients: (asymptotic standard errors) 
#>                  Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)     -1.131233   0.249631 -4.5316 5.853e-06
#> PEXPOSURE        0.218364   0.079301  2.7536  0.005894
#> PCTAGE65P        3.361158   0.654123  5.1384 2.771e-07
#> PCTOWNHOME       0.071903   0.253967  0.2831  0.777085
#> lag.(Intercept)  0.132544   0.056175  2.3595  0.018300
#> lag.PEXPOSURE   -0.035239   0.015536 -2.2681  0.023322
#> lag.PCTAGE65P    0.161685   0.223690  0.7228  0.469798
#> lag.PCTOWNHOME  -0.140681   0.058529 -2.4036  0.016234
#> 
#> Rho: 0.026981, LR test value: 2.558, p-value: 0.10974
#> Asymptotic standard error: 0.016766
#>     z-value: 1.6093, p-value: 0.10755
#> Wald statistic: 2.5899, p-value: 0.10755
#> 
#> Log likelihood: -269.1031 for mixed model
#> ML residual variance (sigma squared): 0.39587, (sigma: 0.62918)
#> Number of observations: 281 
#> Number of parameters estimated: 10 
#> AIC: 558.21, (AIC for lm: 558.76)
#> LM test for residual autocorrelation
#> test value: 4.908, p-value: 0.026732
#> 
res <- MCMCsamp(lsar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                     Mean      SD  Naive SE Time-series SE
#> rho              0.02471 0.01629 0.0002304       0.001219
#> (Intercept)     -1.12989 0.24291 0.0034353       0.017687
#> PEXPOSURE        0.21704 0.08665 0.0012254       0.007356
#> PCTAGE65P        3.40551 0.61299 0.0086689       0.043441
#> PCTOWNHOME       0.04146 0.25683 0.0036321       0.019859
#> lag.(Intercept)  0.12868 0.05650 0.0007991       0.004535
#> lag.PEXPOSURE   -0.03472 0.01732 0.0002449       0.001526
#> lag.PCTAGE65P    0.17214 0.20763 0.0029364       0.015671
#> lag.PCTOWNHOME  -0.13618 0.06186 0.0008749       0.005326
#> 
#> 2. Quantiles for each variable:
#> 
#>                      2.5%      25%      50%      75%     97.5%
#> rho             -0.008938  0.01339  0.02554  0.03617  0.054650
#> (Intercept)     -1.604105 -1.28894 -1.12504 -0.97119 -0.660270
#> PEXPOSURE        0.053725  0.15881  0.21725  0.27236  0.389018
#> PCTAGE65P        2.213387  2.99625  3.42281  3.85452  4.551724
#> PCTOWNHOME      -0.451879 -0.14346  0.04721  0.22476  0.545593
#> lag.(Intercept)  0.019705  0.09016  0.13030  0.16767  0.243720
#> lag.PEXPOSURE   -0.067460 -0.04719 -0.03527 -0.02388  0.003145
#> lag.PCTAGE65P   -0.238603  0.04557  0.17728  0.29777  0.562834
#> lag.PCTOWNHOME  -0.246111 -0.18197 -0.13460 -0.09357 -0.018730
#> 
ssar0 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(ssar0)
#> 
#> Call:sacsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.468382 -0.375687 -0.034996  0.314714  3.833950 
#> 
#> Type: sac 
#> Coefficients: (asymptotic standard errors) 
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.386572   0.123188 -3.1381  0.001701
#> PEXPOSURE    0.026684   0.024013  1.1112  0.266479
#> PCTAGE65P    3.089824   0.562851  5.4896 4.029e-08
#> PCTOWNHOME  -0.323052   0.137449 -2.3503  0.018756
#> 
#> Rho: 0.089451
#> Asymptotic standard error: 0.019427
#>     z-value: 4.6046, p-value: 4.1325e-06
#> Lambda: -0.08192
#> Asymptotic standard error: 0.033201
#>     z-value: -2.4674, p-value: 0.01361
#> 
#> LR test value: 10.114, p-value: 0.0063661
#> 
#> Log likelihood: -273.672 for sac model
#> ML residual variance (sigma squared): 0.3766, (sigma: 0.61368)
#> Number of observations: 281 
#> Number of parameters estimated: 7 
#> AIC: 561.34, (AIC for lm: 567.46)
#> 
res <- MCMCsamp(ssar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean      SD Naive SE Time-series SE
#> rho         -0.04897 0.07268 0.001028        0.02171
#> lambda       0.07158 0.07109 0.001005        0.01911
#> (Intercept) -0.79253 0.30171 0.004267        0.05929
#> PEXPOSURE    0.12518 0.08005 0.001132        0.01428
#> PCTAGE65P    3.33373 0.62239 0.008802        0.04114
#> PCTOWNHOME  -0.24558 0.24700 0.003493        0.03231
#> 
#> 2. Quantiles for each variable:
#> 
#>                  2.5%      25%      50%        75%    97.5%
#> rho         -0.148205 -0.10220 -0.07484 -0.0003653  0.09879
#> lambda      -0.097739  0.03709  0.10543  0.1218063  0.13809
#> (Intercept) -1.371833 -1.01584 -0.80342 -0.5327538 -0.26199
#> PEXPOSURE   -0.001867  0.05674  0.12102  0.1859304  0.27450
#> PCTAGE65P    2.139989  2.91933  3.33752  3.7143648  4.58803
#> PCTOWNHOME  -0.743597 -0.40260 -0.25567 -0.0792356  0.24964
#> 
ssar1 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="sacmixed")
summary(ssar1)
#> 
#> Call:sacsarlm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = nydata, 
#>     listw = listw_NY, type = "sacmixed")
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -1.633958 -0.363826 -0.019927  0.348238  3.655509 
#> 
#> Type: sacmixed 
#> Coefficients: (asymptotic standard errors) 
#>                  Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)     -1.133298   0.247495 -4.5791 4.670e-06
#> PEXPOSURE        0.206963   0.074480  2.7788  0.005456
#> PCTAGE65P        3.083983   0.671081  4.5955 4.316e-06
#> PCTOWNHOME       0.174800   0.256280  0.6821  0.495196
#> lag.(Intercept)  0.153427   0.050817  3.0192  0.002534
#> lag.PEXPOSURE   -0.033400   0.013817 -2.4173  0.015634
#> lag.PCTAGE65P   -0.079738   0.222144 -0.3589  0.719634
#> lag.PCTOWNHOME  -0.102502   0.056760 -1.8059  0.070940
#> 
#> Rho: 0.092495
#> Asymptotic standard error: 0.023829
#>     z-value: 3.8817, p-value: 0.00010375
#> Lambda: -0.091069
#> Asymptotic standard error: 0.038431
#>     z-value: -2.3697, p-value: 0.017804
#> 
#> LR test value: 22.379, p-value: 0.0010335
#> 
#> Log likelihood: -267.5392 for sacmixed model
#> ML residual variance (sigma squared): 0.35617, (sigma: 0.5968)
#> Number of observations: 281 
#> Number of parameters estimated: 11 
#> AIC: 557.08, (AIC for lm: 567.46)
#> 
res <- MCMCsamp(ssar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#> 
#> Iterations = 1:5000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 5000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                      Mean      SD  Naive SE Time-series SE
#> rho             -0.005935 0.06247 0.0008835       0.012391
#> lambda           0.025149 0.06578 0.0009302       0.012564
#> (Intercept)     -1.104966 0.25946 0.0036694       0.020691
#> PEXPOSURE        0.220056 0.08072 0.0011415       0.006564
#> PCTAGE65P        3.448064 0.68477 0.0096841       0.053253
#> PCTOWNHOME       0.002352 0.25423 0.0035954       0.018914
#> lag.(Intercept)  0.107623 0.06651 0.0009406       0.007925
#> lag.PEXPOSURE   -0.034004 0.01609 0.0002276       0.001352
#> lag.PCTAGE65P    0.279208 0.31780 0.0044943       0.037387
#> lag.PCTOWNHOME  -0.134012 0.06451 0.0009123       0.005533
#> 
#> 2. Quantiles for each variable:
#> 
#>                     2.5%      25%       50%      75%     97.5%
#> rho             -0.10755 -0.05711 -0.014110  0.05341  0.099104
#> lambda          -0.09950 -0.03206  0.037324  0.08249  0.118177
#> (Intercept)     -1.60962 -1.28203 -1.096564 -0.93087 -0.597130
#> PEXPOSURE        0.05429  0.16818  0.217623  0.27874  0.368060
#> PCTAGE65P        2.11661  3.00219  3.455553  3.93036  4.829667
#> PCTOWNHOME      -0.46367 -0.19124 -0.002169  0.17287  0.544393
#> lag.(Intercept) -0.02153  0.05985  0.109887  0.15841  0.231289
#> lag.PEXPOSURE   -0.06423 -0.04610 -0.033598 -0.02331 -0.001144
#> lag.PCTAGE65P   -0.29223  0.02822  0.266936  0.50401  0.926973
#> lag.PCTOWNHOME  -0.26326 -0.17692 -0.130161 -0.08812 -0.009731
#> 
# }
```
