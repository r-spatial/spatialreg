# Prediction for spatial simultaneous autoregressive linear model objects

`predict.Sarlm()` calculates predictions as far as is at present
possible for for spatial simultaneous autoregressive linear model
objects, using Haining's terminology for decomposition into trend,
signal, and noise, or other types of predictors — see references.

## Usage

``` r
# S3 method for class 'Sarlm'
predict(object, newdata = NULL, listw = NULL, pred.type = "TS", all.data = FALSE,
 zero.policy = NULL, legacy = TRUE, legacy.mixed = FALSE, power = NULL, order = 250,
 tol = .Machine$double.eps^(3/5), spChk = NULL, ...)
#\method{predict}{SLX}(object, newdata, listw, zero.policy=NULL, ...)
# S3 method for class 'Sarlm.pred'
print(x, ...)
# S3 method for class 'Sarlm.pred'
as.data.frame(x, ...)
```

## Arguments

- object:

  `Sarlm` object returned by `lagsarlm`, `errorsarlm` or `sacsarlm`, the
  method for SLX objects takes the output of `lmSLX`

- newdata:

  data frame in which to predict — if NULL, predictions are for the data
  on which the model was fitted. Should have row names corresponding to
  region.id. If row names are exactly the same than the ones used for
  training, it uses in-sample predictors for forecast. See ‘Details’

- listw:

  a `listw` object created for example by `nb2listw`. In the
  out-of-sample prediction case (ie. if newdata is not NULL), if
  `legacy.mixed=FALSE` or if `pred.type!="TS"`, it should include both
  in-sample and out-of-sample spatial units. In this case, if regions of
  the listw are not in the correct order, they are reordered. See
  ‘Details’

- pred.type:

  predictor type — default “TS”, use decomposition into trend, signal,
  and noise ; other types available depending on `newdata`. If
  `newdata=NULL` (in-sample prediction), “TS”, “trend”, “TC” and “BP”
  are available. If `newdata` is not NULL and its row names are the same
  than the `data` used to fit the model (forecast case), “TS”, “trend”
  and “TC” are available. In other cases (out-of-sample prediction),
  “TS”, “trend”, “KP1”, “KP2”, “KP3”, “KP4”, “KP5”, “TC”, “BP”, “BPW”,
  “BPN”, “TS1”, “TC1”, “BP1”, “BPW1” and “BPN1” are available. See
  ‘Details’ and references

- all.data:

  (only applies to `pred.type="TC"` and newdata is not NULL) default
  FALSE: return predictions only for newdata units, if TRUE return
  predictions for all data units. See ‘Details’

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE (default) assign
  NA - causing the function to terminate with an error

- legacy:

  (only applies to lag and Durbin (mixed) models for `pred.type="TS"`)
  default TRUE: use ad-hoc predictor, if FALSE use DGP-based predictor

- legacy.mixed:

  (only applies to mixed models if newdata is not NULL) default FALSE:
  compute lagged variables from both in-sample and out-of-sample units
  with \\\[W X\]\_O\\ and \\\[W X\]\_S\\ where `X=cbind(Xs, Xo)`, if
  TRUE compute lagged variables independantly between in-sample and
  out-of-sample units with \\W\_{OO} X_O\\ and \\W\_{SS} X_S\\

- power:

  (only applies to lag and Durbin (mixed) models for “TS”, “KP1”, “KP2”,
  “KP3”, “TC”, “TC1”, “BP”, “BP1”, “BPN”, “BPN1”, “BPW” and “BPW1”
  types) use `powerWeights`, if default NULL, set FALSE if
  `object$method` is “eigen”, otherwise TRUE

- order:

  power series maximum limit if `power` is TRUE

- tol:

  tolerance for convergence of power series if `power` is TRUE

- spChk:

  should the row names of data frames be checked against the spatial
  objects for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)

- x:

  the object to be printed

&nbsp;

- ...:

  further arguments passed through

## Details

The function supports three types of prediction. In-sample prediction is
the computation of predictors on the data used to fit the model
(`newdata=NULL`). Prevision, also called forecast, is the computation of
some predictors (“trend”, in-sample “TC” and out-of-sample “TS”) on the
same spatial units than the ones used to fit the model, but with
different observations of the variables in the model (row names of
`newdata` should have the same row names than the data frame used to fit
the model). And out-of-sample prediction is the computation of
predictors on other spatial units than the ones used to fit the model
(`newdata` has different row names). For extensive definitions, see
Goulard et al. (2017).

`pred.type` of predictors are available according to the model of
`object` an to the type of prediction. In the two following tables,
“yes” means that the predictor can be used with the model, “no” means
that `predict.Sarlm()` will stop with an error, and “yes\*” means that
the predictor is not designed for the specified model, but it can be
used with `predict.Sarlm()`. In the last case, be careful with the
computation of a inappropriate predictor.

*In-sample predictors by models*

|           |             |             |             |
|-----------|-------------|-------------|-------------|
| pred.type | sem (mixed) | lag (mixed) | sac (mixed) |
|           |             |             |             |
| “trend”   | yes         | yes         | yes         |
| “TS”      | yes         | yes         | no          |
| “TC”      | no          | yes         | yes\*       |
| “BP”      | no          | yes         | yes\*       |

Note that only “trend” and “TC” are available for prevision.

*Out-of-sample predictors by models*

|                |             |             |             |
|----------------|-------------|-------------|-------------|
| pred.type      | sem (mixed) | lag (mixed) | sac (mixed) |
|                |             |             |             |
| “trend”        | yes         | yes         | yes         |
| “TS”           | yes         | yes         | no          |
| “TS1” or “KP4” | no          | yes         | yes         |
| “TC”           | no          | yes         | yes\*       |
| “TC1” or “KP1” | yes         | yes         | yes         |
| “BP”           | no          | yes         | yes\*       |
| “BP1”          | no          | yes         | yes\*       |
| “BPW”          | no          | yes         | yes\*       |
| “BPW1”         | no          | yes         | yes\*       |
| “BN”           | no          | yes         | yes\*       |
| “BPN1”         | no          | yes         | yes\*       |
| “KP2”          | yes         | yes         | yes         |
| “KP3”          | yes         | yes         | yes         |
| “KP5”          | yes         | no          | yes\*       |

Values for `pred.type=` include “TS1”, “TC”, “TC1”, “BP”, “BP1”, “BPW”,
“BPW1”, “BPN”, “BPN1”, following the notation in Goulard et al. (2017),
and for `pred.type=` “KP1”, “KP2”, “KP3”, “KP4”, “KP5”, following the
notation in Kelejian et al. (2007). `pred.type="TS"` is described bellow
and in Bivand (2002).

In the following, the trend is the non-spatial smooth, the signal is the
spatial smooth, and the noise is the residual. The fit returned by
`pred.type="TS"` is the sum of the trend and the signal.

When `pred.type="TS"`, the function approaches prediction first by
dividing invocations between those with or without newdata. When no
newdata is present, the response variable may be reconstructed as the
sum of the trend, the signal, and the noise (residuals). Since the
values of the response variable are known, their spatial lags are used
to calculate signal components (Cressie 1993, p. 564). For the error
model, trend = \\X \beta\\, and signal = \\\lambda W y - \lambda W X
\beta\\. For the lag and mixed models, trend = \\X \beta\\, and signal =
\\\rho W y\\.

This approach differs from the design choices made in other software,
for example GeoDa, which does not use observations of the response
variable, and corresponds to the newdata situation described below.

When however newdata is used for prediction, no observations of the
response variable being predicted are available. Consequently, while the
trend components are the same, the signal cannot take full account of
the spatial smooth. In the error model and Durbin error model, the
signal is set to zero, since the spatial smooth is expressed in terms of
the error: \\(I - \lambda W)^{-1} \varepsilon\\.

In the lag model, the signal can be expressed in the following way (for
legacy=TRUE):

\$\$(I - \rho W) y = X \beta + \varepsilon\$\$ \$\$y = (I - \rho W)^{-1}
X \beta + (I - \rho W)^{-1} \varepsilon\$\$

giving a feasible signal component of:

\$\$\rho W y = \rho W (I - \rho W)^{-1} X \beta\$\$

For legacy=FALSE, the trend is computed first as:

\$\$X \beta\$\$

next the prediction using the DGP:

\$\$(I - \rho W)^{-1} X \beta\$\$

and the signal is found as the difference between prediction and trend.
The numerical results for the legacy and DGP methods are identical.

setting the error term to zero. This also means that predictions of the
signal component for lag and mixed models require the inversion of an
n-by-n matrix.

Because the outcomes of the spatial smooth on the error term are
unobservable, this means that the signal values for newdata are
incomplete. In the mixed model, the spatially lagged RHS variables
influence both the trend and the signal, so that the root mean square
prediction error in the examples below for this case with newdata is
smallest, although the model was not the best fit.

If `newdata` has more than one row, leave-one-out predictors
(`pred.type=` include “TS1”, “TC1”, “BP1”, “BPW1”, “BPN1”, “KP1”, “KP2”,
“KP3”, “KP4”, “KP5”) are computed separatly on each out-of-sample unit.

`listw` should be provided except if `newdata=NULL` and `pred.type=`
include “TS”, “trend”, or if `newdata` is not `NULL`,
`pred.type="trend"` and `object` is not a mixed model.

`all.data` is useful when some out-of-sample predictors return different
predictions for in-sample units, than the same predictor type computed
only on in-sample data.

## Value

`predict.Sarlm()` returns a vector of predictions with three attribute
vectors of trend, signal (only for `pred.type="TS"`) and region.id
values and two other attributes of pred.type and call with class
`Sarlm.pred`.

`print.Sarlm.pred()` is a print function for this class, printing and
returning a data frame with columns: "fit", "trend" and "signal" (when
available) and with region.id as row names.

## References

Haining, R. 1990 *Spatial data analysis in the social and environmental
sciences*, Cambridge: Cambridge University Press, p. 258; Cressie, N. A.
C. 1993 *Statistics for spatial data*, Wiley, New York; Michel Goulard,
Thibault Laurent & Christine Thomas-Agnan, 2017 *About predictions in
spatial autoregressive models: optimal and almost optimal strategies*,
Spatial Economic Analysis Volume 12, Issue 2–3, 304–325
[doi:10.1080/17421772.2017.1300679](https://doi.org/10.1080/17421772.2017.1300679)
, ; Kelejian, H. H. and Prucha, I. R. 2007 *The relative efficiencies of
various predictors in spatial econometric models containing spatial
lags*, Regional Science and Urban Economics, Volume 37, Issue 3,
363–374; Bivand, R. 2002 *Spatial econometrics functions in R: Classes
and methods*, Journal of Geographical Systems, Volume 4, No. 4, 405–421

## Author

Roger Bivand <Roger.Bivand@nhh.no> and Martin Gubri

## See also

[`errorsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`sacsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md)

## Examples

``` r
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)

COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
  type="mixed")
print(p1 <- predict(COL.mix.eig))
#> This method assumes the response is known - see manual page
#>            fit      trend    signal
#> 1001 26.044311 14.8543508 11.189960
#> 1002 44.034234 29.2632112 14.771023
#> 1003 43.511934 25.8193818 17.692553
#> 1004 37.656561 16.4555583 21.201002
#> 1005 10.902976  0.3664066 10.536570
#> 1006 36.829798 24.2905246 12.539274
#> 1007 44.290467 27.0386615 17.251806
#> 1008 38.853571 21.5342393 17.319331
#> 1009 50.870854 29.5092783 21.361576
#> 1010 16.401300  5.6029104 10.798389
#> 1011 36.354390 28.6415353  7.712855
#> 1012 20.452836 12.4607277  7.992108
#> 1013 20.324088 14.4173433  5.906745
#> 1014 19.243496 10.2606419  8.982854
#> 1015 19.747775 12.2556861  7.492089
#> 1016  6.962527 -2.0137491  8.976276
#> 1017  7.452143 -6.3808928 13.833036
#> 1018 28.481587 14.2125594 14.269028
#> 1019 43.351392 28.0442064 15.307186
#> 1020 50.359682 30.6608153 19.698867
#> 1021 38.905226 24.7490977 14.156128
#> 1022 44.724478 28.8314299 15.893048
#> 1023 37.888974 23.7778863 14.111087
#> 1024 45.527017 26.9163190 18.610698
#> 1025 32.429571 17.1892401 15.240331
#> 1026 26.490842 14.8893980 11.601444
#> 1027 35.629158 23.4577209 12.171437
#> 1028 35.574326 21.9001006 13.674226
#> 1029 38.598639 23.1818442 15.416795
#> 1030 36.602053 14.8614072 21.740646
#> 1031 50.320031 30.1013982 20.218633
#> 1032 53.698863 31.2094168 22.489447
#> 1033 49.364208 26.5151201 22.849088
#> 1034 46.262357 25.5538226 20.708534
#> 1035 39.177121 15.7689329 23.408188
#> 1036 54.984344 32.6590841 22.325260
#> 1037 51.611458 33.1290203 18.482438
#> 1038 51.998831 30.7428313 21.256000
#> 1039 43.651605 27.4107880 16.240817
#> 1040 44.196841 25.9409252 18.255916
#> 1041 49.310592 29.5106497 19.799943
#> 1042 37.995310 15.7039024 22.291408
#> 1043 46.908709 28.2687603 18.639948
#> 1044 28.976789 19.5389223  9.437867
#> 1045 25.343793 17.1838200  8.159973
#> 1046 24.006252 16.0103703  7.995882
#> 1047 25.034907 18.4616473  6.573260
#> 1048 10.478529  3.3573578  7.121171
#> 1049 13.495623  5.3356502  8.159973
print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE))
#>            fit      trend    signal
#> 1001 29.038788 14.8543508 14.184437
#> 1002 46.227075 29.2632112 16.963864
#> 1003 45.640479 25.8193818 19.821097
#> 1004 36.643520 16.4555583 20.187962
#> 1005 14.819940  0.3664066 14.453533
#> 1006 38.764777 24.2905246 14.474252
#> 1007 45.715716 27.0386615 18.677055
#> 1008 37.514611 21.5342393 15.980372
#> 1009 49.324228 29.5092783 19.814950
#> 1010 17.510607  5.6029104 11.907696
#> 1011 34.973608 28.6415353  6.332072
#> 1012 21.079100 12.4607277  8.618372
#> 1013 19.704134 14.4173433  5.286791
#> 1014 16.365521 10.2606419  6.104879
#> 1015 17.063856 12.2556861  4.808170
#> 1016  6.190282 -2.0137491  8.204031
#> 1017  5.967260 -6.3808928 12.348153
#> 1018 29.250462 14.2125594 15.037902
#> 1019 41.530036 28.0442064 13.485830
#> 1020 49.344770 30.6608153 18.683954
#> 1021 39.508818 24.7490977 14.759720
#> 1022 42.772692 28.8314299 13.941262
#> 1023 37.114901 23.7778863 13.337015
#> 1024 43.622499 26.9163190 16.706180
#> 1025 33.247197 17.1892401 16.057957
#> 1026 30.301331 14.8893980 15.411933
#> 1027 38.316063 23.4577209 14.858342
#> 1028 36.886068 21.9001006 14.985967
#> 1029 38.970564 23.1818442 15.788720
#> 1030 33.014615 14.8614072 18.153208
#> 1031 48.209875 30.1013982 18.108477
#> 1032 50.808064 31.2094168 19.598647
#> 1033 44.555996 26.5151201 18.040876
#> 1034 43.232773 25.5538226 17.678951
#> 1035 35.009061 15.7689329 19.240128
#> 1036 52.113364 32.6590841 19.454280
#> 1037 52.189015 33.1290203 19.059995
#> 1038 51.631805 30.7428313 20.888973
#> 1039 46.543565 27.4107880 19.132776
#> 1040 45.036095 25.9409252 19.095170
#> 1041 45.907835 29.5106497 16.397185
#> 1042 35.337110 15.7039024 19.633208
#> 1043 43.948398 28.2687603 15.679638
#> 1044 32.091257 19.5389223 12.552334
#> 1045 29.647005 17.1838200 12.463185
#> 1046 26.375304 16.0103703 10.364934
#> 1047 27.235807 18.4616473  8.774160
#> 1048 14.785518  3.3573578 11.428160
#> 1049 17.798835  5.3356502 12.463185
AIC(COL.mix.eig)
#> [1] 376.787
sqrt(deviance(COL.mix.eig)/length(COL.nb))
#> [1] 9.580773
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
#> [1] 9.580773
sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))
#> [1] 10.35029

COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
AIC(COL.err.eig)
#> [1] 376.7609
sqrt(deviance(COL.err.eig)/length(COL.nb))
#> [1] 9.776221
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
#> This method assumes the response is known - see manual page
#> [1] 9.776221
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))
#> [1] 11.61744

COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
 etype="emixed")
AIC(COL.SDerr.eig)
#> [1] 377.1693
sqrt(deviance(COL.SDerr.eig)/length(COL.nb))
#> [1] 9.619298
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig)))^2)/length(COL.nb))
#> This method assumes the response is known - see manual page
#> [1] 9.619298
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))
#> Warning: only legacy.mixed=TRUE is supported for pred.type='TS' and mixed models. legacy.mixed is forced
#> [1] 10.40472

AIC(COL.lag.eig)
#> [1] 374.7809
sqrt(deviance(COL.lag.eig)/length(COL.nb))
#> [1] 9.772129
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
#> This method assumes the response is known - see manual page
#> [1] 9.772129
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))
#> [1] 10.72654

p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
all.equal(p2, p3, check.attributes=FALSE)
#> [1] TRUE
p4 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, power=TRUE, legacy.mixed = TRUE)
all.equal(p2, p4, check.attributes=FALSE)
#> [1] TRUE
p5 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=TRUE, power=TRUE, legacy.mixed = TRUE)
all.equal(p2, p5, check.attributes=FALSE)
#> [1] TRUE
```
