# Moran eigenvector GLM filtering

The Moran eigenvector filtering function is intended to remove spatial
autocorrelation from the residuals of generalised linear models. It uses
brute force eigenvector selection to reach a subset of such vectors to
be added to the RHS of the GLM model to reduce residual autocorrelation
to below the specified alpha value. Since eigenvector selection only
works on symmetric weights, the weights are made symmetric before the
eigenvectors are found (from spdep 0.5-50).

## Usage

``` r
ME(formula, data=list(), family = gaussian, weights, offset,
 na.action=na.fail,listw=NULL, alpha=0.05, nsim=99, verbose=NULL,
 stdev=FALSE, zero.policy=NULL)
```

## Arguments

- formula:

  a symbolic description of the model to be fit

- data:

  an optional data frame containing the variables in the model

- family:

  a description of the error distribution and link function to be used
  in the model

- weights:

  an optional vector of weights to be used in the fitting process

- offset:

  this can be used to specify an a priori known component to be included
  in the linear predictor during fitting

- na.action:

  a function (default `options("na.action")`), can also be `na.omit` or
  `na.exclude` with consequences for residuals and fitted values - in
  these cases the spatial weights list will be subsetted to remove NAs
  in the data. It may be necessary to set zero.policy to TRUE because
  this subsetting may create no-neighbour observations. Note that only
  weights lists created without using the glist argument to `nb2listw`
  may be subsetted.

- listw:

  a `listw` object created for example by `nb2listw`

- alpha:

  used as a stopping rule to choose all eigenvectors up to and including
  the one with a p-value exceeding alpha

- nsim:

  number of permutations for permutation bootstrap for finding p-values

- verbose:

  default NULL, use global option value; if TRUE report eigenvectors
  selected

- stdev:

  if TRUE, p-value calculated from bootstrap permutation standard
  deviate using `pnorm` with alternative="greater", if FALSE the
  Hope-type p-value

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

## Details

The eigenvectors for inclusion are chosen by calculating the empirical
Moran's I values for the initial model plus each of the doubly centred
symmetric spatial weights matrix eigenvectors in turn. Then the first
eigenvector is chosen as that with the lowest Moran's I value. The
procedure is repeated until the lowest remaining Moran's I value has a
permutation-based probability value above alpha. The probability value
is either Hope-type or based on using the mean and standard deviation of
the permutations to calculate ZI based on the stdev argument.

## Value

An object of class `Me_res`:

- selection:

  a matrix summarising the selection of eigenvectors for inclusion, with
  columns:

  Eigenvector

  :   number of selected eigenvector

  ZI

  :   permutation-based standardized deviate of Moran's I if stdev=TRUE

  pr(ZI)

  :   probability value: if stdev=TRUE of the permutation-based
      standardized deviate, if FALSE the Hope-type probability value, in
      both cases on-sided

  The first row is the value at the start of the search

- vectors:

  a matrix of the selected eigenvectors in order of selection

## References

Dray S, Legendre P and Peres-Neto PR (2005) Spatial modeling: a
comprehensive framework for principle coordinate analysis of neigbbor
matrices (PCNM), Ecological Modelling; Griffith DA and Peres-Neto PR
(2006) Spatial modeling in ecology: the flexibility of eigenfunction
spatial analyses.

## Author

Roger Bivand and Pedro Peres-Neto

## See also

[`SpatialFiltering`](https://r-spatial.github.io/spatialreg/reference/SpatialFiltering.md),
[`glm`](https://rdrr.io/r/stats/glm.html)

## Examples

``` r
#require("spdep", quietly=TRUE)
data(hopkins, package="spData")
hopkins_part <- hopkins[21:36,36:21]
hopkins_part[which(hopkins_part > 0, arr.ind=TRUE)] <- 1
hopkins.rook.nb <- spdep::cell2nb(16, 16, type="rook")
glmbase <- glm(c(hopkins_part) ~ 1, family="binomial")
lw <- spdep::nb2listw(hopkins.rook.nb, style="B")
set.seed(123)
system.time(MEbinom1 <- ME(c(hopkins_part) ~ 1, family="binomial",
 listw=lw, alpha=0.05, verbose=TRUE, nsim=49))
#> eV[,1], I: 0.08273151 ZI: NA, pr(ZI): 0.04
#> eV[,9], I: 0.06266473 ZI: NA, pr(ZI): 0.14
#>    user  system elapsed 
#>   1.237   0.001   1.248 
glmME <- glm(c(hopkins_part) ~ 1 + fitted(MEbinom1), family="binomial")
#anova(glmME, test="Chisq")
coef(summary(glmME))
#>                       Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept)          -1.149323  0.1546365 -7.432422 1.066270e-13
#> fitted(MEbinom1)vec1  8.328982  2.4564041  3.390721 6.970892e-04
#> fitted(MEbinom1)vec9  5.493335  2.4015917  2.287372 2.217409e-02
anova(glmbase, glmME, test="Chisq")
#> Analysis of Deviance Table
#> 
#> Model 1: c(hopkins_part) ~ 1
#> Model 2: c(hopkins_part) ~ 1 + fitted(MEbinom1)
#>   Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#> 1       255     292.23                          
#> 2       253     274.81  2   17.413 0.0001655 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# \dontrun{
require("sf", quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
lw <- spdep::nb2listw(col.gal.nb)
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", alpha=0.1, verbose=TRUE)
#> Step 0 SelEvec 0 MinMi 0.2123742 ZMinMi 2.681 Pr(ZI) 0.007340246 
#> Step 1 SelEvec 6 MinMi 0.1178225 ZMinMi 1.84512 Pr(ZI) 0.06502014 
#> Step 2 SelEvec 4 MinMi 0.06242664 ZMinMi 1.494821 Pr(ZI) 0.1349611 
lagcol
#>   Step SelEvec      Eval      MinMi   ZMinMi      Pr(ZI)        R2     gamma
#> 0    0       0 0.0000000 0.21237415 2.681000 0.007340246 0.5524040   0.00000
#> 1    1       6 0.7161123 0.11782248 1.845120 0.065020139 0.6038801  25.46181
#> 2    2       4 0.8682938 0.06242664 1.494821 0.134961136 0.6531288 -26.68319
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
anova(lmbase, lmlag)
#> Analysis of Variance Table
#> 
#> Model 1: CRIME ~ INC + HOVAL
#> Model 2: CRIME ~ INC + HOVAL + fitted(lagcol)
#>   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
#> 1     46 6014.9                                
#> 2     44 4661.3  2    1353.6 6.3884 0.003666 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
set.seed(123)
system.time(lagcol1 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=lw, alpha=0.1, verbose=TRUE))
#> eV[,6], I: 0.1178225 ZI: NA, pr(ZI): 0.08
#> eV[,4], I: 0.06242664 ZI: NA, pr(ZI): 0.27
#>    user  system elapsed 
#>   0.566   0.000   0.569 
lagcol1
#>   Eigenvector ZI pr(ZI)
#> 0          NA NA   0.01
#> 1           6 NA   0.08
#> 2           4 NA   0.27
lmlag1 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol1), data=columbus)
anova(lmbase, lmlag1)
#> Analysis of Variance Table
#> 
#> Model 1: CRIME ~ INC + HOVAL
#> Model 2: CRIME ~ INC + HOVAL + fitted(lagcol1)
#>   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
#> 1     46 6014.9                                
#> 2     44 4661.3  2    1353.6 6.3884 0.003666 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

set.seed(123)
lagcol2 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=lw, alpha=0.1, stdev=TRUE, verbose=TRUE)
#> eV[,6], I: 0.1178225 ZI: 1.5509, pr(ZI): 0.06046283
#> eV[,4], I: 0.06242664 ZI: 0.681174, pr(ZI): 0.2478807
lagcol2
#>   Eigenvector       ZI      pr(ZI)
#> 0          NA 2.351591 0.009346653
#> 1           6 1.550900 0.060462832
#> 2           4 0.681174 0.247880696
lmlag2 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol2), data=columbus)
anova(lmbase, lmlag2)
#> Analysis of Variance Table
#> 
#> Model 1: CRIME ~ INC + HOVAL
#> Model 2: CRIME ~ INC + HOVAL + fitted(lagcol2)
#>   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
#> 1     46 6014.9                                
#> 2     44 4661.3  2    1353.6 6.3884 0.003666 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
NA.columbus <- columbus
NA.columbus$CRIME[20:25] <- NA
COL.ME.NA <- ME(CRIME ~ INC + HOVAL, data=NA.columbus, family="gaussian",
 listw=lw, alpha=0.1, stdev=TRUE, verbose=TRUE,
 na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> eV[,8], I: 0.1426723 ZI: 1.483169, pr(ZI): 0.06901474
#> eV[,1], I: 0.09838877 ZI: 0.9862904, pr(ZI): 0.1619953
COL.ME.NA$na.action
#> 20 21 22 23 24 25 
#> 20 21 22 23 24 25 
#> attr(,"class")
#> [1] "exclude"
summary(lm(CRIME ~ INC + HOVAL + fitted(COL.ME.NA), data=NA.columbus,
 na.action=na.exclude))
#> 
#> Call:
#> lm(formula = CRIME ~ INC + HOVAL + fitted(COL.ME.NA), data = NA.columbus, 
#>     na.action = na.exclude)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -30.1382  -6.0105   0.4095   7.1504  19.9399 
#> 
#> Coefficients:
#>                    Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)        66.92248    5.28663  12.659 3.33e-15 ***
#> INC                -1.40484    0.35678  -3.938  0.00034 ***
#> HOVAL              -0.30446    0.09831  -3.097  0.00366 ** 
#> fitted(COL.ME.NA)1 29.69422   10.58481   2.805  0.00788 ** 
#> fitted(COL.ME.NA)2 26.61612   11.29187   2.357  0.02367 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 10.48 on 38 degrees of freedom
#>   (6 observations deleted due to missingness)
#> Multiple R-squared:  0.6294, Adjusted R-squared:  0.5904 
#> F-statistic: 16.13 on 4 and 38 DF,  p-value: 8.353e-08
#> 
nc.sids <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
rn <- as.character(nc.sids$FIPS)
ncCC89_nb <- spdep::read.gal(system.file("weights/ncCC89.gal", package="spData")[1],
 region.id=rn)
#> Warning: neighbour object has 3 sub-graphs
ncCR85_nb <- spdep::read.gal(system.file("weights/ncCR85.gal", package="spData")[1],
 region.id=rn)
glmbase <- glm(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson")
set.seed(123)
MEpois1 <- ME(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson", listw=spdep::nb2listw(ncCR85_nb, style="B"), alpha=0.2, verbose=TRUE)
#> eV[,1], I: 0.1327384 ZI: NA, pr(ZI): 0.03
#> eV[,8], I: 0.06936385 ZI: NA, pr(ZI): 0.12
#> eV[,4], I: 0.03584503 ZI: NA, pr(ZI): 0.3
MEpois1
#>   Eigenvector ZI pr(ZI)
#> 0          NA NA   0.01
#> 1           1 NA   0.03
#> 2           8 NA   0.12
#> 3           4 NA   0.30
glmME <- glm(SID74 ~ 1 + fitted(MEpois1), data=nc.sids, offset=log(BIR74),
 family="poisson")
anova(glmME, test="Chisq")
#> Analysis of Deviance Table
#> 
#> Model: poisson, link: log
#> 
#> Response: SID74
#> 
#> Terms added sequentially (first to last)
#> 
#> 
#>                 Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#> NULL                               99     203.34              
#> fitted(MEpois1)  3   32.499        96     170.84 4.108e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(glmbase, glmME, test="Chisq")
#> Analysis of Deviance Table
#> 
#> Model 1: SID74 ~ 1
#> Model 2: SID74 ~ 1 + fitted(MEpois1)
#>   Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#> 1        99     203.34                          
#> 2        96     170.84  3   32.499 4.108e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
