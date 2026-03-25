# Approximate profile-likelihood estimator (APLE)

The Approximate profile-likelihood estimator (APLE) of the simultaneous
autoregressive model's spatial dependence parameter was introduced in Li
et al. (2007). It employs a correction term using the eigenvalues of the
spatial weights matrix, and consequently should not be used for large
numbers of observations. It also requires that the variable has a mean
of zero, and it is assumed that it has been detrended. The spatial
weights object is assumed to be row-standardised, that is using default
`style="W"` in `nb2listw`.

## Usage

``` r
aple(x, listw, override_similarity_check=FALSE, useTrace=TRUE)
```

## Arguments

- x:

  a zero-mean detrended continuous variable

- listw:

  a `listw` object from for example
  [`spdep::nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html)

- override_similarity_check:

  default FALSE, if TRUE - typically for row-standardised weights with
  asymmetric underlying general weights - similarity is not checked

- useTrace:

  default TRUE, use trace of sparse matrix `W %*% W` (Li et al. (2010)),
  if FALSE, use crossproduct of eigenvalues of `W` as in Li et al.
  (2007)

## Details

This implementation has been checked with Hongfei Li's own
implementation using her data; her help was very valuable.

## Value

A scalar APLE value.

## References

Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I:
testing for spatial dependence based on the spatial autoregressive
model. Geographical Analysis 39, 357-375; Li, H, Calder, C. A. and
Cressie N. A. C. (2012) One-step estimation of spatial dependence
parameters: Properties and extensions of the APLE statistic, Journal of
Multivariate Analysis 105, 68-84.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html),
[`aple.mc`](https://r-spatial.github.io/spatialreg/reference/aple.mc.md),
[`aple.plot`](https://r-spatial.github.io/spatialreg/reference/aple.plot.md)

## Examples

``` r
wheat <- st_read(system.file("shapes/wheat.gpkg", package="spData")[1], quiet=TRUE)
library(spdep)
#> 
#> Attaching package: ‘spdep’
#> The following objects are masked from ‘package:spatialreg’:
#> 
#>     get.ClusterOption, get.VerboseOption, get.ZeroPolicyOption,
#>     get.coresOption, get.mcOption, set.ClusterOption,
#>     set.VerboseOption, set.ZeroPolicyOption, set.coresOption,
#>     set.mcOption
nbr1 <- spdep::poly2nb(wheat, queen=FALSE)
nbrl <- spdep::nblag(nbr1, 2)
#> Warning: lag 2 neighbour object has 2 sub-graphs
nbr12 <- spdep::nblag_cumul(nbrl)
cms0 <- with(as.data.frame(wheat), tapply(yield, c, median))
cms1 <- c(model.matrix(~ factor(c) -1, data=wheat) %*% cms0)
wheat$yield_detrend <- wheat$yield - cms1
isTRUE(all.equal(c(with(as.data.frame(wheat),
 tapply(yield_detrend, c, median))), rep(0.0, 25),
 check.attributes=FALSE))
#> [1] TRUE
spdep::moran.test(wheat$yield_detrend, spdep::nb2listw(nbr12, style="W"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  wheat$yield_detrend  
#> weights: spdep::nb2listw(nbr12, style = "W")    
#> 
#> Moran I statistic standard deviate = 10.305, p-value < 2.2e-16
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      0.1935469610     -0.0020040080      0.0003600869 
#> 
aple(as.vector(scale(wheat$yield_detrend, scale=FALSE)), spdep::nb2listw(nbr12, style="W"))
#> [1] 0.6601805
# \dontrun{
errorsarlm(yield_detrend ~ 1, wheat, spdep::nb2listw(nbr12, style="W"))
#> 
#> Call:
#> errorsarlm(formula = yield_detrend ~ 1, data = wheat, listw = spdep::nb2listw(nbr12, 
#>     style = "W"))
#> Type: error 
#> 
#> Coefficients:
#>      lambda (Intercept) 
#>  0.60189686 -0.00251772 
#> 
#> Log likelihood: -192.9519 
# }
```
