<!-- badges: start -->
[![R-CMD-check](https://github.com/r-spatial/spatialreg/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/r-spatial/spatialreg/actions/workflows/check-standard.yaml)
[![CRAN](http://www.r-pkg.org/badges/version/spatialreg)](https://cran.r-project.org/package=spatialreg)
<!-- badges: end -->

# spatialreg

### spatialreg: spatial models estimation and testing

A collection of all the estimation functions for spatial cross-sectional models (on lattice/areal data using spatial weights matrices) contained up to now in **spdep**. These model fitting functions include maximum likelihood methods for cross-sectional models proposed by Cliff and Ord (1973, ISBN: 0850860369) and (1981, ISBN: 0850860814), fitting methods initially described by Ord (1975) https://doi.org/10.1080/01621459.1975.10480272. The models are further described by Anselin (1988) https://doi.org/10.1007/978-94-015-7799-1. Spatial two stage least squares and spatial general method of moment models initially proposed by Kelejian and Prucha (1998) https://doi.org/10.1023/A:1007707430416 and (1999) https://doi.org/10.1111/1468-2354.00027 are provided. Impact methods and MCMC fitting methods proposed by LeSage and Pace (2009) https://doi.org/10.1201/9781420064254 are implemented for the family of cross-sectional spatial regression models. Methods for fitting the log determinant term in maximum likelihood and MCMC fitting are compared by Bivand et al. (2013) https://doi.org/10.1111/gean.12008, and model fitting methods by Bivand and Piras (2015) https://doi.org/10.18637/jss.v063.i18; both of these articles include extensive lists of references. A recent review is provided by Bivand, Millo and Piras (2021) https://doi.org/10.3390/math9111276. **spatialreg** >= 1.1-1 corresponds to **spdep** = 1.1-1, in which the model fitting functions are deprecated and pass through to **spatialreg**, but will mask those in **spatialreg**. From versions 1.2-1, the functions have been made defunct in **spdep**. From version 1.3-6, add Anselin-Kelejian (1997) test to `stsls` for residual spatial autocorrelation https://doi.org/10.1177/016001769702000109.

Default branch now `main`.
