# Create symmetric similar weights lists

From Ord's 1975 paper, it is known that the Jacobian for SAR models may
be found by "symmetrizing" by similarity (the eigenvalues of similar
matrices are identical, so the Jacobian is too). This applies only to
styles "W" and "S" with underlying symmetric binary neighbour relations
or symmetric general neighbour relations (so no k-nearest neighbour
relations). The function is invoked automatically within the SAR fitting
functions, to call `eigen` on a symmetric matrix for the default eigen
method, or to make it possible to use the Matrix method on weights that
can be "symmetrized" in this way.

## Usage

``` r
similar.listw(listw)
```

## Arguments

- listw:

  a `listw` object created for example by
  [`spdep::nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html)

## Value

a `listw` object

## References

Ord, J. K. 1975 Estimation methods for models of spatial interaction,
*Journal of the American Statistical Association*, 70, 120-126

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`lagsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md),
[`errorsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.md)

## Examples

``` r
#require("spdep", quietly=TRUE)
data(oldcol, package="spdep")
COL.W <- spdep::nb2listw(COL.nb, style="W")
COL.S <- spdep::nb2listw(COL.nb, style="S")
sum(log(1 - 0.5 * eigenw(COL.W)))
#> [1] -1.62766
sum(log(1 - 0.5 * eigenw(similar.listw(COL.W))))
#> [1] -1.62766
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.W)), "CsparseMatrix")
I <- as_dsCMatrix_I(dim(W_J)[1])
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
#> [1] -1.62766
sum(log(1 - 0.5 * eigenw(COL.S)))
#> [1] -1.602757
sum(log(1 - 0.5 * eigenw(similar.listw(COL.S))))
#> [1] -1.602757
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.S)), "CsparseMatrix")
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
#> [1] -1.602757
```
