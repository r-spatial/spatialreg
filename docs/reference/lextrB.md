# Find extreme eigenvalues of binary symmetric spatial weights

The functions find extreme eigenvalues of binary symmetric spatial
weights, when these form planar graphs; general weights are not
permiited. `l_max` finds the largest eigenvalue using Rayleigh quotient
methods of any “listw” object. `lextrB` first calls `l_max`, and uses
its output to find the smallest eigenvalue in addition for binary
symmetric spatial weights. `lextrW` extends these to find the smallest
eigenvalue for intrinsically symmetric row-standardized binary weights
matrices (transformed to symmetric through similarity internally).
`lextrS` does the same for variance-stabilized (“S” style) intrinsically
symmetric binary weights matrices (transformed to symmetric through
similarity internally).

## Usage

``` r
lextrB(lw, zero.policy = TRUE, control = list())
lextrW(lw, zero.policy=TRUE, control=list())
lextrS(lw, zero.policy=TRUE, control=list())
l_max(lw, zero.policy=TRUE, control=list())
```

## Arguments

- lw:

  a binary symmetric `listw` object from, for example, `nb2listw` with
  style “B” for `lextrB`, style “W” for `lextrW` and style “S” for
  `lextrS`; for `l_max`, the object may be asymmetric and does not have
  to be binary

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE assign NA

- control:

  a list of control arguments

## Control arguments

- trace:

  report values in while loops, default NULL assuming FALSE; logical

- tol:

  tolerance for breaking while loops, default
  `.Machine$double.eps^(1/2)`; numeric

- maxiter:

  maximum number of iterations in while loops, default
  `6 * (length(lw$neighbours) - 2`; integer

- useC:

  use C code, default TRUE, logical (not in `l_max`)

## Value

The functions return approximations to the extreme eigenvalues with the
eigenvectors returned as attributes of this object.

## References

Griffith, D. A. (2004). Extreme eigenfunctions of adjacency matrices for
planar graphs employed in spatial analyses. *Linear Algebra and its
Applications*, 388:201–219.

## Author

Roger Bivand, Yongwan Chun, Daniel Griffith

## Note

It may be necessary to modify control arguments if warnings about lack
of convergence are seen.

## Examples

``` r
data(boston, package="spData")
#require(spdep, quietly=TRUE)
ab.listb <- spdep::nb2listw(boston.soi, style="B")
er <- range(eigenw(ab.listb))
er
#> [1] -3.039465  5.306204
res_1 <- lextrB(ab.listb)
c(res_1)
#>  lambda_n  lambda_1 
#> -3.039374  5.306203 
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
B <- as(ab.listb, "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
#> [1] -3.039465
if (run) {
eigs(B, k=1, which="LR")$values
}
#> [1] 5.306204
k5 <- spdep::knn2nb(spdep::knearneigh(boston.utm, k=5))
c(l_max(spdep::nb2listw(k5, style="B")))
#> [1] 5
max(Re(eigenw(spdep::nb2listw(k5, style="B"))))
#> [1] 5
c(l_max(spdep::nb2listw(k5, style="C")))
#> [1] 1
max(Re(eigenw(spdep::nb2listw(k5, style="C"))))
#> [1] 1
ab.listw <- spdep::nb2listw(boston.soi, style="W")
er <- range(eigenw(similar.listw(ab.listw)))
er
#> [1] -0.9708644  1.0000000
res_1 <- lextrW(ab.listw)
c(res_1)
#>   lambda_n   lambda_1 
#> -0.9708644  0.9999991 
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
#> [1] -0.9708644
if (run) {
eigs(B, k=1, which="LR")$values
}
#> [1] 1
# \dontrun{
ab.listw <- spdep::nb2listw(boston.soi, style="S")
er <- range(eigenw(similar.listw(ab.listw)))
er
#> [1] -0.723495  1.110373
res_1 <- lextrS(ab.listw)
c(res_1)
#>   lambda_n   lambda_1 
#> -0.7230376  1.1103694 
# }
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
#> [1] -0.723495
if (run) {
eigs(B, k=1, which="LR")$values
}
#> [1] 1.110373
```
