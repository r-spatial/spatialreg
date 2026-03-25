# Spatial neighbour sparse representation

Interface between Matrix class objects and weights lists. The
`as.spam.listw` method converts a `"listw"` object to a sparse matrix as
defined in the spam package.

## Usage

``` r
as.spam.listw(listw)
listw2U_spam(lw)
listw2U_Matrix(lw)
as_dgRMatrix_listw(listw)
as_dsTMatrix_listw(listw)
as_dsCMatrix_I(n)
as_dsCMatrix_IrW(W, rho)
Jacobian_W(W, rho)
powerWeights(W, rho, order=250, X, tol=.Machine$double.eps^(3/5))
```

## Arguments

- listw, lw:

  a `listw` object from for example `nb2listw`

- W:

  a `dsTMatrix` object created using `as_dsTMatrix_listw` from a
  symmetric `listw` object

- rho:

  spatial regression coefficient

- n:

  length of diagonal for identity matrix

- order:

  Power series maximum limit

- X:

  A numerical matrix

- tol:

  Tolerance for convergence of power series

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html)

## Examples

``` r
# \dontrun{
require(sf, quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
#require(spdep, quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- spdep::nb2listw(col.gal.nb)
if (require("spam", quietly=TRUE)) {
  col.sp <- as.spam.listw(col.listw)
  str(col.sp)
}
#> Spam version 2.11-3 (2026-01-05) is loaded.
#> Type 'help( Spam)' or 'demo( spam)' for a short introduction 
#> and overview of this package.
#> Help for individual functions is also obtained by adding the
#> suffix '.spam' to the function name, e.g. 'help( chol.spam)'.
#> 
#> Attaching package: ‘spam’
#> The following object is masked from ‘package:Matrix’:
#> 
#>     det
#> The following objects are masked from ‘package:base’:
#> 
#>     backsolve, forwardsolve
#> Formal class 'spam' [package "spam"] with 4 slots
#>   ..@ entries    : num [1:230] 0.5 0.5 0.333 0.333 0.333 ...
#>   ..@ colindices : int [1:230] 2 3 1 3 4 1 2 4 5 2 ...
#>   ..@ rowpointers: int [1:50] 1 3 6 10 14 21 23 27 33 41 ...
#>   ..@ dimension  : int [1:2] 49 49
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
nyadjlw <- spdep::mat2listw(nyadjmat)
#> Warning: style is M (missing); style should be set to a valid value
listw_NY <- spdep::nb2listw(nyadjlw$neighbours, style="B")
W_C <- as(listw_NY, "CsparseMatrix")
W_R <- as(listw_NY, "RsparseMatrix")
W_S <- as(listw_NY, "symmetricMatrix")
n <- nrow(W_S)
I <- Diagonal(n)
rho <- 0.1
c(determinant(I - rho * W_S, logarithm=TRUE)$modulus)
#> [1] -9.587255
sum(log(1 - rho * eigenw(listw_NY)))
#> [1] -9.587255
nW <- - W_S
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho))$modulus))
#> [1] 99.8069
# }
nb7rt <- spdep::cell2nb(7, 7, torus=TRUE)
x <- matrix(sample(rnorm(500*length(nb7rt))), nrow=length(nb7rt))
lw <- spdep::nb2listw(nb7rt)
if (FALSE) {
# Only needed in some simulation settings where the input and
# output distributions must agree in all but autocorrelation
e <- eigenw(lw)
x <- apply(x, 2, scale)
st <- apply(x, 2, function(x) shapiro.test(x)$p.value)
x <- x[, (st > 0.2 & st < 0.8)]
x <- apply(x, 2, function(v) residuals(spautolm(v ~ 1, listw=lw,
 method="eigen", control=list(pre_eig=e, fdHess=FALSE))))
x <- apply(x, 2, scale)
}
W <- as(lw, "CsparseMatrix")
system.time(e <- invIrM(nb7rt, rho=0.98, method="solve", feasible=NULL) %*% x)
#>    user  system elapsed 
#>   0.002   0.000   0.002 
system.time(ee <- powerWeights(W, rho=0.98, X=x))
#> Warning: not converged within order iterations
#>    user  system elapsed 
#>   0.197   0.011   0.209 
str(attr(ee, "internal"))
#> List of 5
#>  $ series: num [1:250] 0.287 0.234 0.201 0.178 0.16 ...
#>  $ order : num 250
#>  $ tol   : num 4.05e-10
#>  $ iter  : num 250
#>  $ conv  : logi FALSE
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
#> [1] "Mean relative difference: 0.0060747"
# \dontrun{
system.time(ee <- powerWeights(W, rho=0.9, X=x))
#>    user  system elapsed 
#>   0.145   0.014   0.164 
system.time(ee <- powerWeights(W, rho=0.98, order=1000, X=x))
#>    user  system elapsed 
#>   0.745   0.011   0.772 
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
#> [1] TRUE
nb60rt <- spdep::cell2nb(60, 60, torus=TRUE)
W <- as(spdep::nb2listw(nb60rt), "CsparseMatrix")
set.seed(1)
x <- matrix(rnorm(dim(W)[1]), ncol=1)
system.time(ee <- powerWeights(W, rho=0.3, X=x))
#>    user  system elapsed 
#>   0.011   0.000   0.011 
str(as(ee, "matrix"))
#>  num [1:3600, 1] -0.383 0.207 -0.731 1.552 0.32 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:3600] "1:1" "2:1" "3:1" "4:1" ...
#>   ..$ : NULL
obj <- errorsarlm(as(ee, "matrix")[,1] ~ 1, listw=spdep::nb2listw(nb60rt), method="Matrix")
coefficients(obj)
#>      lambda (Intercept) 
#>  0.30639876  0.01380415 
# }
```
