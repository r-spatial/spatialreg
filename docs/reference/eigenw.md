# Spatial weights matrix eigenvalues

The `eigenw` function returns a numeric vector of eigenvalues of the
weights matrix generated from the spatial weights object `listw`. The
eigenvalues are used to speed the computation of the Jacobian in spatial
model estimation:

\$\$\log(\det\[I - \rho W\]) = \sum\_{i=1}^{n}\log(1 - \rho
\lambda_i)\$\$

where \\W\\ is the n by n spatial weights matrix, and \\\lambda_i\\ are
the eigenvalues of \\W\\.

## Usage

``` r
eigenw(listw, quiet=NULL)
griffith_sone(P, Q, type="rook")
subgraph_eigenw(nb, glist=NULL, style="W", zero.policy=NULL, quiet=NULL)
```

## Arguments

- listw:

  a `listw` object created for example by `nb2listw`

- quiet:

  default NULL, use global !verbose option value; set to FALSE for short
  summary

- P:

  number of columns in the grid (number of units in a horizontal axis
  direction)

- Q:

  number of rows in the grid (number of units in a vertical axis
  direction.)

- type:

  “rook” or “queen”

- nb:

  an object of class `nb`

- glist:

  list of general weights corresponding to neighbours

- style:

  `style` can take values “W”, “B”, “C”, “U”, “minmax” and “S”

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

## Details

The `griffith_sone` function function may be used, following Ord and
Gasim (for references see Griffith and Sone (1995)), to calculate
analytical eigenvalues for binary rook or queen contiguous neighbours
where the data are arranged as a regular P times Q grid. The
`subgraph_eigenw` function may be used when there are multiple graph
components, of which the largest may be handled as a dense matrix. Here
the eigenvalues are computed for each subgraph in turn, and catenated to
reconstruct the complete set. The functions may be used to provide
pre-computed eigenvalues for spatial regression functions.

## Value

a numeric or complex vector of eigenvalues of the weights matrix
generated from the spatial weights object.

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 155; Ord, J.
K. 1975 Estimation methods for models of spatial interaction, Journal of
the American Statistical Association, 70, 120-126.; Griffith, D. A. and
Sone, A. (1995). Trade-offs associated with normalizing constant
computational simplifications for estimating spatial statistical models.
Journal of Statistical Computation and Simulation, 51, 165-183.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`eigen`](https://rdrr.io/r/base/eigen.html)

## Examples

``` r
#require(spdep)
data(oldcol, package="spdep")
W.eig <- eigenw(spdep::nb2listw(COL.nb, style="W"))
1/range(W.eig)
#> [1] -1.536177  1.000000
S.eig <- eigenw(spdep::nb2listw(COL.nb, style="S"))
1/range(S.eig)
#> [1] -1.7364189  0.8918997
B.eig <- eigenw(spdep::nb2listw(COL.nb, style="B"))
1/range(B.eig)
#> [1] -0.3229290  0.1692726
# cases for intrinsically asymmetric weights
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k3 <- spdep::knn2nb(spdep::knearneigh(crds, k=3))
#> Warning: neighbour object has 2 sub-graphs
spdep::is.symmetric.nb(k3)
#> [1] FALSE
k3eig <- eigenw(spdep::nb2listw(k3, style="W"))
is.complex(k3eig)
#> [1] TRUE
rho <- 0.5
Jc <- sum(log(1 - rho * k3eig))
# complex eigenvalue Jacobian
Jc
#> [1] -1.749705+0i
# subgraphs
nc <- attr(k3, "ncomp")
if (is.null(nc)) nc <- spdep::n.comp.nb(k3)
nc$nc
#> [1] 2
table(nc$comp.id)
#> 
#>  1  2 
#> 43  6 
k3eigSG <- subgraph_eigenw(k3, style="W")
all.equal(sort(k3eig), k3eigSG)
#> [1] TRUE
W <- as(spdep::nb2listw(k3, style="W"), "CsparseMatrix")
I <- diag(length(k3))
Jl <- sum(log(abs(diag(slot(lu(I - rho * W), "U")))))
# LU Jacobian equals complex eigenvalue Jacobian
Jl
#> [1] -1.749705
all.equal(Re(Jc), Jl)
#> [1] TRUE
# wrong value if only real part used
Jr <- sum(log(1 - rho * Re(k3eig)))
Jr
#> [1] -1.762734
all.equal(Jr, Jl)
#> [1] "Mean relative difference: 0.007391147"
# construction of Jacobian from complex conjugate pairs (Jan Hauke)
Rev <- Re(k3eig)[which(Im(k3eig) == 0)]
# real eigenvalues
Cev <- k3eig[which(Im(k3eig) != 0)]
pCev <- Cev[Im(Cev) > 0]
# separate complex conjugate pairs
RpCev <- Re(pCev)
IpCev <- Im(pCev)
# reassemble Jacobian
Jc1 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2 + (rho^2)*(IpCev^2)))
all.equal(Re(Jc), Jc1)
#> [1] TRUE
# impact of omitted complex part term in real part only Jacobian
Jc2 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2))
all.equal(Jr, Jc2)
#> [1] TRUE
# trace of asymmetric (WW) and crossprod of complex eigenvalues for APLE
sum(diag(W %*% W))
#> [1] 11.55556
crossprod(k3eig)
#>                        [,1]
#> [1,] 11.55556+1.084755e-17i
# analytical regular grid eigenvalues
rg <- spdep::cell2nb(ncol=7, nrow=7, type="rook")
rg_eig <- eigenw(spdep::nb2listw(rg, style="B"))
rg_GS <- griffith_sone(P=7, Q=7, type="rook")
all.equal(rg_eig, rg_GS)
#> [1] TRUE
# \dontrun{
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
B <- as(spdep::nb2listw(rg, style="B"), "CsparseMatrix")
res1 <- eigs(B, k=1, which="LR")$values
resn <- eigs(B, k=1, which="SR")$values
print(Re(c(resn, res1)))
}
#> [1] -3.695518  3.695518
if (run) {
print(all.equal(range(Re(rg_eig)), c(resn, res1))) 
}
#> [1] TRUE
if (run) {
lw <- spdep::nb2listw(rg, style="W")
rg_eig <- eigenw(similar.listw(lw))
print(range(Re(rg_eig)))
}
#> [1] -1  1
if (run) {
W  <- as(lw, "CsparseMatrix")
print(Re(c(eigs(W, k=1, which="SR")$values, eigs(W, k=1, which="LR")$values)))
}# }
#> [1] -1  1
```
