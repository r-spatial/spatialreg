% Copyright 2001-3 by Roger S. Bivand
\name{as.spam.listw}
\alias{as.spam.listw}
\alias{listw2U_spam}
\alias{as_dgRMatrix_listw}
\alias{as_dsTMatrix_listw}
\alias{as_dsCMatrix_I}
\alias{as_dsCMatrix_IrW}
\alias{listw2U_Matrix}
\alias{Jacobian_W}
\alias{coerce,listw,CsparseMatrix-method}
\alias{coerce,listw,RsparseMatrix-method}
\alias{coerce,listw,symmetricMatrix-method}
\alias{powerWeights}
\title{Spatial neighbour sparse representation}
\description{
Interface between Matrix class objects and weights lists. The \code{as.spam.listw} method converts a \code{"listw"} object to a sparse matrix as defined in the \pkg{spam} package.
}
\usage{
as.spam.listw(listw)
listw2U_spam(lw)
listw2U_Matrix(lw)
as_dgRMatrix_listw(listw)
as_dsTMatrix_listw(listw)
as_dsCMatrix_I(n)
as_dsCMatrix_IrW(W, rho)
Jacobian_W(W, rho)
powerWeights(W, rho, order=250, X, tol=.Machine$double.eps^(3/5))
}
\arguments{
  \item{listw, lw}{a \code{listw} object from for example \code{nb2listw}}
  \item{W}{a \code{dsTMatrix} object created using \code{as_dsTMatrix_listw} from a symmetric \code{listw} object}
  \item{rho}{spatial regression coefficient}
  \item{n}{length of diagonal for identity matrix}
  \item{order}{Power series maximum limit}
  \item{X}{A numerical matrix}
  \item{tol}{Tolerance for convergence of power series}
}



\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[spdep]{nb2listw}}}

\examples{
\dontrun{
require(sf, quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#require(spdep, quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- spdep::nb2listw(col.gal.nb)
if (require("spam", quietly=TRUE)) {
  col.sp <- as.spam.listw(col.listw)
  str(col.sp)
}
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
nyadjlw <- spdep::mat2listw(nyadjmat)
listw_NY <- spdep::nb2listw(nyadjlw$neighbours, style="B")
W_C <- as(listw_NY, "CsparseMatrix")
W_R <- as(listw_NY, "RsparseMatrix")
W_S <- as(listw_NY, "symmetricMatrix")
n <- nrow(W_S)
I <- Diagonal(n)
rho <- 0.1
c(determinant(I - rho * W_S, logarithm=TRUE)$modulus)
sum(log(1 - rho * eigenw(listw_NY)))
nW <- - W_S
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho))$modulus))
}
nb7rt <- spdep::cell2nb(7, 7, torus=TRUE)
x <- matrix(sample(rnorm(500*length(nb7rt))), nrow=length(nb7rt))
lw <- spdep::nb2listw(nb7rt)
if (FALSE) {
# Only needed in some simulation settings where the input and
# output distributions must agree in all but autocorrelation
x <- apply(x, 2, scale)
st <- apply(x, 2, function(x) shapiro.test(x)$p.value)
x <- x[, (st > 0.2 & st < 0.8)]
x <- apply(x, 2, function(v) residuals(spautolm(v ~ 1, listw=lw,
 method="eigen", control=list(pre_eig=e, fdHess=FALSE))))
x <- apply(x, 2, scale)
}
W <- as(lw, "CsparseMatrix")
system.time(e <- spdep::invIrM(nb7rt, rho=0.98, method="solve", feasible=NULL) \%*\% x)
system.time(ee <- powerWeights(W, rho=0.98, X=x))
str(attr(ee, "internal"))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
\dontrun{
system.time(ee <- powerWeights(W, rho=0.9, X=x))
system.time(ee <- powerWeights(W, rho=0.98, order=1000, X=x))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
nb60rt <- spdep::cell2nb(60, 60, torus=TRUE)
W <- as(spdep::nb2listw(nb60rt), "CsparseMatrix")
set.seed(1)
x <- matrix(rnorm(dim(W)[1]), ncol=1)
system.time(ee <- powerWeights(W, rho=0.3, X=x))
str(as(ee, "matrix"))
obj <- errorsarlm(as(ee, "matrix")[,1] ~ 1, listw=spdep::nb2listw(nb60rt), method="Matrix")
coefficients(obj)
}
}
\keyword{spatial}
