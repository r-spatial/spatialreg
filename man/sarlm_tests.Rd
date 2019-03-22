% Copyright 2002-9 by Roger S. Bivand
\name{LR.sarlm}
\alias{LR.sarlm}
\alias{LR1.sarlm}
\alias{Wald1.sarlm}
\alias{Hausman.test}
\alias{Hausman.test.sarlm}
\alias{logLik.sarlm}
\alias{anova.sarlm}
\alias{bptest.sarlm}
\alias{impacts.sarlm}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Likelihood ratio test}
\description{
  The \code{LR.sarlm()} function provides a likelihood ratio test for objects for which a \code{logLik()} function exists for their class, or for objects of class \code{logLik}. \code{LR1.sarlm()} and \code{Wald1.sarlm()} are used internally in \code{summary.sarlm()}, but may be accessed directly; they report the values respectively of LR and Wald tests for the absence of spatial dependence in spatial lag or error models. The spatial Hausman test is available for models fitted with \code{errorsarlm} and \code{GMerrorsar}.
}
\usage{
LR.sarlm(x, y)
\method{logLik}{sarlm}(object, ...)
LR1.sarlm(object)
Wald1.sarlm(object)
\method{Hausman.test}{sarlm}(object, ..., tol=NULL)
\method{anova}{sarlm}(object, ...)
bptest.sarlm(object, varformula=NULL, studentize = TRUE, data=list())
\method{impacts}{sarlm}(obj, \dots, tr, R = NULL, listw = NULL, evalues=NULL,
 useHESS = NULL, tol = 1e-06, empirical = FALSE, Q=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{y}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{object, obj}{a \code{sarlm} object}
  \item{\dots}{further arguments passed to or from other methods}
%  \item{tol}{\code{tol} argument passed to \code{solve}, default NULL}
  \item{varformula}{a formula describing only the potential explanatory variables for the variance (no dependent variable needed). By default the same explanatory variables are taken as in the main regression model}
  \item{studentize}{logical. If set to \code{TRUE} Koenker's studentized
   version of the test statistic will be used.}
  \item{data}{an optional data frame containing the variables in the varformula}  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{listw}{If \code{tr} is not given, a spatial weights object as created by \code{nb2listw}; they must be the same spatial weights as were used in fitting the spatial regression, but do not have to be row-standardised}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects; the objects are used for convenience but are not output by an MCMC process}
  \item{useHESS}{Use the Hessian approximation (if available) even if the asymptotic coefficient covariance matrix is available; used for comparing methods}
  \item{tol}{Argument passed to \code{mvrnorm} and \code{solve}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
}
\value{
  The tests return objects of class \code{htest} with:
  \item{statistic}{value of statistic}
  \item{parameter}{degrees of freedom}
  \item{p.value}{Probability value}
  \item{estimate}{varies with test}
  \item{method}{description of test method}

  \code{logLik.sarlm()} returns an object of class \code{logLik}
  \code{LR1.sarlm}, \code{Hausman.sarlm} and \code{Wald1.sarlm} returm objects of class \code{htest}
}

\note{The numbers of degrees of freedom returned by \code{logLik.sarlm()} include nuisance parameters, that is the number of regression coefficients, plus sigma, plus spatial parameter esitmate(s).}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 61--63; Pace RK and LeSage J (2008) A spatial Hausman test. \emph{Economics Letters} 101, 282--284.
T.S. Breusch & A.R. Pagan (1979),
A Simple Test for Heteroscedasticity and Random Coefficient Variation.
\emph{Econometrica} \bold{47}, 1287--1294

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica.

L. Anselin (1988) \emph{Spatial econometrics: methods and models.}
Dordrecht: Kluwer, pp. 121--122.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, \code{bptest}: Torsten Hothorn and Achim Zeileis, modified by Roger Bivand}

\seealso{\code{\link{logLik.lm}}, \code{\link{anova.sarlm}}}

\examples{
require("sf", quietly=TRUE)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
lm.mod <- lm(CRIME ~ HOVAL + INC, data=columbus)
lag <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb))
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb), type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, spdep::nb2listw(col.gal.nb))
Hausman.test(error)
LR.sarlm(mixed, error)
anova(lag, lm.mod)
anova(lag, error, mixed)
AIC(lag, error, mixed)
bptest.sarlm(error)
bptest.sarlm(error, studentize=FALSE)
\dontrun{
lm.target <- lm(error$tary ~ error$tarX - 1)
if (require(lmtest, quietly=TRUE) && require(sandwich, quietly=TRUE)) {
  print(coeftest(lm.target, vcov=vcovHC(lm.target, type="HC0"), df=Inf))
}
}
}
\keyword{spatial}
