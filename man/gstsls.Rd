\name{gstsls}
\alias{gstsls}
\alias{impacts.Gmsar}
\title{Spatial simultaneous autoregressive SAC model estimation by GMM}
\description{
  An implementation of Kelejian and Prucha's generalised moments estimator for the autoregressive parameter in a spatial model
  with a spatially lagged dependent variable.
}
\usage{
gstsls(formula, data = list(), listw, listw2 = NULL, na.action = na.fail, 
    zero.policy = attr(listw, "zero.policy"), pars=NULL, scaleU=FALSE, control = list(), 
    verbose=NULL, method="nlminb", robust=FALSE, legacy=FALSE, W2X=TRUE) 
\method{impacts}{Gmsar}(obj, \dots, n = NULL, tr = NULL, R = NULL,
 listw = NULL, evalues=NULL, tol = 1e-06, empirical = FALSE, Q=NULL)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{listw2}{a \code{listw} object created for example by \code{nb2listw}, if not given, set to the same spatial weights as the listw argument}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{GMerrorsar()} to terminate with an error}
  \item{pars}{starting values for \eqn{\lambda}{lambda} and \eqn{\sigma^2}{sigma squared} for GMM optimisation, 
  if missing (default), approximated from initial 2sls model as the autocorrelation coefficient corrected for weights style 
  and model sigma squared}
  \item{scaleU}{Default FALSE: scale the OLS residuals before computing the moment matrices; only used if the \code{pars} argument is missing}
  \item{control}{A list of control parameters. See details in \link{optim} or \link{nlminb}}
  \item{verbose}{default NULL, use global option value; if TRUE, reports function values during optimization.}
  \item{method}{default \link{nlminb}, or optionally a method passed to \link{optim} to use an alternative optimizer}
  \item{robust}{see \code{stsls}}
  \item{legacy}{see \code{stsls}}
  \item{W2X}{see \code{stsls}}
  \item{obj}{A spatial regression object created by \code{lagsarlm}, \code{lagmess} or by \code{lmSLX}; in \code{HPDinterval.LagImpact}, a LagImpact object}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects; the objects are used for convenience but are not output by an MCMC process}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
  \item{n}{defaults to \code{length(obj$residuals)}; in the method for \code{Gmsar} objects it may be used in panel settings to compute the impacts for cross-sectional weights only, suggested by Angela Parenti}

}
\details{
  When the control list is set with care, the function will converge to values close to the ML estimator 
  without requiring computation of the Jacobian, the most resource-intensive part of ML estimation. 
}
\value{
  A list object of class \code{Gmsar}
  \item{lambda}{simultaneous autoregressive error coefficient}
  \item{coefficients}{GMM coefficient estimates (including the spatial autocorrelation coefficient)}
  \item{rest.se}{GMM coefficient standard errors}
  \item{s2}{GMM residual variance}
  \item{SSE}{sum of squared GMM errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{NULL}
  \item{call}{the call used to create this object}
  \item{residuals}{GMM residuals}
  \item{lm.target}{NULL}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{formula}{model formula}
  \item{aliased}{NULL}
  \item{zero.policy}{zero.policy for this model}
  \item{LL}{NULL}
  \item{vv}{list of internal bigG and litg components for testing optimisation surface}
  \item{optres}{object returned by optimizer}
  \item{pars}{start parameter values for optimisation}
  \item{Hcov}{NULL}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
}
\references{Kelejian, H. H., and Prucha, I. R., 1999. A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model. International Economic Review, 40, pp. 509--533; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \doi{10.18637/jss.v063.i18}.
}
\author{Gianfranco Piras and Roger Bivand}

\seealso{\code{\link[stats]{optim}}, \code{\link[stats]{nlminb}}, \code{\link{GMerrorsar}}, \code{\link{GMargminImage}}}
\examples{
#require("spdep", quietly=TRUE) 
data(oldcol, package="spdep")
COL.errW.GM <- gstsls(CRIME ~ INC + HOVAL, data=COL.OLD, spdep::nb2listw(COL.nb, style="W"))
summary(COL.errW.GM)
aa <- GMargminImage(COL.errW.GM)
levs <- quantile(aa$z, seq(0, 1, 1/12))
image(aa, breaks=levs, xlab="lambda", ylab="s2")
points(COL.errW.GM$lambda, COL.errW.GM$s2, pch=3, lwd=2)
contour(aa, levels=signif(levs, 4), add=TRUE)
COL.errW.GM <- gstsls(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="W"), scaleU=TRUE)
summary(COL.errW.GM)
listw <- spdep::nb2listw(COL.nb)
W <- as(listw, "CsparseMatrix")
trMat <- trW(W, type="mult")
impacts(COL.errW.GM, tr=trMat)
}
\keyword{spatial}
