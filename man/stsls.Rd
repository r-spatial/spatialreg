\name{stsls}
\alias{stsls}
\alias{print.stsls}
\alias{print.summary.stsls}
\alias{summary.stsls}
\alias{residuals.stsls}
\alias{coef.stsls}
\alias{deviance.stsls}
\alias{impacts.stsls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized spatial two stage least squares}
\description{
  The function fits a spatial lag model by two stage least squares, with the option of adjusting the results for heteroskedasticity.
}
\usage{
stsls(formula, data = list(), listw, zero.policy = NULL,
 na.action = na.fail, robust = FALSE, HC=NULL, legacy=FALSE, W2X = TRUE)
\method{impacts}{stsls}(obj, \dots, tr, R = NULL, listw = NULL, evalues=NULL,
 tol = 1e-06, empirical = FALSE, Q=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{robust}{default FALSE, if TRUE, apply a heteroskedasticity correction to the coefficients covariances}
  \item{HC}{default NULL, if \code{robust} is TRUE, assigned \dQuote{HC0}, may take values \dQuote{HC0} or \dQuote{HC1} for White estimates or MacKinnon-White estimates respectively}
 \item{legacy}{the argument chooses between two implementations of the robustness correction: default FALSE - use the estimate of Omega only in the White consistent estimator of the variance-covariance matrix, if TRUE, use the original implementation which runs a GLS using the estimate of Omega, and yields different coefficient estimates as well - see example below}
 \item{W2X}{default TRUE, if FALSE only WX are used as instruments in the spatial two stage least squares; until release 0.4-60, only WX were used - see example below }
  \item{obj}{A spatial regression object created by \code{lagsarlm}, \code{lagmess} or by \code{lmSLX}; in \code{HPDinterval.lagImpact}, a lagImpact object}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects; the objects are used for convenience but are not output by an MCMC process}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
}
\details{
  The fitting implementation fits a spatial lag model:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

by using spatially lagged X variables as instruments for the spatially lagged dependent variable.
}
\value{
  an object of class "stsls" containing:
  \item{coefficients}{coefficient estimates}
  \item{var}{coefficient covariance matrix}
  \item{sse}{sum of squared errors}
  \item{residuals}{model residuals}
  \item{df}{degrees of freedom}
}
\references{Kelejian, H.H. and I.R. Prucha (1998). A generalized spatial two
stage least squares procedure for estimating a spatial autoregressive
model with autoregressive disturbances. \emph{Journal of Real Estate
Finance and Economics} 17, 99-121.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.
}
\author{Luc Anselin, Gianfranco Piras and Roger Bivand}

\seealso{\code{\link{lagsarlm}}}
\examples{
data(oldcol, package="spdep")
#require(spdep, quietly=TRUE)
lw <- spdep::nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
summary(COL.lag.eig, correlation=TRUE)
COL.lag.stsls <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
(x <- summary(COL.lag.stsls, correlation=TRUE))
coef(x)
W <- as(lw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
loobj1 <- impacts(COL.lag.stsls, R=200, tr=trMatc)
summary(loobj1, zstats=TRUE, short=TRUE)
ev <- eigenw(lw)
loobj2 <- impacts(COL.lag.stsls, R=200, evalues=ev)
summary(loobj2, zstats=TRUE, short=TRUE)
require(coda)
HPDinterval(loobj1)
COL.lag.stslsW <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw, W2X=FALSE)
summary(COL.lag.stslsW, correlation=TRUE)
COL.lag.stslsR <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
robust=TRUE, W2X=FALSE)
summary(COL.lag.stslsR, correlation=TRUE)
COL.lag.stslsRl <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
robust=TRUE, legacy=TRUE, W2X=FALSE)
summary(COL.lag.stslsRl, correlation=TRUE)
data(boston, package="spData")
gp2a <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, spdep::nb2listw(boston.soi))
summary(gp2a)
}
\keyword{spatial}
