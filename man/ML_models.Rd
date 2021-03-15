% Copyright 2002-18 by Roger S. Bivand
\name{ML_models}
\alias{lagsarlm}
\alias{errorsarlm}
\alias{sacsarlm}
\alias{summary.Sarlm}
\alias{print.Sarlm}
\alias{print.summary.Sarlm}
\alias{residuals.Sarlm}
\alias{deviance.Sarlm}
\alias{coef.Sarlm}
\alias{vcov.Sarlm}
\alias{fitted.Sarlm}

\title{Spatial simultaneous autoregressive model estimation by maximum likelihood}

\description{
  The \code{lagsarlm} function provides Maximum likelihood estimation of spatial simultaneous autoregressive lag and spatial Durbin (mixed) models of the form:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

where \eqn{\rho}{rho} is found by \code{optimize()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently (one-dimensional search using optim performs badly on some platforms). In the spatial Durbin (mixed) model, the spatially lagged independent variables are added to X. Note that interpretation of the fitted coefficients should use impact measures, because of the feedback loops induced by the data generation process for this model. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument may need be set when the weights are not row-standardised.

  Maximum likelihood estimation of spatial simultaneous autoregressive
error models of the form:

\deqn{y = X \beta + u, u = \lambda W u + \varepsilon}{y = X beta + u, u = lambda W u + e}

where \eqn{\lambda}{lambda} is found by \code{optimize()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument may need be set when the weights are not row-standardised. When \code{etype} is \dQuote{emixed}, a so-called spatial Durbin error model is fitted.

  Maximum likelihood estimation of spatial simultaneous autoregressive
\dQuote{SAC/SARAR} models of the form:

\deqn{y = \rho W1 y + X \beta + u, u = \lambda W2 u + \varepsilon}{y = rho W1 y + X beta + u, u = lambda W2 u + e}

where \eqn{\rho}{rho} and \eqn{\lambda}{lambda} are found by \code{nlminb} or \code{optim()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently.


}

\usage{
lagsarlm(formula, data = list(), listw, na.action, Durbin, type,
 method="eigen", quiet=NULL, zero.policy=NULL, interval=NULL,
 tol.solve=.Machine$double.eps, trs=NULL, control=list())
errorsarlm(formula, data=list(), listw, na.action, weights=NULL,
 Durbin, etype, method="eigen", quiet=NULL, zero.policy=NULL,
 interval = NULL, tol.solve=.Machine$double.eps, trs=NULL, control=list())
sacsarlm(formula, data = list(), listw, listw2 = NULL, na.action, Durbin, type,
 method="eigen", quiet=NULL, zero.policy=NULL, tol.solve=.Machine$double.eps,
 llprof=NULL, interval1=NULL, interval2=NULL, trs1=NULL, trs2=NULL,
 control = list())
\method{summary}{Sarlm}(object, correlation = FALSE, Nagelkerke = FALSE,
 Hausman=FALSE, adj.se=FALSE, ...)
\method{print}{Sarlm}(x, ...)
\method{print}{summary.Sarlm}(x, digits = max(5, .Options$digits - 3),
 signif.stars = FALSE, ...)
\method{residuals}{Sarlm}(object, ...)
\method{deviance}{Sarlm}(object, ...)
\method{coef}{Sarlm}(object, ...)
\method{vcov}{Sarlm}(object, ...)
\method{fitted}{Sarlm}(object, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw, listw2}{a \code{listw} object created for example by \code{nb2listw}; if \code{nb2listw} not given, set to the same spatial weights as the \code{listw} argument}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{weights}{an optional vector of weights to be used in the fitting process. Non-NULL weights can be used to indicate that different observations have different variances (with the values in weights being inversely proportional to the variances); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations (including the case that there are w_i observations equal to y_i and the data have been summarized) - \code{\link{lm}}}
  \item{Durbin}{default FALSE (spatial lag model); if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
  \item{type}{(use the \sQuote{Durbin=} argument - retained for backwards compatibility only) default "lag", may be set to "mixed"; when "mixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included; \dQuote{Durbin} may be used instead of \dQuote{mixed}}
  \item{etype}{(use the \sQuote{Durbin=} argument - retained for backwards compatibility only) default "error", may be set to "emixed" to include the spatially lagged independent variables added to X; when "emixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix_J" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam or Matrix packages to calculate the determinant; \dQuote{Matrix} and \dQuote{spam_update} provide updating Cholesky decomposition methods; "LU" provides an alternative sparse matrix decomposition approach. In addition, there are "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods; the Smirnov/Anselin (2009) trace approximation is available as "moments". Three methods: "SE_classic", "SE_whichMin", and "SE_interp" are provided experimentally, the first to attempt to emulate the behaviour of Spatial Econometrics toolbox ML fitting functions. All use grids of log determinant values, and the latter two attempt to ameliorate some features of "SE_classic".
}
  \item{quiet}{default NULL, use !verbose global option value; if FALSE, reports function values during optimization.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{interval}{default is NULL, search interval for autoregressive parameter}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{llprof}{default NULL, can either be an integer, to divide the feasible ranges into a grid of points, or a two-column matrix of spatial coefficient values, at which to evaluate the likelihood function}
  \item{trs1, trs2}{default NULL, if given, vectors for each weights object of powered spatial weights matrix traces output by \code{trW}; when given, used in some Jacobian methods}
  \item{interval1, interval2}{default is NULL, search intervals for each weights object for autoregressive parameters}
  \item{trs}{default NULL, if given, a vector of powered spatial weights matrix traces output by \code{trW}; when given, insert the asymptotic analytical values into the numerical Hessian instead of the approximated values; may be used to get around some problems raised when the numerical Hessian is poorly conditioned, generating NaNs in subsequent operations; the use of trs is recommended}
 \item{control}{list of extra control arguments - see section below}
  \item{object}{\code{Sarlm} object from \code{lagsarlm}, \code{errorsarlm} or \code{sacsarlm}}
  \item{correlation}{logical; if 'TRUE', the correlation matrix of the estimated parameters including sigma is returned and printed (default=FALSE)}
  \item{Nagelkerke}{if TRUE, the Nagelkerke pseudo R-squared is reported}
  \item{Hausman}{if TRUE, the results of the Hausman test for error models are reported}
  \item{adj.se}{if TRUE, adjust the coefficient standard errors for the number of fitted coefficients}
  \item{x}{\code{Sarlm} object from \code{lagsarlm}, \code{errorsarlm} or \code{sacsarlm} in \code{print.Sarlm}, summary object from \code{summary.Sarlm} for 
\code{print.summary.Sarlm}}
  \item{digits}{the number of significant digits to use when printing}
  \item{signif.stars}{logical. If TRUE, "significance stars" are printed
          for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The asymptotic standard error of \eqn{\rho}{rho} is only computed when
method=\dQuote{eigen}, because the full matrix operations involved would be costly
for large n typically associated with the choice of method="spam" or "Matrix". The same applies to the coefficient covariance matrix. Taken as the
asymptotic matrix from the literature, it is typically badly scaled, and with the elements involving \eqn{\rho}{rho} (lag model) or \eqn{\lambda}{lambda} (error model) being very small,
while other parts of the matrix can be very large (often many orders
of magnitude in difference). It often happens that the \code{tol.solve}
argument needs to be set to a smaller value than the default, or the RHS variables can be centred or reduced in range.

Versions of the package from 0.4-38 include numerical Hessian values where asymptotic standard errors are not available. This change has been introduced to permit the simulation of distributions for impact measures. The warnings made above with regard to variable scaling also apply in this case.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, 
p. 564). This differs from other software, including GeoDa, which does not use 
knowledge of the response variable in making predictions for the fitting data. 
Refer to the help page of \code{\link{predict.Sarlm}} for discussions and 
references.

Because numerical optimisation is used to find the values of lambda and rho in \code{sacsarlm}, care needs to be shown. It has been found that the surface of the 2D likelihood function often forms a \dQuote{banana trench} from (low rho, high lambda) through (high rho, high lambda) to (high rho, low lambda) values. In addition, sometimes the banana has optima towards both ends, one local, the other global, and conseqently the choice of the starting point for the final optimization becomes crucial. The default approach is not to use just (0, 0) as a starting point, nor the (rho, lambda) values from \code{gstsls}, which lie in a central part of the \dQuote{trench}, but either four values at (low rho, high lambda), (0, 0), (high rho, high lambda), and (high rho, low lambda), and to use the best of these start points for the final optimization. Optionally, nine points can be used spanning the whole (lower, upper) space.
}

\section{Control arguments}{
\describe{
  \item{tol.opt:}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance, a larger root may be used needed, see help(boston) for an example)}
  \item{returnHcov:}{(error model) default TRUE, return the Vo matrix for a spatial Hausman test}
  \item{pWOrder:}{(error model) default 250, if returnHcov=TRUE and the method is not \dQuote{eigen}, pass this order to \code{powerWeights} as the power series maximum limit}
  \item{fdHess:}{default NULL, then set to (method != "eigen") internally; use \code{fdHess} to compute an approximate Hessian using finite differences when using sparse matrix methods; used to make a coefficient covariance matrix when the number of observations is large; may be turned off to save resources if need be}
  \item{optimHess:}{default FALSE, use \code{fdHess} from \pkg{nlme}, if TRUE, use \code{optim} to calculate Hessian at optimum}
  \item{optimHessMethod:}{default \dQuote{optimHess}, may be \dQuote{nlm} or one of the \code{optim} methods}
  \item{compiled_sse:}{default FALSE; logical value used in the log likelihood function to choose compiled code for computing SSE}
  \item{Imult:}{default 2; used for preparing the Cholesky decompositions for updating in the Jacobian function}
  \item{super:}{if NULL (default), set to FALSE to use a simplicial decomposition for the sparse Cholesky decomposition and method \dQuote{Matrix_J}, set to  \code{as.logical(NA)} for method \dQuote{Matrix}, if TRUE, use a supernodal decomposition}
  \item{cheb_q:}{default 5; highest power of the approximating polynomial for the Chebyshev approximation}
  \item{MC_p:}{default 16; number of random variates}
  \item{MC_m:}{default 30; number of products of random variates matrix and spatial weights matrix}
  \item{spamPivot:}{default \dQuote{MMD}, alternative \dQuote{RCM}}
  \item{in_coef}{default 0.1, coefficient value for initial Cholesky decomposition in \dQuote{spam_update}}
  \item{type}{default \dQuote{MC}, used with method \dQuote{moments}; alternatives \dQuote{mult} and \dQuote{moments}, for use if \code{trs} is missing, \code{\link{trW}}}
  \item{correct}{default TRUE, used with method \dQuote{moments} to compute the Smirnov/Anselin correction term}
  \item{trunc}{default TRUE, used with method \dQuote{moments} to truncate the Smirnov/Anselin correction term}
  \item{SE_method}{default \dQuote{LU}, may be \dQuote{MC}}
  \item{nrho}{default 200, as in SE toolbox; the size of the first stage lndet grid; it may be reduced to for example 40}
  \item{interpn}{default 2000, as in SE toolbox; the size of the second stage lndet grid}
  \item{small_asy}{default TRUE; if the method is not \dQuote{eigen}, use asymmetric covariances rather than numerical Hessian ones if n <= small}
  \item{small}{default 1500; threshold number of observations for asymmetric covariances when the method is not \dQuote{eigen}}
  \item{SElndet}{default NULL, may be used to pass a pre-computed SE toolbox style matrix of coefficients and their lndet values to the "SE_classic" and "SE_whichMin" methods}
  \item{LU_order}{default FALSE; used in \dQuote{LU_prepermutate}, note warnings given for \code{lu} method}
  \item{pre_eig}{default NULL; may be used to pass a pre-computed vector of eigenvalues}
  \item{OrdVsign}{default 1; used to set the sign of the final component to negative if -1 (alpha times ((sigma squared) squared) in Ord (1975) equation B.1).}
  \item{opt_method:}{default \dQuote{nlminb}, may be set to \dQuote{L-BFGS-B} to use box-constrained optimisation in \code{optim}}
  \item{opt_control:}{default \code{list()}, a control list to pass to \code{nlminb} or \code{optim}}
  \item{pars:}{default \code{NULL}, for which five trial starting values spanning the lower/upper range are tried and the best selected, starting values of \eqn{\rho}{rho} and \eqn{\lambda}{lambda}}
  \item{npars}{default integer \code{4L}, four trial points; if not default value, nine trial points}
  \item{pre_eig1, pre_eig2}{default NULL; may be used to pass pre-computed vectors of eigenvalues}
}}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126;
Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}
(Dordrecht: Kluwer); Anselin, L. 1995 SpaceStat, a software program for
the analysis of spatial data, version 1.80. Regional Research Institute,
West Virginia University, Morgantown, WV;
Anselin L, Bera AK (1998) Spatial dependence in linear regression models
with an introduction to spatial econometrics. In: Ullah A, Giles DEA
(eds) Handbook of applied economic statistics. Marcel Dekker, New York,
pp. 237-289; Nagelkerke NJD (1991) A note on a general definition of the 
coefficient of determination. Biometrika 78: 691-692; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York; LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the Jacobian in Gaussian spatial autoregressive models: An illustrated comparison of available methods. \emph{Geographical Analysis}, 45(2), 150-179.

}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, with thanks to Andrew 
Bernat for contributions to the asymptotic standard error code.}

\seealso{\code{\link{lm}}, \code{\link{impacts}}}

\examples{
data(oldcol, package="spdep")
listw <- spdep::nb2listw(COL.nb, style="W")
ev <- eigenw(listw)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", quiet=FALSE, control=list(pre_eig=ev, OrdVsign=1))
(x <- summary(COL.lag.eig, correlation=TRUE))
coef(x)
\dontrun{
COL.lag.eig$fdHess
COL.lag.eig$resvar
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", control=list(pre_eig=ev, OrdVsign=-1))
summary(COL.lag.eigb)
COL.lag.eigb$fdHess
COL.lag.eigb$resvar
# force numerical Hessian
COL.lag.eig1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25))
summary(COL.lag.eig1)
COL.lag.eig1$fdHess
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25), trs=trMatc)
summary(COL.lag.eig1a)
COL.lag.eig1a$fdHess
COL.lag.eig$resvar[2,2]
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb$resvar[2,2]
# force numerical Hessian
COL.lag.eig1$fdHess[1,1]
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a$fdHess[2,2]
}
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb), method="Matrix", quiet=FALSE))
summary(COL.lag.M)
impacts(COL.lag.M, listw=spdep::nb2listw(COL.nb))
\dontrun{
system.time(COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="spam", quiet=FALSE))
summary(COL.lag.sp)
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"), control=list(pre_eig=ev))
summary(COL.lag.B)
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9,
 control=list(pre_eig=ev))
summary(COL.mixed.B)
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, type="mixed", control=list(pre_eig=ev))
summary(COL.mixed.W)
COL.mixed.D00 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=TRUE, control=list(pre_eig=ev))
summary(COL.mixed.D00)
COL.mixed.D01 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=FALSE, control=list(pre_eig=ev))
summary(COL.mixed.D01)
COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC + HOVAL, control=list(pre_eig=ev))
summary(COL.mixed.D1)
f <- CRIME ~ INC + HOVAL
COL.mixed.D2 <- lagsarlm(f, data=COL.OLD, listw,
 Durbin=as.formula(delete.response(terms(f))),
 control=list(pre_eig=ev))
summary(COL.mixed.D2)
COL.mixed.D1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC, control=list(pre_eig=ev))
summary(COL.mixed.D1a)
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ inc + HOVAL, control=list(pre_eig=ev)))
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ DISCBD + HOVAL, control=list(pre_eig=ev)))
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, na.action=na.exclude)
COL.lag.NA$na.action
COL.lag.NA
resid(COL.lag.NA)
COL.lag.NA1 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, Durbin=~INC) # https://github.com/r-spatial/spatialreg/issues/10
COL.lag.NA1$na.action
COL.lag.NA2 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, Durbin=~INC, na.action=na.exclude)
COL.lag.NA2$na.action
# https://github.com/r-spatial/spatialreg/issues/11
COL.lag.NA3 <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 listw, control=list(pre_eig=ev))
COL.lag.NA3$na.action
}

\dontrun{
data(boston, package="spData")
gp2mM <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, spdep::nb2listw(boston.soi), type="mixed", method="Matrix")
summary(gp2mM)
W <- as(spdep::nb2listw(boston.soi), "CsparseMatrix")
trMatb <- trW(W, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, spdep::nb2listw(boston.soi), type="mixed", method="Matrix", 
trs=trMatb)
summary(gp2mMi)
}
lw <- spdep::nb2listw(COL.nb, style="W")
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, quiet=FALSE, control=list(pre_eig=ev))
summary(COL.errW.eig)
COL.errW.eig_ev <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, control=list(pre_eig=ev))
all.equal(coefficients(COL.errW.eig), coefficients(COL.errW.eig_ev))
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb, style="B"))
summary(COL.errB.eig)
W <- as(spdep::nb2listw(COL.nb), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.errW.M <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", quiet=FALSE, trs=trMatc)
summary(COL.errW.M)
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, etype="emixed", control=list(pre_eig=ev))
summary(COL.SDEM.eig)
\dontrun{
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, Durbin=TRUE, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
COL.SDEM.eig <- errorsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 lw, Durbin=~INC, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
summary(impacts(COL.SDEM.eig))
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 spdep::nb2listw(COL.nb), na.action=na.exclude)
COL.err.NA$na.action
COL.err.NA
resid(COL.err.NA)
lw <- spdep::nb2listw(COL.nb, style="W")
print(system.time(ev <- eigenw(similar.listw(lw))))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev))))
ocoef <- coefficients(COL.errW.eig)
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev, LAPACK=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev, compiled_sse=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=as.logical(NA)))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=as.logical(NA)))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="MMD"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="RCM"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="MMD"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="RCM"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
}
listw <- spdep::nb2listw(COL.nb, style="W")
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.sacW.eig)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
set.seed(1)
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 type="sacmixed", control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW.eig)
set.seed(1)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW1.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=TRUE, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW1.eig)
set.seed(1)
summary(impacts(COL.msacW1.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW2.eig <- sacsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD, 
 listw, Durbin= ~ INC, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW2.eig)
summary(impacts(COL.msacW2.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
\dontrun{
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb), type="mixed", method="eigen")
summary(COL.mix.eig, correlation=TRUE, Nagelkerke=TRUE)
COL.mix.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 spdep::nb2listw(COL.nb), type="mixed", method="Matrix")
summary(COL.mix.M, correlation=TRUE, Nagelkerke=TRUE)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
  spdep::nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, correlation=TRUE, Nagelkerke=TRUE, Hausman=TRUE)
}
}

