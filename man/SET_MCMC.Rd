\name{spBreg_lag}
\alias{spBreg_lag}
\alias{spBreg_err}
\alias{spBreg_sac}
\alias{impacts.MCMC_sar_g}
\alias{impacts.MCMC_sem_g}
\alias{impacts.MCMC_sac_g}

\title{Bayesian MCMC spatial simultaneous autoregressive model estimation}
\description{
The \code{spBreg_lag} function is an early-release version of the Matlab Spatial Econometrics Toolbox function \code{sar_g.m}, using drawing by inversion, and not accommodating heteroskedastic disturbances.
}
\usage{
spBreg_lag(formula, data = list(), listw, na.action, Durbin, type,
    zero.policy=NULL, control=list())
spBreg_sac(formula, data = list(), listw, listw2=NULL, na.action, 
    Durbin, type, zero.policy=NULL, control=list())
spBreg_err(formula, data = list(), listw, na.action, Durbin, etype,
    zero.policy=NULL, control=list())
\method{impacts}{MCMC_sar_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
\method{impacts}{MCMC_sem_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
\method{impacts}{MCMC_sac_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
}


\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw, listw2}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{Durbin}{default FALSE (spatial lag model); if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
  \item{type, etype}{(use the \sQuote{Durbin=} argument - retained for backwards compatibility only) default "lag", may be set to "mixed"; when "mixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included; \dQuote{Durbin} may be used instead of \dQuote{mixed}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA}

 \item{control}{list of extra control arguments - see section below}
  \item{obj}{A spatial regression object}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
}


\section{Control arguments}{
\describe{
  \item{tol.opt:}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance, a larger root may be used needed, see help(boston) for an example)}
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
}}

\section{Extra Bayesian control arguments}{
\describe{
  \item{ldet_method}{default \dQuote{SE_classic}; equivalent to the \code{method} argument in \code{lagsarlm}}
  \item{interval}{default \code{c(-1, 1)}; used unmodified or set internally by \code{jacobianSetup}}
  \item{ndraw}{default \code{2500L}; integer total number of draws}
  \item{nomit}{default \code{500L}; integer total number of omitted burn-in draws}
  \item{thin}{default \code{1L}; integer thinning proportion}
  \item{verbose}{default \code{FALSE}; inverse of \code{quiet} argument in \code{lagsarlm}}
  \item{detval}{default \code{NULL}; not yet in use, precomputed matrix of log determinants}
  \item{prior}{a list with the following components:
    \describe{
      \item{rhoMH, lambdaMH}{default FALSE; use Metropolis or griddy Gibbs}
      \item{Tbeta}{default \code{NULL}; values of the betas variance-covariance matrix, set to \code{diag(k)*1e+12} if \code{NULL}}
      \item{c_beta}{default \code{NULL}; values of the betas set to 0 if \code{NULL}}
      \item{rho}{default \code{0.5}; value of the autoregressive coefficient}
      \item{sige}{default \code{1}; value of the residual variance}
      \item{nu}{default \code{0}; informative Gamma(nu,d0) prior on sige}
      \item{d0}{default \code{0}; informative Gamma(nu,d0) prior on sige}
      \item{a1}{default \code{1.01}; parameter for beta(a1,a2) prior on rho}
      \item{a2}{default \code{1.01}; parameter for beta(a1,a2) prior on rho}
      \item{cc}{default \code{0.2}; initial tuning parameter for M-H sampling}
      \item{gG_sige}{default TRUE; include sige in lambda griddy Gibbs update}
      \item{cc1}{default \code{0.2}; initial tuning parameter for M-H sampling}
      \item{cc2}{default \code{0.2}; initial tuning parameter for M-H sampling}
    }}
}}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}, with thanks to Abhirup Mallik and Virgilio Gómez-Rubio for initial coding GSoC 2011}

\examples{
#require("spdep", quietly=TRUE)
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb, style="W")
ev <- eigenw(lw)
W <- as(lw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
require("coda", quietly=TRUE)
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
print(summary(COL.err.Bayes))
print(raftery.diag(COL.err.Bayes, r=0.01))
\dontrun{
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE)
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC)
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.sacW.B0 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=FALSE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B0))
print(summary(impacts(COL.sacW.B0, tr=trMatc), zstats=TRUE, short=TRUE))
set.seed(1)
COL.sacW.B1 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B1))
print(summary(impacts(COL.sacW.B1, tr=trMatc), zstats=TRUE, short=TRUE))
set.seed(1)
COL.lag.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=lw)
print(summary(COL.lag.Bayes))
print(summary(impacts(COL.lag.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
print(summary(impacts(COL.lag.Bayes, evalues=ev), short=TRUE, zstats=TRUE))
set.seed(1)
COL.D0.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=lw, Durbin=TRUE)
print(summary(COL.D0.Bayes))
print(summary(impacts(COL.D0.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
set.seed(1)
COL.D1.Bayes <- spBreg_lag(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 listw=lw, Durbin= ~ INC)
print(summary(COL.D1.Bayes))
print(summary(impacts(COL.D1.Bayes, tr=trMatc), short=TRUE, zstats=TRUE))
#data(elect80, package="spData")
#lw <- spdep::nb2listw(e80_queen, zero.policy=TRUE)
#el_ml <- lagsarlm(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
# + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE, method="LU")
#print(summary(el_ml))
#set.seed(1)
#el_B <- spBreg_lag(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
# + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE)
#print(summary(el_B))
#print(el_ml$timings)
#print(attr(el_B, "timings"))
}
}
