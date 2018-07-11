\name{gstslshet}
\alias{gstslshet}
\alias{impacts.gstsls}


\title{GM estimation of a Cliff-Ord type model with Heteroskedastic Innovations }

\description{Multi step GM/IV estimation of a linear Cliff and Ord -type of model
of the form:

\deqn{y=\lambda W y + X \beta + u  }
\deqn{u=\rho W u + e}  with \deqn{e ~ N(0,\sigma^2_i) }

The model allows for spatial lag in the dependent variable
and disturbances. The innovations in the disturbance process are assumed 
heteroskedastic of an unknown form.
}

\usage{
gstslshet(formula, data=list(), listw, na.action=na.fail, 
zero.policy=NULL,initial.value=0.2, abs.tol=1e-20, 
rel.tol=1e-10, eps=1e-5, inverse=T,sarar=T)
\method{impacts}{gstsls}(obj, \dots, tr, R = NULL, listw = NULL,
tol = 1e-06, empirical = FALSE, Q=NULL)
}

\arguments{
  \item{formula}{a description of the model to be fit }
  \item{data}{an object of class \link{data.frame}. An optional data frame containing the variables
  in the model.}
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw} }
  \item{na.action}{a function which indicates what should happen when the data contains missing values.
  See \link{lm} for details.}
  \item{zero.policy}{See \code{lagsarlm} for details}
  \item{initial.value}{The initial value for \eqn{\rho}. It can be either numeric (default is 0.2) or
  set to \code{'SAR'}, in which case the optimization will start from the estimated coefficient of a regression of the 2SLS 
  residuals over their spatial lag (i.e. a spatial AR model)}
  \item{abs.tol}{Absolute tolerance. See \link{nlminb} for details.}
    \item{rel.tol}{Relative tolerance. See \link{nlminb} for details.}
    \item{eps}{Tolerance level for the approximation. See Details.}
    \item{inverse}{\code{TRUE}. If \code{FALSE}, an appoximated inverse is calculated. See Details.}
    \item{sarar}{\code{TRUE}. If \code{FALSE}, a spatial error model is estimated. }
  \item{obj}{A gstsls spatial regression object created by \code{gstslshet}}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
}

\details{
The procedure consists of two steps alternating GM and IV estimators. Each step consists of sub-steps.
In step one \eqn{\delta = [\beta',\lambda]'} is estimated by 2SLS. The 2SLS residuals are first employed
to obtain an initial (consistent but not efficient) GM estimator of \eqn{\rho} and then a consistent and efficient 
estimator (involving the variance-covariance matrix of the limiting distribution of the normalized sample moments). 
In step two, the spatial Cochrane-Orcutt transformed model is estimated by 2SLS. This corresponds to a GS2SLS procedure. 
The GS2SLS residuals are used to obtain a consistent and efficient GM estimator for \eqn{\rho}. 

The initial value for the optimization in step 1b is taken to be \code{initial.value}. The initial value in step 1c is the 
optimal parameter of step 1b. Finally, the initial value for the optimization of step 2b is the optimal parameter of step 1c.

Internally, the object of class \code{listw} is transformed into a \link{Matrix} 
using the function \link{listw2dgCMatrix}.


The expression of the estimated variance covariance matrix of the limiting 
distribution of the normalized sample moments based on 2SLS residuals 
involves the inversion of \eqn{I-\rho W'}.
When \code{inverse} is \code{FALSE}, the inverse is calculated using the approximation 
\eqn{I +\rho W' + \rho^2 W'^2 + ...+ \rho^n W'^n}. 
The powers considered depend on a condition. 
The 
function will keep adding terms until the absolute value of the \code{sum} of all elements 
of the matrix \eqn{\rho^i W^i} is greater than a fixed \eqn{\epsilon} (\code{eps}). By default \code{eps}
is set to 1e-5.}





\value{
A list object of class \code{sphet}
\item{coefficients}{Generalized Spatial two stage least squares coefficient estimates of \eqn{\delta} and GM estimator for \eqn{\rho}. }
\item{var}{variance-covariance matrix of the estimated coefficients}
\item{s2}{GS2SLS residuals variance}
\item{residuals}{GS2SLS residuals}
\item{yhat}{difference between GS2SLS residuals and response variable}
\item{call}{the call used to create this object}
\item{model}{the model matrix of data}
\item{method}{\code{'gs2slshac'}}
\item{W}{Wald test for both \eqn{\rho} and \eqn{\lambda} are zero}
}

\seealso{ \code{\link{stslshac}}}

\references{

Arraiz, I. and Drukker, M.D. and Kelejian, H.H. and Prucha, I.R. (2007) 
A spatial Cliff-Ord-type Model with Heteroskedastic Innovations: Small and Large Sample Results,
\emph{Department of Economics, University of Maryland}'

      Kelejian, H.H. and Prucha, I.R. (2007) 
Specification and Estimation of Spatial Autoregressive Models with Autoregressive and Heteroskedastic Disturbances,
    \emph{Journal of Econometrics}, forthcoming.

  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1998) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.


}

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\examples{
data(columbus)
listw<-nb2listw(col.gal.nb)
res<-gstslshet(CRIME~HOVAL + INC, data=columbus, listw=listw)
summary(res)
}
\keyword{models}
