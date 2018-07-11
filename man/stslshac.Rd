\name{stslshac}
\alias{stslshac}
\alias{tsls}

\title{Spatial two stages least square with HAC standard errors}

\description{Non-parametric heteroskedasticity and autocorrelation 
consistent (HAC) estimator of the variance-covariance (VC) for a vector of sample moments within a 
spatial context. The disturbance vector is generated as follows:
\deqn{ u = R  \epsilon }
where \eqn{R} is a non-stochastic matrix.
}

\usage{
stslshac(formula, data = list(),listw, na.action = na.fail, zero.policy = NULL, 
HAC=TRUE, distance=NULL,type = c("Epanechnikov", "Triangular", "Bisquare",
"Parzen", "QS","TH","Rectangular"), 
bandwidth = "variable", W2X = TRUE)
}

\arguments{
  \item{formula}{a description of the model to be fit }
  \item{data}{an object of class \link{data.frame}. An optional data frame containing the variables
  in the model.}
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw} }
  \item{distance}{an object of class \code{distance} created for example by \link{read.gwt2dist} 
  The object contains the specification of the distance measure 
  to be employed in the estimation of the VC matrix. See Details. }
  \item{type}{One of \code{c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH")}. 
  The type of Kernel to be used. See Details. }
  \item{na.action}{a function which indicates what should happen when the data contains missing values.
  See \link{lm} for details.}
  \item{zero.policy}{See \code{lagsarlm} for details}
  \item{bandwidth}{"variable" (default) - or numeric when a fixed bandwidth is specified by the user.}
  \item{HAC}{if FALSE traditional standard errors are provided.}
  \item{W2X}{default TRUE. if FALSE only WX are used as instruments in the spatial two stage least squares.}
}


\details{
The default sets the bandwith for each observation to the maximum distance for that observation (i.e.
the max of each element of the list of distances). 

Six different kernel functions are implemented:
\itemize{
    \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
    \item \code{'Rectangular'}: \eqn{K(z) = 1}
    \item \code{'Triangular'}: \eqn{K(z) = 1-z} 
    \item \code{'Bisquare'}: \eqn{K(z) = (1-z^2)^2}
    \item \code{'Parzen'}: \eqn{K(z) = 1-6z^2+6 |z|^3} if  \eqn{z \leq 0.5} and  
\eqn{ K(z) = 2(1-|z|)^3} if \eqn{0.5 < z \leq 1} 
	\item \code{'TH'} (Tukey - Hanning):  \eqn{ K(z) = \frac{1+ \cos(\pi z)}{2}}
	\item \code{'QS'} (Quadratic Spectral): \eqn{K(z) = \frac{25}{12\pi^2z^2} 
  (\frac{\sin(6\pi z)/5)}{6\pi z/5} - \cos(6\pi z)/5)}). 
}

If the kernel type is not one of the six implemented, the function will terminate with an error message. 
The spatial two stage least square estimator is based on the matrix of instruments \eqn{H=[X,WX,W^2X^2]}.

}


\value{
A list object of class \code{sphet}
\item{coefficients}{Spatial two stage least squares coefficient estimates}
\item{vcmat}{variance-covariance matrix of the estimated coefficients}
\item{s2}{S2sls residulas variance}
\item{residuals}{S2sls residuals}
\item{yhat}{difference between residuals and response variable}
\item{call}{the call used to create this object}
\item{model}{the model matrix of data}
\item{type}{the kernel employed in the estimation}
\item{bandwidth}{the type of bandwidth}
\item{method}{\code{'s2slshac'}}


}

\references{
      Kelejian, H.H. and Prucha, I.R. (2007) 
HAC estimation in a spatial framework,
    \emph{Journal of Econometrics}, \bold{140}, pages 131--154.


  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1998) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.


}

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{\code{\link{gstslshet}}, \code{\link{distance}}, \code{\link{distance}}}

\examples{
data(columbus)
listw<-nb2listw(col.gal.nb)
data(coldis)
res<-stslshac(CRIME~HOVAL + INC, data=columbus,listw=listw,  distance=coldis, type='Triangular')
summary(res)
}

\keyword{models}
