\name{kpjtest}
\alias{kpjtest}
\title{Kelejian and Piras J-test}

\description{ The function calculate the Kelejian and Piras J-test for spatial models.
Both models (under the null and under the alternative) 
can be specified with additional endogenous variables, and additional instruments. 
The model under the null allows for heteroskedasticity as well as spatial autocorrelation:

\deqn{y=\lambda W y + X \beta + u  }
\deqn{u=Re}  with \deqn{e ~ N(0,\sigma^2_i) }

Note that when R reduces to an identity matrix, the error term, while still heteroskedastic, is not spatially autocorrelated. 
On the other hand, when the \eqn{\sigma^2_i} are all the same (and R is an identity matrix) than the error term is neither heteroskedastic nor autocorrelated.  
}

\usage{
kpjtest(H0model, H1model, data = list(), listw0 = NULL, listw1 = NULL, endogH0 = NULL, 
endogH1 = NULL, instrumentsH0 = NULL, instrumentsH1 = NULL, lag.instr = FALSE,
model = c("lag", "sarar"), het = FALSE, HAC = F, distance = NULL, 
type = c("Epanechnikov", "Triangular", "Bisquare", 
"Parzen", "QS", "TH", "Rectangular"), 
bandwidth = "variable",  na.action = na.fail)
}

\arguments{
  \item{H0model}{Formula object for the specification of the model under the null}
  \item{H1model}{Formula object for the specification of the model under the alternative}
    \item{data}{an object of class \link{data.frame}. An optional data frame containing the variables
  in the model.}
    \item{listw0}{an object of class \code{listw}, \code{matrix}, or \code{Matrix}. The spatial weighting matrix under the null model}
   \item{listw1}{an object of class \code{listw}, \code{matrix}, or \code{Matrix}. The spatial weighting matrix under the alternative model}
\item{endogH0}{additional endogenous variables under the null model. Default \code{NULL}. If not \code{NULL} should be specified as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression. }
\item{endogH1}{additional endogenous variables under the alternative model. Default \code{NULL}. If not \code{NULL} should be specified as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression. }
\item{instrumentsH0}{external instruments for the null model. Default \code{NULL}. If not \code{NULL} should be specified 
as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression.  }
\item{instrumentsH1}{external instruments for the alternative model. Default \code{NULL}. If not \code{NULL} should be specified 
as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression.  }
\item{lag.instr}{should the external instruments be spatially lagged?}
 \item{model}{one of \code{lag}, or \code{sarar}. The current version of the function only implements the lag model. } 
  \item{het}{default FALSE: if TRUE uses the methods developed for heteroskedasticity}
   \item{na.action}{a function which indicates what should happen when the data contains missing values. See \link{lm} for details.}
 \item{HAC}{perform the HAC estimator of Kelejian and Prucha, 2007 on the null (and augmented) model.}
   \item{distance}{an object of class \code{distance} created for example by \link{read.gwt2dist} 
  The object contains the specification of the distance measure 
  to be employed in the estimation of the VC matrix. See Details. }
    \item{type}{One of \code{c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular")}. 
  The type of Kernel to be used. See Details. }
    \item{bandwidth}{"variable" (default) - or numeric when a fixed bandwidth is specified by the user.}
}


\details{
In order to calculate the J-test, the function follows a few steps:
\itemize{
	\item The alternative model is estimated by S2SLS. 
	\item Based on the estimated parameters in the previous step,  obtain a prediction based on the alternative models of the dependent vector in the null model. The predictor is based on the right hand side of the model. 
	\item Use these predicted values of the dependent variable based on the alternative models into the null model to obtain the augmented model. 
	\item Estimate the augmented model by 2SLS using all of the instruments relating to the null model as well as all of the instruments relating to the alternative models. 
	\item Test for the statistical significance of the predicted value. If it is not significant, accept the null model. If it is  significant, reject the null and conclude that the true model is the alternative models. 
	}
The output is an object of class \code{sphet} where the last row of the table of coefficients is the prediction. 

When the model is heteroskedastic as well as spatially autocorrelated, an HAC procedure is employed. 
The default sets the bandwith for each observation to the maximum distance for that observation (i.e.
the max of each element of the list of distances). 

Six different kernel functions are implemented:
\itemize{
    \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
    \item \code{'Triangular'}: \eqn{K(z) = 1-z} 
    \item \code{'Bisquare'}: \eqn{K(z) = (1-z^2)^2}
    \item \code{'Parzen'}: \eqn{K(z) = 1-6z^2+6 |z|^3} if  \eqn{z \leq 0.5} and  
\eqn{ K(z) = 2(1-|z|)^3} if \eqn{0.5 < z \leq 1} 
	\item \code{'TH'} (Tukey - Hanning):  \eqn{ K(z) = \frac{1+ \cos(\pi z)}{2}}
	\item \code{'Rectangular'}:  \eqn{ K(z) = 1}
	\item \code{'QS'} (Quadratic Spectral): \eqn{K(z) = \frac{25}{12\pi^2z^2} 
  (\frac{\sin(6\pi z)/5)}{6\pi z/5} - \cos(6\pi z)/5)}). 
}

If the kernel type is not one of the six implemented, the function will terminate with an error message. 
The spatial two stage least square estimator is based on the matrix of instruments \eqn{H=[X,WX,W^2X^2]}.
}


\value{
A list object of class \code{sphet}
\item{coefficients}{Generalized Spatial two stage least squares coefficient estimates of \eqn{\delta} and GM estimator for \eqn{\rho}. }
\item{var}{variance-covariance matrix of the estimated coefficients}
\item{s2}{GS2SLS residuals variance}
\item{residuals}{GS2SLS residuals}
\item{yhat}{difference between GS2SLS residuals and response variable}
\item{call}{the call used to create this object}
\item{model}{the model matrix of data}
\item{method}{\code{'s2slshac'}}
}
\references{

Kelejian and Piras (2017). \emph{Spatial Econometrics}. Academic Press. ISBN: 978-0-12-813387-3 


Gianfranco Piras (2010). sphet: Spatial Models with Heteroskedastic Innovations in R. \emph{Journal of Statistical Software}, 35(1), 1-21. \url{http://www.jstatsoft.org/v35/i01/}.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{http://www.jstatsoft.org/v63/i18/}.
}


\author{ Gianfranco Piras \email{gpiras@mac.com}}

\examples{
data(boston)
boslw <- nb2listw(boston.soi)
Bos.Knn <- knearneigh(boston.utm, k = 5)
bos.nb <- knn2nb(Bos.Knn)
boslw2 <- nb2listw(bos.nb)
fm <- log(MEDV) ~ CRIM + ZN + INDUS + CHAS
fm2 <- log(MEDV) ~ CRIM + ZN + INDUS + RM + AGE 
test <- kpjtest(fm, fm2, data = boston.c, listw0 = boslw, listw1 = boslw2, model = "lag")
}
\keyword{spatial}
