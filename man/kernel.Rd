\name{kernel}
\alias{kernel}
\alias{bisq}
\alias{epan}
\alias{parzen}
\alias{qs}
\alias{th}
\alias{triangular}
\alias{rectangular}

\title{Spatial two stages least square with HAC standard errors}

\description{Main functions to calculate the kernels in \code{stslshac}
}




\details{
Six different kernel functions are implemented:
\itemize{
    \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
    \item \code{'Rectangular'}: \eqn{K(z) = 1-z^2}
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



\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{\code{\link{gstslshet}}, \code{\link{distance}}, \code{\link{distance}}}


\keyword{models}
