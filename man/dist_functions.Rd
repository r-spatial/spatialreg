\name{dist_functions}
\alias{dist.euclidean}
\alias{dist.chebyshev}
\alias{dist.braycur}
\alias{dist.canberra}
\alias{dist.gcircle}


\title{Distance measures available in distance}

\description{
The great circle distance is calculated by the function \code{rdist.earth} in library \code{fields}.
}




\details{
The distance measures implemented in \code{sphet} are:
\itemize{
    \item \code{'euclidean'}: \eqn{\sqrt{\sum{(x_i - y_i)^2}}}
    \item \code{'chebyshev'}: \eqn{max(|x_i - y_i|)} 
    \item \code{'braycur'}: \eqn{ \frac{\sum{|x_i - y_i|}}{\sum{|x_i + y_i|}}}
	\item \code{'canberra'}:  \eqn{  \frac{\sum{|x_i - y_i|}}{\sum{|x_i| + |y_i|}} }
	\item \code{'gcircle'}: see \code{\link[sp]{spDists}}, which uses an approximation to the WGS84 spheroid. 
}
}



\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{\code{\link{gstslshet}}, \code{\link{distance}}, \code{\link{distance}}}


\keyword{models}
