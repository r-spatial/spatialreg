\name{print.sphet}
\alias{print.sphet}
\alias{print.summary.sphet}
\title{print method for class sphet}

\description{
Method used to print objects of class \code{'summary.sphet'} and \code{'sphet'}
}

\usage{
\method{print}{sphet}(x, digits = max(3,getOption("digits") -3), ...)
}

\arguments{
\item{x}{an object of class \code{'summary.sphet'} and \code{sphet}}
\item{digits}{minimal number of significant digits, see \code{\link{print.default}}}
\item{...}{additional argument to be passed}
}


\details{
The summary function summary.sphet 
returns an objects of class 'sphet' 
organized in a coefficient matrix.

}

\seealso{
  \code{\link{gstslshet}}, \code{\link{stslshac}}}

\author{Gianfranco Piras\email{gpiras@mac.com}}

\examples{
data(columbus)
listw<-nb2listw(col.gal.nb)
res<-gstslshet(CRIME~HOVAL + INC, data=columbus, listw=listw)
summary(res)

}
