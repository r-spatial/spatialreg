\name{summary.sphet}
\alias{summary.sphet}
\alias{sumres}

\title{print method for class sphet}

\description{
Method used to print objects of class \code{'summary.sphet'} and \code{'sphet'}
}

\usage{
\method{summary}{sphet}(object,width=getOption("width"), digits=getOption("digits"),obsinfo=FALSE,...)
}

\arguments{
\item{object}{an object of class 'sphet'}
\item{width}{controls the maximum number of columns on a line used in printing}
\item{digits}{minimal number of significant digits, see \code{print.default}}
\item{obsinfo}{for objects of class \code{distance}: if \code{TRUE} prints observation-wise information}
\item{...}{additional arguments to be passed}
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
