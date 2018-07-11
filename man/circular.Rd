\name{circular}
\alias{circular}
\title{Generate "circular" weighting matrices }

\description{ The function generates k-ahead and k-behind weighting matrices as in Kelejian and Prucha (1999)
}

\usage{
circular(nrow, ncol, ab)
}

\arguments{
  \item{nrow}{number of rows}
  \item{ncol}{number of columns}
  \item{ab}{ahead - behind}
}



\value{
An object of class \code{nb}.
}


\author{ Gianfranco Piras \email{gpiras@mac.com}}

\examples{
ab.4<-circular(10, 10, 3)
}
\keyword{spatial}
