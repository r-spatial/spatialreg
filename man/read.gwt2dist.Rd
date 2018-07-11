\name{read.gwt2dist}
\alias{read.gwt2dist}
\title{Read distance ojbects}

\description{ The function reads "GWT" files (i.e. generated using \link{distance}. It will read also other more general formats
(as for example .txt files). 
}

\usage{
read.gwt2dist(file, region.id=NULL, skip=1)
}

\arguments{
  \item{file}{name of file to be read}
  \item{region.id}{variable that defines the ordering of the observations }
  \item{skip}{skip number of lines}
}


\details{
The first line of a \code{'GWT'} file generally contains some information (e.g. the name of the shape file, the number of observations),
in which case, \code{skip} should be equal to 1.
When the \code{file} has a \code{'GWT'} extension, the number of observations is generally retrived from the first line. Alternatively, it is fixed to the length
of the \link{unique} \code{region.id} variable.   
}


\value{
An object of class \code{distance}.
}


\author{ Gianfranco Piras \email{gpiras@mac.com}}

\examples{
\dontrun{
data(columbus)
dist<-read.gwt2dist(file=system.file('extdata/knn10columbus.GWT',
 package="sphet"), region.id=columbus$POLYID)
}
}
\keyword{spatial}
