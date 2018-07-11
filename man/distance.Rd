\name{distance}
\alias{distance}
\alias{distance.matrix}
\title{Writes distance matrices}

\description{Reads points coordinates and generates objects of class \code{distance.matrix}
}

\usage{
distance(coord,region.id=NULL,output=TRUE,
type=c("NN","distance","inverse"),
measure=c("euclidean","gcircle","chebyshev","braycur","canberra"),
nn=6, cutoff=FALSE, miles=TRUE,R=NULL, shape.name=NULL,region.id.name=NULL,
firstline=FALSE,file.name=NULL)
}

\arguments{
  \item{coord}{a matrix with the (X,Y)-coordinates of the points. The first column can be the region.id variable giving the 
  ordering of the observations}
  \item{region.id}{variable that defines the ordering of the observations}
  \item{output}{when TRUE (default) writes the object to a file}
  \item{type}{one of \code{("NN","distance","inverse")}. Nearest neighbors, distance or inverse distance }
  \item{measure}{one of \code{("euclidean","gcircle","chebyshev","braycur","canberra")}.The distance measure to be employed in the calculations.}
  \item{nn}{the number of nearest neighbors}
  \item{cutoff}{If type is \code{distance} or \code{inverse}. Assumes values 1, 2 or 3. 
  When 1, the cutoff is set to the first quantile of the distribution of distances. When 2 to the median, and when 3 to 
  the third quantile. Only observations with distance less than cutoff distance are neighbors. }
  \item{miles}{If TRUE (default), distances are in miles, otherwise in Km. (See \code{\link[sp]{spDists}} wcich returns km, and are converted if required) }
  \item{R}{ deprecated, \code{\link[sp]{spDists}} uses an approximation to the WGS84 spheroid}
  \item{shape.name}{The name of the shape file. See Details}
   \item{region.id.name}{The name of the \code{region.id} variable. See Details}
   \item{firstline}{If \code{TRUE}, a first line is added to the output file. See Details}
    \item{file.name}{If \code{output}, the name of the output file. See Details}
}


\details{
The object created is similar to the content of a \code{'GWT'} file. The output file can be of any format. 
In particular, it could be a \code{'GWT'} file. When \code{firstline} is \code{TRUE}, an header line is added to the \code{'GWT'} file. 
The first element is simply a place holder, the second is the number of observations. The name of the shape file 
and of the id variable can be specified by the options \code{shape.name} and \code{region.id.name} respectively. 
The function performs a series of test on the \code{region.id} variable. If a \code{region.id} variable is not specified and 
\code{coord} only has two columns, a sequence from 1 to the number of observations is generated and 
used as identification variable. If \code{region.id} is specified 
and the first column of \code{coord} contains an id variable they should be the same. 

}


\value{
A \code{matrix} of three columns: \code{from}, \code{to}, and \code{distance}
}


\author{ Gianfranco Piras \email{gpiras@mac.com}}



\examples{
X<-runif(100,0,70)
Y<-runif(100,-30,20)
coord1<-cbind(seq(1,100),X,Y)
thm2 <- distance(coord1,region.id=NULL,output=FALSE,type="NN", nn=6)
thm2 <- distance(coord1,region.id=NULL,output=FALSE,type="distance", cutoff=1)
}
\keyword{spatial}
