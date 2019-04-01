\name{set.ZeroPolicyOption}
\alias{set.VerboseOption}
\alias{get.VerboseOption}
\alias{set.ZeroPolicyOption}
\alias{get.ZeroPolicyOption}
%\alias{set.listw_is_CsparseMatrix_Option}
%\alias{get.listw_is_CsparseMatrix_Option}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Control checking of spatial object IDs}
\description{
  Provides support for checking the mutual integrity of spatial neighbour weights and spatial data; similar mechanisms are used for passing global verbose and zero.policy options, and for providing access to a running cluster for embarrassingly parallel tasks. 
}
\usage{
set.VerboseOption(check)
get.VerboseOption()
set.ZeroPolicyOption(check)
get.ZeroPolicyOption()
#set.listw_is_CsparseMatrix_Option(check)
#get.listw_is_CsparseMatrix_Option()
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{check}{a logical value, TRUE or FALSE}
}
\details{
  Analysis functions will have an spChk argument by default set to NULL, and will call \code{get.spChkOption()} to get the global spatial option for whether to check or not --- this is initialised to FALSE, and consequently should not break anything. It can be changed to TRUE using \code{set.spChkOption(TRUE)}, or the spChk argument can be assigned in analysis functions. \code{spNamedVec()} is provided to ensure that rownames are passed on to single columns taken from two-dimensional arrays and data frames.
}
\value{
  \code{set.spChkOption()} returns the old logical value, \code{get.spChkOption()} returns the current logical value, and \code{chkIDs()} returns a logical value for the test lack of difference. \code{spNamedVec()} returns the selected column with the names set to the row names of the object from which it has been extracted.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
get.VerboseOption()
get.ZeroPolicyOption()
}
 
\keyword{ spatial }
