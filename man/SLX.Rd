\name{lmSLX}
\alias{lmSLX}
\alias{create_WX}
\alias{impacts.SlX}
\alias{print.WXimpact}
\alias{summary.WXimpact}
\alias{print.summary.WXimpact}
\alias{predict.SlX}

\title{Spatial Durbin linear (SLX, spatially lagged X) model}
\description{
\code{lmSLX} fits an \code{lm} model augmented with the spatially lagged RHS variables, including the lagged intercept when the spatial weights are not row-standardised. \code{create_WX} creates spatially lagged RHS variables, and is exposed for use in model fitting functions.
}
\usage{
lmSLX(formula, data = list(), listw, na.action, weights=NULL,
 Durbin=TRUE, zero.policy=NULL)
create_WX(x, listw, zero.policy=NULL, prefix="")
\method{impacts}{SlX}(obj, ...)
\method{print}{WXimpact}(x, ...)
\method{summary}{WXimpact}(object, ..., adjust_k=(attr(object, "type") == "SDEM"))
\method{predict}{SlX}(object, newdata, listw, zero.policy=NULL, ...)

}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the spatial weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{weights}{an optional vector of weights to be used in the fitting process. Non-NULL weights can be used to indicate that different observations have different variances (with the values in weights being inversely proportional to the variances); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations (including the case that there are w_i observations equal to y_i and the data have been summarized) - \code{\link{lm}}}
  \item{Durbin}{default TRUE for \code{lmSLX} (Durbin model including WX); if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE assign NA}
  \item{obj}{A spatial regression object created by \code{lmSLX}}
  \item{\dots}{Arguments passed through}
  \item{prefix}{default empty string, may be \dQuote{lag} in some cases}
  \item{x, object}{model matrix to be lagged; lagImpact objects created by \code{impacts} methods}
  \item{adjust_k}{default TRUE if SDEM else FALSE, adjust internal OLS SDEM standard errors by dividing by n rather than (n-k) (default changed and bug fixed after 0.7-8; standard errors now ML in SDEM summary and impacts summary and identical - for SLX use FALSE)}
  \item{newdata}{data frame in which to predict --- if NULL, predictions are
for the data on which the model was fitted. Should have row names corresponding to region.id. If row names are exactly the same than the ones used for training, it uses in-sample predictors for forecast.}
}

\value{
The \code{lmSLX} function returns an \dQuote{lm} object with a \dQuote{mixedImps} list of three impact matrixes (impacts and standard errors) for direct, indirect and total impacts; total impacts calculated using \code{gmodels::estimable()}.
}

\seealso{\code{\link{lm}}, \code{\link[gmodels]{estimable}}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\examples{
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb, style="W")
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(COL.SLX)
summary(impacts(COL.SLX))
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=TRUE)
summary(impacts(COL.SLX))
summary(COL.SLX)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=~INC)
summary(impacts(COL.SLX))
summary(COL.SLX)
COL.SLX <- lmSLX(CRIME ~ INC, data=COL.OLD, listw=lw)
summary(COL.SLX)
summary(impacts(COL.SLX))
\dontrun{
crds <- cbind(COL.OLD$X, COL.OLD$Y)
mdist <- sqrt(sum(diff(apply(crds, 2, range))^2))
dnb <- spdep::dnearneigh(crds, 0, mdist)
dists <- spdep::nbdists(dnb, crds)
f <- function(x, form, data, dnb, dists, verbose) {
  glst <- lapply(dists, function(d) 1/(d^x))
  lw <- spdep::nb2listw(dnb, glist=glst, style="B")
  res <- logLik(lmSLX(form=form, data=data, listw=lw))
  if (verbose) cat("power:", x, "logLik:", res, "\n")
  res
}
opt <- optimize(f, interval=c(0.1, 4), form=CRIME ~ INC + HOVAL,
 data=COL.OLD, dnb=dnb, dists=dists, verbose=TRUE, maximum=TRUE)
glst <- lapply(dists, function(d) 1/(d^opt$maximum))
lw <- spdep::nb2listw(dnb, glist=glst, style="B")
SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(SLX)
summary(impacts(SLX))
}
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
pslx0 <- predict(COL.SLX)
pslx1 <- predict(COL.SLX, newdata=COL.OLD, listw=lw)
all.equal(pslx0, pslx1)
COL.OLD1 <- COL.OLD
COL.OLD1$INC <- COL.OLD1$INC + 1
pslx2 <- predict(COL.SLX, newdata=COL.OLD1, listw=lw)
sum(coef(COL.SLX)[c(2,4)])
mean(pslx2-pslx1)
}

