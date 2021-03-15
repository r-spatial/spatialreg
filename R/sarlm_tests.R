logLik.Sarlm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- length(residuals(object))
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

NK.Sarlm <- function(obj) {
     n <- length(residuals(obj))
     nullLL <- obj$LLNullLlm
     if (is.null(nullLL)) return(nullLL)
     c(1 - exp(-(2/n)*(logLik(obj) - nullLL)))
}


LR.Sarlm <- function(x, y)
{
	if (!inherits(x, "logLik")) LLx <- logLik(x)
	else LLx <- x
	if (!inherits(y, "logLik")) LLy <- logLik(y)
	else LLy <- y
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
	if (parameter < 1) 
		stop("non-positive degrees of freedom: no test possible")
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	attr(estimate, "names") <- c(paste("Log likelihood of",
		deparse(substitute(x))), paste("Log likelihood of",
		deparse(substitute(y))))
	method <- "Likelihood ratio for spatial linear models"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}


LR1.Sarlm <- function(object)
{
	if (!inherits(object, "Sarlm")) stop("Not a Sarlm object")
	LLx <- logLik(object)
#	LLy <- logLik(object$lm.model)
        LLy <- object$logLik_lm.model
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
	if (parameter < 1) 
		stop("non-positive degrees of freedom: no test possible")
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	if (object$type == "error") alt <- "spatial error model"
	else alt <- "spatial lag model"
	attr(estimate, "names") <- c(paste("Log likelihood of",
		alt), paste("Log likelihood of OLS fit",
		deparse(substitute(y))))
	method <- "Likelihood Ratio diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

Wald1.Sarlm <- function(object) {
	if (!inherits(object, "Sarlm")) stop("Not a Sarlm object")
#	if (!object$ase) 
#		stop("Cannot compute Wald statistic: parameter a.s.e. missing")
	LLx <- logLik(object)
#	LLy <- logLik(object$lm.model)
        LLy <- object$logLik_lm.model
	if (object$type == "lag" || object$type == "mixed") {
		estimate <- object$rho
                rse <- object$rho.se
                if (is.null(rse)) return(rse)
		statistic <- (object$rho / rse)^2
		attr(statistic, "names") <- ifelse(is.logical(object$fdHess), 
                    "Wald statistic", "Approximate Wald statistic")
	} else {
		estimate <- object$lambda
                lse <- object$lambda.se
                if (is.null(lse)) return(lse)
		statistic <- (object$lambda / lse)^2
		attr(statistic, "names") <- ifelse(is.logical(object$fdHess), 
                    "Wald statistic", "Approximate Wald statistic")
	}
	parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
	if (parameter < 1) 
		stop("non-positive degrees of freedom: no test possible")
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	method <- "Wald diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res

}

Hausman.test.Sarlm <- function(object, ..., tol=NULL) {
    if (!inherits(object, "Sarlm")) stop("not a Sarlm object")
    if (object$type != "error") stop("not a spatial error model")
    fmeth <- ifelse(object$method != "eigen", "(approximate)", "(asymptotic)") 
    if (is.null(object$Hcov)) stop("Vo not available")
    s2 <- object$s2
    Vo <- s2 * object$Hcov
    Vs <- s2 * object$Vs
    d <- object$coef_lm.model - object$coefficients
    if (!is.null(tol)) VV <- try(solve((Vo - Vs), tol=tol))
    else VV <- try(solve(Vo - Vs))
    if (inherits(VV, "try.error")) {
        warning("(Vo - Vs) inversion failure")
        return(NULL)
    }
    statistic <- t(d) %*% VV %*% d
    attr(statistic, "names") <- "Hausman test"
    parameter <- length(d)
    attr(parameter, "names") <- "df"
    p.value <- 1 - pchisq(abs(statistic), parameter)
    method <- paste("Spatial Hausman test", fmeth)
    data.name <- strwrap(deparse(object$formula), exdent=4)
    if (length(data.name) > 1L) 
        data.name <- paste(data.name, collapse="\n    ")
    res <- list(statistic = statistic, parameter = parameter, 
        p.value = p.value, method = method, data.name=data.name)
    class(res) <- "htest"
    res
}

# Copyright 2004-2011 by Roger Bivand (original taken from bptest() in the lmtest
# package, Copyright (C) 2001 Torsten Hothorn and Achim Zeileis and released
# under GNU General Public License, Version 2 or 3.
#

bptest.Sarlm <- function (object, varformula=NULL, studentize = TRUE, data=list()) 
{
    if(!inherits(object, "Sarlm")) stop("not Sarlm object")
    Z <- object$tarX
    if (!is.null(varformula)) Z <- model.matrix(varformula, data = data)
    k <- ncol(Z)
    n <- nrow(Z)
    resi <- object$residuals
    if (length(resi) != nrow(Z))
        stop("number of residuals differs from varformula matrix rows")
    sigma2 <- sum(resi^2)/n
    if (studentize) {
        w <- resi^2 - sigma2
        fv <- lm.fit(Z, w)$fitted
        bp <- n * sum(fv^2)/sum(w^2)
        method <- "studentized Breusch-Pagan test"
    }
    else {
        f <- resi^2/sigma2 - 1
        fv <- lm.fit(Z, f)$fitted
        bp <- 0.5 * sum(fv^2)
        method <- "Breusch-Pagan test"
    }
    names(bp) <- "BP"
    df <- k - 1
    names(df) <- "df"
    RVAL <- list(statistic = bp, parameter = df, method = method, 
        p.value = 1 - pchisq(bp, df))
    class(RVAL) <- "htest"
    return(RVAL)
}

anova.Sarlm <- function(object, ...) {
    if (length(list(object, ...)) > 1L) {
        getResponseFormula <- function (object) 
        {
            form <- formula(object$call)
            if (!(inherits(form, "formula") && (length(form) == 3L))) {
                stop("\"Form\" must be a two sided formula")
            }
            eval(parse(text = paste("~", deparse(form[[2]]))))
        }

        object <- list(object, ...)
        ancall <- sys.call()
        nmodels <- length(object)
        if (nmodels == 1) return(anova(object))
        termsClass <- unlist(lapply(object, data.class))
        if (!all(match(termsClass, c("lm", "Sarlm"), 0))) {
            stop(paste("Objects must inherit from classes \"Sarlm\" or \"lm\""))
        }
        resp <- unlist(lapply(object, 
	    function(el) deparse(getResponseFormula(el)[[2]])))
        subs <- as.logical(match(resp, resp[1], FALSE))
        if (!all(subs)) 
            warning(paste("Some fitted objects deleted because", 
                "response differs from the first model"))
        if (sum(subs) == 1) 
            stop("First model has a different response from the rest")
        object <- object[subs]
        aux <- lapply(object, logLik)
        if (length(unique(unlist(lapply(object, 
	    function(el) length(residuals(el)))))) > 1L) {
            stop("All fitted objects must use the same number of observations")
        }
        dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
        logLik <- unlist(lapply(aux, function(el) c(el)))
        AIC <- unlist(lapply(aux, AIC))
        aod <- data.frame(Model = (1:nmodels), df = dfModel, 
        AIC = AIC, logLik = logLik, check.names = FALSE)
        ddf <- diff(dfModel)
        if (sum(abs(ddf)) > 0) {
	   effects <- rep("", nmodels)
	   for (i in 2:nmodels) {
                if (ddf[i - 1] != 0) {
            	    effects[i] <- paste(i - 1, i, sep = " vs ")
                }
            }
            pval <- rep(NA, nmodels - 1)
            ldf <- as.logical(ddf)
            lratio <- 2 * abs(diff(logLik))
            lratio[!ldf] <- NA
            pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
	    aod <- data.frame(aod, Test = effects, L.Ratio = c(NA, 
                lratio), "p-value" = c(NA, pval), check.names = FALSE)
    	}
    	row.names(aod) <- unlist(lapply(as.list(ancall[-1]), 
            deparse))
    	attr(aod, "nmodels") <- nmodels
    	class(aod) <- c("anova", "data.frame")
    	return(aod)

    } else {
    	if (!inherits(object, "Sarlm")) 
            stop("object not a fitted simultaneous autoregressive model")
        LL <- logLik(object)
        AIC <- AIC(LL)
        res <- data.frame("AIC"=AIC, "Log likelihood"=LL, "df"=attr(LL, "df"),
	    row.names=deparse(substitute(object)))
	class(res) <- c("anova", "data.frame")
        return(res)
    }
}

# Copyright 2002-12 by Roger Bivand, 2015 Martin Gubri
#

residuals.Sarlm <- function(object, ...) {
  if (is.null(object$na.action))
    object$residuals
  else napredict(object$na.action, object$residuals)
}

deviance.Sarlm <- function(object, ...) {
  object$SSE
}

coef.Sarlm <- function(object, ...) {
  ret <- NULL
  #	ret <- sqrt(object$s2)
  #	names(ret) <- "sigma"
  if (object$type == "error") ret <- c(ret, object$lambda)
  else if (object$type == "lag" || object$type == "mixed")
    ret <- c(ret, object$rho)
  else if (object$type == "sac" || object$type == "sacmixed")
    ret <- c(ret, object$rho, object$lambda)
  ret <- c(ret, object$coefficients)
  
  ret
}

vcov.Sarlm <- function(object, ...) {
  if (object$ase) res <- object$resvar[-1,-1]
  else {
    if (!is.null(object$fdHess)) {
      if (object$insert) res <- object$resvar[-1,-1]
      else res <- object$resvar
    } else {
      stop("vcov not available for this model")
    }
  }
  res
}


fitted.Sarlm <- function(object, ...) {
  message("This method assumes the response is known - see manual page")
# thanks to Philipp Otto, email 2019-11-29
  if (is.null(object$na.action))
    object$fitted.values
  else napredict(object$na.action, object$fitted.values)
}

impacts.Sarlm <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, evalues=NULL,
  useHESS=NULL, tol=1e-6, empirical=FALSE, Q=NULL) {
    if (obj$type == "error") {
        if (obj$etype == "emixed") {
            return(impactSDEM(obj))
        } else {
            stop("impact measures not for error models")
        }
    }
    if (is.null(listw) && !is.null(obj$listw_style) && 
            obj$listw_style != "W")
            stop("Only row-standardised weights supported")
    rho <- obj$rho
    beta <- obj$coefficients
    s2 <- obj$s2
    if (obj$type == "sac" || obj$type == "sacmixed") lambda <- obj$lambda
    usingHESS <- NULL
    iNsert <- obj$insert
    if (!is.null(R)) {
        resvar <- obj$resvar
        usingHESS <- FALSE
        irho <- 2
        drop2beta <- 1:2
        if (obj$type == "sac" || obj$type == "sacmixed")
            drop2beta <- c(drop2beta, 3)
        if (is.logical(resvar)) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("coefficient covariance matrix not available")
            usingHESS <- TRUE
            if (!iNsert) {
                irho <- 1
                drop2beta <- 1
                if (obj$type == "sac" || obj$type == "sacmixed")
                    drop2beta <- c(drop2beta, 2)
            }
        }
        if (!is.null(useHESS) && useHESS) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("Hessian matrix not available")
            usingHESS <- TRUE
            if (!iNsert) {
                irho <- 1
                drop2beta <- 1
                if (obj$type == "sac" || obj$type == "sacmixed")
                    drop2beta <- c(drop2beta, 2)
            }
        }
        interval <- obj$interval
        if (is.null(interval)) interval <- c(-1,0.999)
    }
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    zero_fill <- NULL
    dvars <- obj$dvars
    if (obj$type == "lag" || obj$type == "sac") {
      if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
      } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
      }
      p <- length(beta)
    } else if (obj$type == "mixed" || obj$type == "sacmixed") {
      if (!is.null(dvars)) zero_fill <- attr(dvars, "zero_fill")
      if (iicept) {
        b1 <- beta[-icept]
      } else {
        b1 <- beta
      }
      if (!is.null(zero_fill)) {
        if (length(zero_fill) > 0L) {
          inds <- attr(dvars, "inds")
          b1_long <- rep(0, 2*(dvars[1]-1))
          b1_long[1:(dvars[1]-1L)] <- b1[1:(dvars[1]-1)]
          names(b1_long)[1:(dvars[1]-1L)] <- names(b1)[1:(dvars[1]-1)]
          for (i in seq(along=inds)) {
            b1_long[(dvars[1]-1L)+(inds[i]-1L)] <- b1[(dvars[1]-1L)+i]
          }
          b1 <- b1_long
#          for (i in s_zero_fill) {
#            b1 <- append(b1, values=as.numeric(NA), after=i-1L)
#          }
        }
      }
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
    }
    n <- length(obj$residuals)
    mu <- NULL
    Sigma <- NULL
    if (!is.null(R)) {
        if (usingHESS && !iNsert) {
            mu <- c(rho, beta)
            if (obj$type == "sac" || obj$type == "sacmixed")
                mu <- c(rho, lambda, beta)
            Sigma <- fdHess
        } else {
            mu <- c(s2, rho, beta)
            if (obj$type == "sac" || obj$type == "sacmixed")
                mu <- c(s2, rho, lambda, beta)
            if (usingHESS) {
                Sigma <- fdHess
            } else {
                Sigma <- resvar
            }
        }
    }
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=interval,
        type=obj$type, tr=tr, R=R, listw=listw, evalues=evalues, tol=tol,
        empirical=empirical,Q=Q, icept=icept, iicept=iicept, p=p,
        zero_fill=zero_fill, dvars=dvars)
    attr(res, "useHESS") <- usingHESS
    attr(res, "insert") <- iNsert
    attr(res, "iClass") <- class(obj)
    res
}



