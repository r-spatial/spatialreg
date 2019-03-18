logLik.sarlm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- length(residuals(object))
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

NK.sarlm <- function(obj) {
     n <- length(residuals(obj))
     nullLL <- obj$LLNullLlm
     if (is.null(nullLL)) return(nullLL)
     c(1 - exp(-(2/n)*(logLik(obj) - nullLL)))
}


LR.sarlm <- function(x, y)
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


LR1.sarlm <- function(object)
{
	if (!inherits(object, "sarlm")) stop("Not a sarlm object")
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

Wald1.sarlm <- function(object) {
	if (!inherits(object, "sarlm")) stop("Not a sarlm object")
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

Hausman.test.sarlm <- function(object, ..., tol=NULL) {
    if (!inherits(object, "sarlm")) stop("not a sarlm object")
    if (object$type != "error") stop("not a spatial error model")
    fmeth <- ifelse(object$method != "eigen", "(approximate)", "(asymptotic)") 
    if (is.null(object$Hcov)) stop("Vo not available")
    s2 <- object$s2
    Vo <- s2 * object$Hcov
    Vs <- s2 * object$Vs
    d <- object$coef_lm.model - object$coefficients
    if (!is.null(tol)) VV <- try(solve((Vo - Vs), tol=tol))
    else VV <- try(solve(Vo - Vs))
    if (class(VV) == "try.error") {
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

