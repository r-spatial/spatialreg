# Copyright 2009-2013 by Roger Bivand


is.formula <- function(x){
   inherits(x,"formula")
}


# Copyright 1998-2011 by Roger Bivand (Wald test suggested by Rein Halbersma,
# output of correlations suggested by Michael Tiefelsdorf)
#

print.Sarlm <- function(x, ...)
{
#FIXME
       if (x$type == "error") if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
       if (x$type == "lag" || x$type == "mixed")
            if (isTRUE(all.equal(x$rho, x$interval[1])) ||
            isTRUE(all.equal(x$rho, x$interval[2]))) 
            warning("rho on interval bound - results should not be used")
	cat("\nCall:\n")
	print(x$call)
	cat("Type:", x$type, "\n")
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nLog likelihood:", logLik(x), "\n")
	invisible(x)
}

summary.Sarlm <- function(object, correlation = FALSE, Nagelkerke=FALSE,
 Hausman=FALSE, adj.se=FALSE, ...)
{
#FIXME
        if (Hausman && object$type == "error" && !is.null(object$Hcov)) {
                object$Haus <- Hausman.test(object)
        }
	if (object$type == "error") {
		object$Wald1 <- Wald1.Sarlm(object)
		if (correlation) {
                        oresvar <- object$resvar
                        ctext <- "Correlation of coefficients"
                        if (is.null(oresvar) || is.logical(oresvar) || 
                            inherits(oresvar, "try-error")) {
                            oresvar <- object$fdHess
                            ctext <- ifelse(object$insert,
                                "Approximate correlation of coefficients",
                                "** Guesswork correlation of coefficients **")
                        }
			object$correlation <- diag((diag(oresvar))
				^(-1/2)) %*% oresvar %*% 
				diag((diag(oresvar))^(-1/2))
			dimnames(object$correlation) <- dimnames(oresvar)
                        object$correltext <- ctext
		}
	} else if (object$type != "error") {
		object$Wald1 <- Wald1.Sarlm(object)
		if (correlation) {
                        oresvar <- object$resvar
                        ctext <- "Correlation of coefficients"
                        if (is.null(oresvar) || is.logical(oresvar) || 
                            inherits(oresvar, "try-error")) {
                            oresvar <- object$fdHess
                            ctext <- "Approximate correlation of coefficients"
                        }
			object$correlation <- diag((diag(oresvar))
				^(-1/2)) %*% oresvar %*% 
				diag((diag(oresvar))^(-1/2))
			dimnames(object$correlation) <- dimnames(oresvar)
                        object$correltext <- ctext
		}
        }
	object$LR1 <- LR1.Sarlm(object)

	adj <- NULL
	if (object$type == "error" || ((object$type == "lag" || 
		object$type == "mixed" || object$type == "sac" || 
                object$type == "sacmixed") && object$ase)) {
		object$coeftitle <- "(asymptotic standard errors)"
                SE <- object$rest.se
                if (adj.se) {
                    N <- length(residuals(object))
                    adj <- N/(N-(length(object$coefficients)))
                    SE <- sqrt((SE^2) * adj)
                }
#                varnames <- names(object$coefficients)
		object$Coef <- cbind(object$coefficients, SE, 
			object$coefficients/SE,
			2*(1-pnorm(abs(object$coefficients/SE))))
		colnames(object$Coef) <- c("Estimate", "Std. Error", 
			ifelse(adj.se, "t value", "z value"), "Pr(>|z|)")
	        rownames(object$Coef) <- names(object$coefficients)
	} else {
	    # intercept-only bug fix Larry Layne 20060404
            if (!is.null(object$rest.se)) {
		object$coeftitle <- "(numerical Hessian approximate standard errors)"
                SE <- object$rest.se
                if (adj.se) {
                    N <- length(residuals(object))
                    adj <- N/(N-(length(object$coefficients)))
                    SE <- sqrt((SE^2) * adj)
                }
#                varnames <- names(object$coefficients)
		object$Coef <- cbind(object$coefficients, SE, 
			object$coefficients/SE,
			2*(1-pnorm(abs(object$coefficients/SE))))
		colnames(object$Coef) <- c("Estimate", "Std. Error", 
			ifelse(adj.se, "t value", "z value"), "Pr(>|z|)")
	        rownames(object$Coef) <- names(object$coefficients)
              }
	}
# temporary fix for broom 210312
#        object$Coef <- object$coefficients
        object$adj.se <- adj

        if (Nagelkerke) {
            nk <- NK.Sarlm(object)
            if (!is.null(nk)) object$NK <- nk
        }
	structure(object, class=c("summary.Sarlm", class(object)))
}

print.summary.Sarlm <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
{
	cat("\nCall:", deparse(x$call),	sep = "", fill=TRUE)
       if (x$type == "error") if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
       if (x$type == "lag" || x$type == "mixed")
            if (isTRUE(all.equal(x$rho, x$interval[1])) ||
            isTRUE(all.equal(x$rho, x$interval[2]))) 
            warning("rho on interval bound - results should not be used")
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2L) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	cat("\nType:", x$type, "\n")
	if (x$zero.policy) {
		zero.regs <- attr(x, "zero.regs")
		if (!is.null(zero.regs))
			cat("Regions with no neighbours included:\n",
			zero.regs, "\n")
	}
        if (!is.null(x$coeftitle)) {
	    cat("Coefficients:", x$coeftitle, "\n")
	    coefs <- x$Coef
	    if (!is.null(aliased <- x$aliased) && any(x$aliased)){
		cat("    (", table(aliased)["TRUE"], 
			" not defined because of singularities)\n", sep = "")
		cn <- names(aliased)
		coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                	colnames(x$Coef)))
            	coefs[!aliased, ] <- x$Coef
	    }
	    printCoefmat(coefs, signif.stars=signif.stars, digits=digits,
		na.print="NA")
	}

	res <- x$LR1
        pref <- ifelse(x$ase, "Asymptotic", "Approximate (numerical Hessian)")
	if (x$type == "error") {
		cat("\nLambda: ", format(signif(x$lambda, digits)),
			", LR test value: ", format(signif(res$statistic,
                        digits)), ", p-value: ", format.pval(res$p.value,
                        digits), "\n", sep="")
		if (!is.null(x$lambda.se)) {
                    if (!is.null(x$adj.se)) {
                        x$lambda.se <- sqrt((x$lambda.se^2)*x$adj.se)   
                    }
		    cat(pref, " standard error: ", 
		        format(signif(x$lambda.se, digits)),
			ifelse(is.null(x$adj.se), "\n    z-value: ",
                               "\n    t-value: "), format(signif((x$lambda/
				x$lambda.se), digits)),
			", p-value: ", format.pval(2*(1-pnorm(abs(x$lambda/
				x$lambda.se))), digits), "\n", sep="")
		    cat("Wald statistic: ", format(signif(x$Wald1$statistic, 
			digits)), ", p-value: ", format.pval(x$Wald1$p.value, 
			digits), "\n", sep="")
		}
	} else if (x$type == "sac" || x$type == "sacmixed") {
		cat("\nRho: ", format(signif(x$rho, digits)), "\n",
                    sep="")
                if (!is.null(x$rho.se)) {
                    if (!is.null(x$adj.se)) {
                        x$rho.se <- sqrt((x$rho.se^2)*x$adj.se)   
                    }
		  cat(pref, " standard error: ", 
			format(signif(x$rho.se, digits)), 
                        ifelse(is.null(x$adj.se), "\n    z-value: ",
                               "\n    t-value: "), 
			format(signif((x$rho/x$rho.se), digits)),
			", p-value: ", format.pval(2 * (1 - pnorm(abs(x$rho/
				x$rho.se))), digits), "\n", sep="")
                }
		cat("Lambda: ", format(signif(x$lambda, digits)), "\n", sep="")
		if (!is.null(x$lambda.se)) {
                    pref <- ifelse(x$ase, "Asymptotic",
                        "Approximate (numerical Hessian)")
                    if (!is.null(x$adj.se)) {
                        x$lambda.se <- sqrt((x$lambda.se^2)*x$adj.se)   
                    }
		    cat(pref, " standard error: ", 
		        format(signif(x$lambda.se, digits)),
			ifelse(is.null(x$adj.se), "\n    z-value: ",
                               "\n    t-value: "), format(signif((x$lambda/
				x$lambda.se), digits)),
			", p-value: ", format.pval(2*(1-pnorm(abs(x$lambda/
				x$lambda.se))), digits), "\n", sep="")
                }
                cat("\nLR test value: ", format(signif(res$statistic, digits)),
		    ", p-value: ", format.pval(res$p.value, digits), "\n",
                    sep="")
        } else {
		cat("\nRho: ", format(signif(x$rho, digits)), 
                    ", LR test value: ", format(signif(res$statistic, digits)),
		    ", p-value: ", format.pval(res$p.value, digits), "\n",
                    sep="")
                if (!is.null(x$rho.se)) {
                    if (!is.null(x$adj.se)) {
                        x$rho.se <- sqrt((x$rho.se^2)*x$adj.se)   
                    }
		  cat(pref, " standard error: ", 
			format(signif(x$rho.se, digits)), 
                        ifelse(is.null(x$adj.se), "\n    z-value: ",
                               "\n    t-value: "), 
			format(signif((x$rho/x$rho.se), digits)),
			", p-value: ", format.pval(2 * (1 - pnorm(abs(x$rho/
				x$rho.se))), digits), "\n", sep="")
                }
		if (!is.null(x$Wald1)) {
		    cat("Wald statistic: ", format(signif(x$Wald1$statistic, 
			digits)), ", p-value: ", format.pval(x$Wald1$p.value, 
			digits), "\n", sep="")
		}

	}
	cat("\nLog likelihood:", logLik(x), "for", x$type, "model\n")
	cat("ML residual variance (sigma squared): ", 
		format(signif(x$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$s2), digits)), ")\n", sep="")
        if (!is.null(x$NK)) cat("Nagelkerke pseudo-R-squared:",
            format(signif(x$NK, digits)), "\n")
	cat("Number of observations:", length(x$residuals), "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	cat("AIC: ", format(signif(AIC(x), digits)), ", (AIC for lm: ",
		format(signif(x$AIC_lm.model, digits)), ")\n", sep="")
	if (x$type == "error") {
		if (!is.null(x$Haus)) {
		    cat("Hausman test: ", format(signif(x$Haus$statistic, 
			digits)), ", df: ", format(x$Haus$parameter),
                        ", p-value: ", format.pval(x$Haus$p.value, digits),
                        "\n", sep="")
		}
        }
	if ((x$type == "lag" || x$type ==  "mixed") && x$ase) {
		cat("LM test for residual autocorrelation\n")
		cat("test value: ", format(signif(x$LMtest, digits)),
			", p-value: ", format.pval((1 - pchisq(x$LMtest, 1)), 
			digits), "\n", sep="")
	}
        if (x$type != "error" && !is.null(x$LLCoef)) {
		cat("\nCoefficients: (log likelihood/likelihood ratio)\n")
		printCoefmat(x$LLCoef, signif.stars=signif.stars,
			digits=digits, na.print="NA")
        }
    	correl <- x$correlation
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\n", x$correltext, "\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)
}

coef.summary.Sarlm <- function(object, ...) object$Coef

getVmate <- function(coefs, env, s2, trs, tol.solve=1.0e-10, optim=FALSE,
    optimM="optimHess") {
    if (optim) {
      if (optimM == "nlm") {
           options(warn=-1)
           opt <- nlm(f=f_laglm_hess_nlm, p=coefs, env=env, hessian=TRUE)
           options(warn=0)
           mat <- opt$hessian
#        opt <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
#        mat <- opt
       } else if (optimM == "optimHess") {
           mat <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
       } else {
           opt <- optim(par=coefs, fn=f_laglm_hess, env=env, method=optimM,
           hessian=TRUE)
           mat <- opt$hessian
      }
#        opt <- optimHess(par=coefs, fn=f_errlm_hess, env=env)
#        mat <- opt
    } else {
        fd <- fdHess(coefs, f_errlm_hess, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asye(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

sar_error_hess_sse <- function(lambda, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml1_sse_env", env, lambda, beta, PACKAGE="spatialreg")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
        res <- get("sw", envir=env) * (yl - (xl %*% beta))
        SSE <- c(crossprod(res))
    }
    SSE
}

f_errlm_hess <- function(coefs, env) {
    lambda <- coefs[1]
    int <- get("interval", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-1]
    SSE <- sar_error_hess_sse(lambda, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    det <- do_ldet(lambda, env)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2) * log(2 * pi)) - 
        (n/2) * log(s2) - (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
        " Jacobian:", det, " SSE:", SSE, "\n")
    assign("hf_calls", get("hf_calls", envir=env)+1L, envir=env)
    if (!is.finite(ret)) return(-Inf)
    ret
}

insert_asye <- function(coefs, env, s2, mat, trs) {
    lambda <- coefs[1]
    p <- length(coefs)-1L
    p2 <- p+2
    omat <- matrix(0, nrow=p2, ncol=p2)
    LX <- get("sw", envir=env) * (get("x", envir=env) - lambda *
        get("WX", envir=env))
#    omat[3:p2, 3:p2] <- -crossprod(LX)*s2
#    omat[3:p2, 3:p2] <- -crossprod(LX)
    omat[3:p2, 3:p2] <- -crossprod(LX)/s2
    omat[2, 2] <- mat[1, 1]
    n <- get("n", envir=env)
    omat[1, 1] <- -n/(2*(s2^2))
#    omat[1, 1] <- -n/(2*(s2))
    omat[1, 2] <- omat[2, 1] <- -trB(lambda, trs)/s2
#    omat[1, 2] <- omat[2, 1] <- -trB(lambda, trs)
    omat
}


getVmatl <- function(coefs, env, s2, trs, tol.solve=1.0e-10, optim=FALSE,
    optimM="optimHess") {
    if (optim) {
      if (optimM == "nlm") {
           options(warn=-1)
           opt <- nlm(f=f_laglm_hess_nlm, p=coefs, env=env, hessian=TRUE)
           options(warn=0)
           mat <- opt$hessian
#        opt <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
#        mat <- opt
       } else if (optimM == "optimHess") {
           mat <- optimHess(par=coefs, fn=f_laglm_hess, env=env)
       } else {
           opt <- optim(par=coefs, fn=f_laglm_hess, env=env, method=optimM,
           hessian=TRUE)
           mat <- opt$hessian
      }
    } else {
        fd <- fdHess(coefs, f_laglm_hess, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}


sar_lag_hess_sse <- function(rho, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spatialreg")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        res <- (get("y", envir=env) - rho * get("wy", envir=env)) - 
            get("x", envir=env) %*% beta
        SSE <- c(crossprod(res))
    }
    SSE
}

f_laglm_hess <- function(coefs, env) {
    rho <- coefs[1]
    int <- get("interval", envir=env)
    if (rho <= int[1] || rho >= int[2]) return(-Inf)
    beta <- coefs[-1]
    SSE <- sar_lag_hess_sse(rho, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    det <- do_ldet(rho, env)
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("Hessian: rho:\t", rho, "\tfunction value:\t", ret, "\n")
    assign("hf_calls", get("hf_calls", envir=env)+1L, envir=env)
    if (!is.finite(ret)) return(-Inf)
    ret
}

f_laglm_hess_nlm <- function(coefs, env) {
    ret <- f_laglm_hess(coefs, env)
    -ret
}


trB <- function(rho, tr)  sum(sapply(0:(length(tr)-1L),
    function(i) rho^i * tr[i+1]))

insert_asy <- function(coefs, env, s2, mat, trs) {
    p <- length(coefs)-1L
    p2 <- p+2
    n <- get("n", envir=env)
    omat <- matrix(0, nrow=p2, ncol=p2)
    omat[3:p2, 3:p2] <- -crossprod(get("x", envir=env))/s2
    omat[2, 2] <- mat[1, 1]
    omat[2, 3:p2] <- omat[3:p2, 2] <- -c(crossprod(get("wy", envir=env),
        get("x", envir=env))/s2)
    omat[1, 1] <- -n/(2*(s2^2))
    omat[1, 2] <- omat[2, 1] <- -trB(coefs[1], trs)/s2
    omat
}
sar_sac_hess_sse <- function(rho, lambda, beta, env) {
    yl <- get("y", envir=env) - rho * get("wy", envir=env) - 
        lambda * get("w2y", envir=env) + rho * lambda * get("w2wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
    SSE <- c(crossprod(res))
    SSE
}


getVmatsac <- function(coefs, env, tol.solve=1.0e-10) {
    fd <- fdHess(coefs, f_sac_hess, env)
    mat <- fd$Hessian
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_sac_hess <- function(coefs, env) {
    rho <- coefs[1]
    int <- get("interval1", envir=env)
    if (rho <= int[1] || rho >= int[2]) return(-Inf)
    lambda <- coefs[2]
    int <- get("interval2", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-(1:2)]
    SSE <- sar_sac_hess_sse(rho, lambda, beta, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet1 <- do_ldet(rho, env, which=1)
    ldet2 <- do_ldet(lambda, env, which=2)
    ret <- (ldet1 + ldet2 - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("rho:", rho, "lambda:", lambda,
        " function:", ret, " Jacobian1:", ldet1, " Jacobian2:",
        ldet2, " SSE:", SSE, "\n")
    if (!is.finite(ret)) return(-Inf)
   ret
}


sar_sac_hess_sse <- function(rho, lambda, beta, env) {
    yl <- get("y", envir=env) - rho * get("wy", envir=env) - 
        lambda * get("w2y", envir=env) + rho * lambda * get("w2wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
    SSE <- c(crossprod(res))
    SSE
}


