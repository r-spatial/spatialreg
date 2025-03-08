# Copyright 1998-2016 by Roger Bivand (non-W styles Rein Halbersma)
#

errorsarlm <- function(formula, data = list(), listw, na.action, weights=NULL,
        Durbin, etype, method="eigen", quiet=NULL, zero.policy=NULL,
        interval=NULL, tol.solve=.Machine$double.eps, trs=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5, returnHcov=TRUE,
            pWOrder=250, fdHess=NULL, optimHess=FALSE,
            optimHessMethod="optimHess", LAPACK=FALSE,
            compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
            super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
            correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
            interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
            LU_order=FALSE, pre_eig=NULL, return_impacts=TRUE)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spatialregOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
        stopifnot(is.logical(zero.policy))
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        stopifnot(is.logical(con$optimHess))
        stopifnot(is.logical(con$LAPACK))
#        stopifnot(is.logical(con$super))
        stopifnot(is.logical(con$compiled_sse))
        stopifnot(is.character(con$spamPivot))
        stopifnot(is.logical(con$return_impacts))
        if (!inherits(formula, "formula")) formula <- as.formula(formula)
#	mt <- terms(formula, data = data)
#	mf <- lm(formula, data, na.action=na.action, method="model.frame")
#
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")
        if (attr(mt, "intercept") == 1 && !any(attr(mt, "factors") == 1) &&
            (!missing(Durbin)) && (is.formula(Durbin) || isTRUE(Durbin))) {
            warning("intercept-only model, Durbin invalid and set FALSE")
            Durbin <- FALSE
        }
        have_factor_preds <- have_factor_preds_mf(mf)
#
	na.act <- attr(mf, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
            if (!is.null(con$pre_eig)) {
                warning("NAs found, precomputed eigenvalues ignored")
                con$pre_eig <- NULL
            }
	}
        if (missing(etype)) etype <- "error"
        if (etype == "Durbin") etype <- "emixed"
        if (missing(Durbin)) Durbin <- ifelse(etype == "error", FALSE, TRUE)
# FIXME does this hold?
#        if (listw$style != "W" && is.formula(Durbin)) {
#            Durbin <- TRUE
#            warning("formula Durbin requires row-standardised weights; set TRUE")
#        }
        if (is.logical(Durbin) && isTRUE(Durbin)) {
            etype <- "emixed"
            if (have_factor_preds) warn_factor_preds(have_factor_preds)
        }
        if (is.formula(Durbin)) etype <- "emixed"
        if (is.logical(Durbin) && !isTRUE(Durbin)) etype <- "error"
        
	switch(etype, error = if (!quiet)
                cat("\nSpatial autoregressive error model\n"),
	    emixed = if (!quiet)
                cat("\nSpatial mixed autoregressive error model\n"),
	    stop("\nUnknown model type\n"))

	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
        n <- nrow(x)
	if (n != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
#
        weights <- as.vector(model.extract(mf, "weights"))
        if (!is.null(weights) && !is.numeric(weights)) 
            stop("'weights' must be a numeric vector")
        if (is.null(weights)) weights <- rep(as.numeric(1), n)
        if (any(is.na(weights))) stop("NAs in weights")
        if (any(weights < 0)) stop("negative weights")
#
	m <- NCOL(x)
	xxcolnames <- colnames(x)

        stopifnot(is.logical(con$small_asy))
        if (method != "eigen") {
            if (con$small >= n && con$small_asy) do_asy <- TRUE
            else do_asy <- FALSE
        } else do_asy <- TRUE
        if (is.null(con$fdHess)) {
            con$fdHess <- method != "eigen" && !do_asy
            fdHess <- NULL
        }
        stopifnot(is.logical(con$fdHess))
        dvars <- c(NCOL(x), 0L)

	if (is.formula(Durbin) || isTRUE(Durbin)) {
                prefix <- "lag"
                if (isTRUE(Durbin)) {
                    WX <- create_WX(x, listw, zero.policy=zero.policy,
                        prefix=prefix)
                } else {
                    data1 <- data
                    if (!is.null(na.act) && (inherits(na.act, "omit") ||
                        inherits(na.act, "exclude"))) {
                        data1 <- data1[-c(na.act),]
                    }
	            dmf <- lm(Durbin, data1, na.action=na.fail, 
		        method="model.frame")
	           formula_durbin_factors <- have_factor_preds_mf(dmf)
                   if (formula_durbin_factors) 
                       warn_factor_preds(formula_durbin_factors)
#	            dmf <- lm(Durbin, data, na.action=na.action, 
#		        method="model.frame")
                    fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
                    if (inherits(fx, "try-error")) 
                        stop("Durbin variable mis-match")
                    WX <- create_WX(fx, listw, zero.policy=zero.policy,
                        prefix=prefix)
                    inds <- match(substring(colnames(WX), 5,
	                nchar(colnames(WX))), colnames(x))
                    if (anyNA(inds)) stop("WX variables not in X: ",
                        paste(substring(colnames(WX), 5,
                        nchar(colnames(WX)))[is.na(inds)], collapse=" "))
                    icept <- grep("(Intercept)", colnames(x))
                    iicept <- length(icept) > 0L
                    if (iicept) {
                        xn <- colnames(x)[-1]
                    } else {
                        xn <- colnames(x)
                    }
                    wxn <- substring(colnames(WX), nchar(prefix)+2,
                        nchar(colnames(WX)))
                    zero_fill <- NULL
                    if (length((which(!(xn %in% wxn)))) > 0L) 
                        zero_fill <- length(xn) + (which(!(xn %in% wxn)))
                }
                dvars <- c(NCOL(x), NCOL(WX))
                if (is.formula(Durbin)) {
                    attr(dvars, "f") <- Durbin
                    attr(dvars, "inds") <- inds
                    attr(dvars, "zero_fill") <- zero_fill
                }
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1, weights=weights)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
	    x <- x[,-nacoef]
	    m <- NCOL(x)
	    xxcolnames <- colnames(x)
          }
	}
#
        sw <- sqrt(weights)
#
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1,
            weights=weights))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
# added no intercept Guillaume Blanchet 091103
	if (m > 1 || (m == 1 && K == 1)) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)

	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
        sum_lw <- sum(log(weights))

#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("family", "SAR", envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("first_time", TRUE, envir=env)
        assign("LAPACK", con$LAPACK, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("listw", listw, envir=env)
        assign("similar", FALSE, envir=env)
        assign("f_calls", 0L, envir=env)
        assign("hf_calls", 0L, envir=env)
        assign("sum_lw", sum_lw, envir=env)
        assign("sw", sw, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (!quiet) cat(paste("\nJacobian calculated using "))

        interval <- jacobianSetup(method, env, con, pre_eig=con$pre_eig,
            trs=trs, interval=interval)
        assign("interval", interval, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (con$compiled_sse) {
             ptr <- .Call("opt_error_init", PACKAGE="spatialreg")
             assign("ptr", ptr, envir=env)
        }
	opt <- optimize(sar.error.f, interval=interval, 
		maximum=TRUE, tol=con$tol.opt, env=env)
	lambda <- opt$maximum
        if (isTRUE(all.equal(lambda, interval[1])) ||
            isTRUE(all.equal(lambda, interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	names(lambda) <- "lambda"
	LL <- opt$objective
        if (con$compiled_sse) {
             .Call("opt_error_free", get("ptr", envir=env), PACKAGE="spatialreg")
        }
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1,
            weights=weights)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
        sum_lm_target <- summary.lm(lm.target, correlation = FALSE)
        emixedImps <- NULL
        if (any(sum_lm_target$aliased)) warning("aliased variables found")
	if (con$return_impacts && etype == "emixed") {
          if (isTRUE(Durbin)) {
            m.1 <- m > 1
            if (m.1 && K == 2) {
                m2 <- (m-1)/2
            } else {
                m2 <- m/2
            }
                cm <- matrix(0, ncol=m, nrow=m2)
                if (K == 2) {
                    if (m.1) {
                        rownames(cm) <- xxcolnames[2:(m2+1)]
                    } else {
                        rownames(cm) <- xxcolnames[1:m2]
                    }
                    LI <- ifelse(!(listw$style == "W"), 1, 0) #TR: lagged intercept
                    for (i in 1:m2) cm[i, c(i+1, i+(m2+1 + LI)) ] <- 1 #TR: Add to index
# drop bug fix 2016-09-21 Philipp Hunziker
                    dirImps <- sum_lm_target$coefficients[2:(m2+1), 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_lm_target$coefficients[((m2 + 2):m + LI), 1:2, drop=FALSE] #TR: Add to index
                    rownames(indirImps) <- rownames(cm)
                } else {
                    rownames(cm) <- xxcolnames[1:m2]
                    for (i in 1:m2) cm[i, c(i, i+m2)] <- 1
                    dirImps <- sum_lm_target$coefficients[1:m2, 1:2, drop=FALSE]
                    rownames(dirImps) <- rownames(cm)
                    indirImps <- sum_lm_target$coefficients[(m2+1):m, 1:2, drop=FALSE]
                    rownames(indirImps) <- rownames(cm)
                }
            suppressWarnings(lc <- summary(multcomp::glht(lm.target, linfct=cm)))
            totImps <- cbind("Estimate"=lc$test$coefficients, "Std. Error"=lc$test$sigma)
          } else if (is.formula(Durbin)) {
#FIXME
            LI <- ifelse(listw$style != "W" 
                         && attr(terms(Durbin), "intercept") == 1, 1, 0)
              m <- sum(dvars)
              KIL <- max((LI - (K - 1)), 0)
              m2 <- dvars[2] - KIL
              cm <- matrix(0, ncol=m, nrow=m2)
              for (i in 1:m2) {
                  cm[i, c(inds[i], i+dvars[1] + KIL)] <- 1
              }
              if (LI == 1 && K == 1) { 
                rownames(cm) <- wxn[!grepl("Intercept", wxn)]
              } else {
                rownames(cm) <- wxn
              }
              dirImps <- sum_lm_target$coefficients[K:dvars[1], 1:2,
                drop=FALSE]
              rownames(dirImps) <- xn
              indirImps <- sum_lm_target$coefficients[(dvars[1] + 1 + KIL):m, 1:2,
                drop=FALSE]
              if (!is.null(zero_fill)) {
                if (length(zero_fill) > 0L) {
                 lres <- vector(mode="list", length=2L)
                 for (j in 1:2) {
                  jindirImps <- rep(as.numeric(NA), (dvars[1] + (1-K)))
                  for (i in seq(along=inds)) {
                    jindirImps[(inds[i]+(1-K))] <- indirImps[i, j]
                  }
                  lres[[j]] <- jindirImps
                 }
                 indirImps <- do.call("cbind", lres)
                }
              }
              rownames(indirImps) <- xn
              suppressWarnings(lc <- summary(multcomp::glht(lm.target, linfct=cm)))
              totImps <- cbind("Estimate"=lc$test$coefficients, "Std. Error"=lc$test$sigma)
              if (!is.null(zero_fill)) {
                if (length(zero_fill) > 0L) {
                 lres <- vector(mode="list", length=2L)
                 for (j in 1:2) {
                  jtotImps <- dirImps[, j]
                  for (i in seq(along=inds)) {
                    jtotImps[(inds[i]+(1-K))] <- totImps[i, j]
                  }
                  lres[[j]] <- jtotImps
                 }
                 totImps <- do.call("cbind", lres)
                }
              }
              rownames(totImps) <- xn
          } else stop("undefined emixed state")
          emixedImps <- list(dirImps=dirImps, indirImps=indirImps,
              totImps=totImps)
        }
        Vs <- sum_lm_target$cov.unscaled
        tarX <- model.matrix(lm.target)
        tary <- model.response(model.frame(lm.target))
	lm.model <- lm(y ~ x - 1, weights=weights)
        logLik_lm.model <- logLik(lm.model)
        AIC_lm.model <- AIC(lm.model)
	ase <- FALSE
	lambda.se <- NULL
	LMtest <- NULL
	asyvar1 <- FALSE
        force_assign_eigen <- FALSE
        Hcov <- NULL
        pWinternal <- NULL
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
	if (method == "eigen" || do_asy) {
		tr <- function(A) sum(diag(A))
		W <- listw2mat(listw)
		A <- solve(diag(n) - lambda*W)
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, (sw*x)) %*% A
#                    A <- solve(diag(n) - lambda*t(W))
                    C <- A %*% (sw*x) %*% R
                    Hcov <- B %*% C
                    attr(Hcov, "method") <- method
                    timings[["eigen_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }

		WA <- W %*% A
		asyvar <- matrix(0, nrow=2+p, ncol=2+p)
		asyvar[1,1] <- n / (2*(s2^2))
#		asyvar[1,1] <- n / (2*(s2))
		asyvar[2,1] <- asyvar[1,2] <- tr(WA) / s2
#		asyvar[2,1] <- asyvar[1,2] <- tr(WA)
		asyvar[2,2] <- tr(WA %*% WA) + tr(crossprod(WA))
# bug found 100224 German Muchnik Izon
#		asyvar[3:(p+2),3:(p+2)] <- s2*(t(x - lambda*WX) %*% 
                xl <- sw * (x - lambda*WX)
#		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)
		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)/s2
		asyvar1 <- try(solve(asyvar, tol=tol.solve), silent=TRUE)
                if (inherits(asyvar1, "try-error")) {
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                    con$fdHess <- TRUE
                    force_assign_eigen <- TRUE
                    warning(paste("inversion of asymptotic covariance",
                        "matrix failed for tol.solve =", tol.solve,
                        "\n", strsplit(attr(asyvar1, "condition")$message,
                            ":")[[1]][2], "- using numerical Hessian."))
                } else {
		    rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "lambda", xcolnames)
		
		    lambda.se <- sqrt(asyvar1[2,2])
#		lambda.se <- sqrt(s2*asyvar1[2,2])
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
		    ase <- TRUE
                }
	} else {
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, (sw*x))
                    W <- as(get("listw", envir=env), "CsparseMatrix")
                    B0 <- powerWeights(W=W, rho=lambda, order=con$pWOrder,
                        X=B, tol=tol.solve)
                    if (!is.null(attr(B0, "internal")) &&
                        !attr(B0, "internal")$conv)
                        pWinternal <- c(pWinternal, attr(B0, "internal"))
                    B1 <- as(B0, "matrix")
                    C <- (sw*x) %*% R
                    C0 <- powerWeights(W=t(W), rho=lambda, order=con$pWOrder,
                        X=C, tol=tol.solve)
                    if (!is.null(attr(C0, "internal")) &&
                        !attr(C0, "internal")$conv)
                        pWinternal <- c(pWinternal, attr(C0, "internal"))
                    C1 <- as(C0, "matrix")
                    Hcov <- B1 %*% C1
                    attr(Hcov, "method") <- method
                    timings[["sparse_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
        }
        if (con$fdHess) {
            coefs <- c(lambda, coef.lambda)
            if (con$compiled_sse) {
                ptr <- .Call("hess_error_init", PACKAGE="spatialreg")
                assign("ptr", ptr, envir=env)
            }
            fdHess <- getVmate(coefs, env, s2, trs, tol.solve=tol.solve,
                optim=con$optimHess, optimM=con$optimHessMethod)
            if (con$compiled_sse) {
                .Call("hess_error_free", get("ptr", envir=env),
                    PACKAGE="spatialreg")
            }
            if (is.null(trs)) {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("lambda", colnames(x))
                if (method != "eigen" || force_assign_eigen) {
                    lambda.se <- sqrt(fdHess[1, 1])
                }
            } else {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("sigma2", "lambda", colnames(x))
                if (method != "eigen" || force_assign_eigen) {
                    lambda.se <- sqrt(fdHess[2, 2])
                }
            }
            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
        } else fdHess <- FALSE
	call <- match.call()
        if (method=="SE_classic") {
            iC <- get("intern_classic", envir=env)
        } else iC <- NULL
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type="error", etype=etype, lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), #lm.model=lm.model, 
                logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
                coef_lm.model=coef(lm.model),
                tarX=tarX, tary=tary, y=y, X=x,
		method=method, call=call, residuals=r, #lm.target=lm.target,
		opt=opt, fitted.values=fit, ase=ase, #formula=formula,
		se.fit=NULL, resvar=asyvar1, similar=get("similar", envir=env),
		lambda.se=lambda.se, LMtest=LMtest, zero.policy=zero.policy, 
		aliased=aliased, LLNullLlm=LL_null_lm, Hcov=Hcov, Vs=Vs,
                interval=interval, fdHess=fdHess,
                optimHess=con$optimHess, insert=!is.null(trs), trs=trs,
                timings=do.call("rbind", timings)[, c(1, 3)], 
                f_calls=get("f_calls", envir=env),
                hf_calls=get("hf_calls", envir=env), intern_classic=iC,
                pWinternal=pWinternal, weights=weights, emixedImps=emixedImps,
                dvars=dvars), class=c("Sarlm"))
        rm(env)
        GC <- gc()
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0L)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
        if (is.formula(Durbin)) attr(ret, "Durbin") <- deparse(Durbin)
	ret
}

sar_error_sse <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spatialreg")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        yl <- get("sw", envir=env) * yl
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
        xl <- get("sw", envir=env) * xl
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
#        xl.q.yl <- qr.qty(qr(xl, LAPACK=get("LAPACK", envir=env)), yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    SSE
}


sar.error.f <- function(lambda, env) {
    SSE <- sar_error_sse(lambda, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    ret <- (ldet + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
    assign("f_calls", get("f_calls", envir=env)+1L, envir=env)
    ret
}


lagsarlm <- function(formula, data = list(), listw, 
	na.action, Durbin, type, method="eigen", quiet=NULL, 
	zero.policy=NULL, interval=NULL, tol.solve=.Machine$double.eps, 
	trs=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5,
            fdHess=NULL, optimHess=FALSE, optimHessMethod="optimHess",
            compiled_sse=FALSE, Imult=2,
            cheb_q=5, MC_p=16L, MC_m=30L, super=NULL, spamPivot="MMD",
            in_coef=0.1, type="MC", correct=TRUE, trunc=TRUE,
            SE_method="LU", nrho=200, interpn=2000, small_asy=TRUE,
            small=1500, SElndet=NULL, LU_order=FALSE, pre_eig=NULL,
            OrdVsign=1)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spatialregOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
        if (!inherits(formula, "formula")) formula <- as.formula(formula)
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, 
		method="model.frame")
        have_factor_preds <- have_factor_preds_mf(mf)
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
            if (!is.null(con$pre_eig)) {
                warning("NAs found, precomputed eigenvalues ignored")
                con$pre_eig <- NULL
            }
	}

        if (missing(type)) type <- "lag"
        if (type == "Durbin") type <- "mixed"
        if (missing(Durbin)) Durbin <- ifelse(type == "lag", FALSE, TRUE)
        if (listw$style != "W" && is.formula(Durbin)) {
            Durbin <- TRUE
            warning("formula Durbin requires row-standardised weights; set TRUE")
        }
        if (is.logical(Durbin) && isTRUE(Durbin)) {
            type <- "mixed"
            if (have_factor_preds) warn_factor_preds(have_factor_preds)
        }
        if (is.formula(Durbin)) type <- "mixed"
        if (is.logical(Durbin) && !isTRUE(Durbin)) type <- "lag"
	switch(type, lag = if (!quiet) cat("\nSpatial lag model\n"),
	    mixed = if (!quiet) cat("\nSpatial mixed autoregressive model\n"),
	    stop("\nUnknown model type\n"))
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)
	if (NROW(x) != length(listw$neighbours))
		stop("Input data and weights have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
        stopifnot(is.logical(con$small_asy))
        if (method != "eigen") {
            if (con$small >= n && con$small_asy) do_asy <- TRUE
            else do_asy <- FALSE
        } else do_asy <- TRUE
        if (is.null(con$fdHess)) {
            con$fdHess <- method != "eigen" && !do_asy
            fdHess <- NULL
        }
        stopifnot(is.logical(con$fdHess))
        stopifnot(is.numeric(con$OrdVsign))
        stopifnot(length(con$OrdVsign) == 1)
        stopifnot(abs(con$OrdVsign) == 1)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy))) stop("NAs in lagged dependent variable")
        dvars <- c(NCOL(x), 0L)
#FIXME
	if (is.formula(Durbin) || isTRUE(Durbin)) {
                prefix <- "lag"
                if (isTRUE(Durbin)) {
                    WX <- create_WX(x, listw, zero.policy=zero.policy,
                        prefix=prefix)
                } else {
                    data1 <- data
                    if (!is.null(na.act) && (inherits(na.act, "omit") ||
                        inherits(na.act, "exclude"))) {
                        data1 <- data1[-c(na.act),]
                    }
	            dmf <- lm(Durbin, data1, na.action=na.fail, 
		        method="model.frame")
	            formula_durbin_factors <- have_factor_preds_mf(dmf)
                    if (formula_durbin_factors) 
                        warn_factor_preds(formula_durbin_factors)
                    fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
                    if (inherits(fx, "try-error")) 
                        stop("Durbin variable mis-match")
                    WX <- create_WX(fx, listw, zero.policy=zero.policy,
                        prefix=prefix)
                    inds <- match(substring(colnames(WX), 5,
	                nchar(colnames(WX))), colnames(x))
                    if (anyNA(inds)) stop("WX variables not in X: ",
                        paste(substring(colnames(WX), 5,
                        nchar(colnames(WX)))[is.na(inds)], collapse=" "))
                    icept <- grep("(Intercept)", colnames(x))
                    iicept <- length(icept) > 0L
                    if (iicept) {
                        xn <- colnames(x)[-1]
                    } else {
                        xn <- colnames(x)
                    }
                    wxn <- substring(colnames(WX), nchar(prefix)+2,
                        nchar(colnames(WX)))
                    zero_fill <- NULL
                    if (length((which(!(xn %in% wxn)))) > 0L)
                        zero_fill <- length(xn) + (which(!(xn %in% wxn)))
                }
                dvars <- c(NCOL(x), NCOL(WX))
                if (is.formula(Durbin)) {
                    attr(dvars, "f") <- Durbin
                    attr(dvars, "inds") <- inds
                    attr(dvars, "zero_fill") <- zero_fill
                }
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
#FIXME
	if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
		x <- x[,-nacoef]
          }
	}
	LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	similar <- FALSE
	lm.null <- lm(y ~ x - 1)
        logLik_lm.model <- logLik(lm.null)
        AIC_lm.model <- AIC(lm.null)
	lm.w <- lm.fit(x, wy)
	e.null <- lm.null$residuals
	e.w <- lm.w$residuals
	e.a <- t(e.null) %*% e.null
	e.b <- t(e.w) %*% e.null
	e.c <- t(e.w) %*% e.w
#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("wy", wy, envir=env)
        assign("x", x, envir=env)
        assign("n", n, envir=env)
        assign("m", m, envir=env)
        assign("K", K, envir=env)
        assign("e.a", e.a, envir=env)
        assign("e.b", e.b, envir=env)
        assign("e.c", e.c, envir=env)
        assign("family", "SAR", envir=env)
        assign("verbose", !quiet, envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("listw", listw, envir=env)
        assign("similar", FALSE, envir=env)
        assign("f_calls", 0L, envir=env)
        assign("hf_calls", 0L, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	if (!quiet) cat("Jacobian calculated using ")

        interval <- jacobianSetup(method, env, con, pre_eig=con$pre_eig,
            trs=trs, interval=interval)
        assign("interval", interval, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	opt <- optimize(sar.lag.mixed.f, interval=interval, 
		maximum=TRUE, tol=con$tol.opt, env=env)
	rho <- opt$maximum
        if (isTRUE(all.equal(rho, interval[1])) ||
            isTRUE(all.equal(rho, interval[2]))) 
            warning("rho on interval bound - results should not be used")
	names(rho) <- "rho"
	LL <- opt$objective
	optres <- opt
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	lm.lag <- lm((y - rho*wy) ~ x - 1)
	r <- residuals(lm.lag)
	fit <- y - r
	names(r) <- names(fit)
	coef.rho <- coefficients(lm.lag)
        tarX <- model.matrix(lm.lag)
        tary <- model.response(model.frame(lm.lag))
	names(coef.rho) <- colnames(x)
	SSE <- deviance(lm.lag)
	s2 <- SSE/n
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
        LMtest <- NULL
	varb <- FALSE
	ase <- FALSE
        force_assign_eigen <- FALSE
	if (method == "eigen" || do_asy) {
		rest.se <- NULL
		rho.se <- NULL
		tr <- function(A) sum(diag(A))
# beware of complex eigenvalues!
                if (do_asy && method != "eigen") eigen_setup(env)
                eig <- get("eig", envir=env)
		O <- (eig/(1-rho*eig))^2
		omega <- sum(O)
		if (is.complex(omega)) omega <- Re(omega)
		W <- listw2mat(listw)
		A <- solve(diag(n) - rho*W)
		AW <- A %*% W
		zero <- rbind(rep(0,length(coef.rho)))
		xtawxb <- s2*(t(x) %*% AW %*% x %*% coef.rho)
#		V <- s2*(s2*tr(t(AW) %*% AW) +
#			t(AW %*% x %*% coef.rho) %*%
#			(AW %*% x %*% coef.rho)) + omega*s2^2
		V <- s2*(s2*tr(crossprod(AW)) +
			crossprod(AW %*% x %*% coef.rho)) +
                        sign(con$OrdVsign)*omega*s2^2
		inf1 <- rbind(n/2, s2*tr(AW), t(zero))
		inf2 <- rbind(s2*tr(AW), V, xtawxb)
#		xtx <- s2*t(x) %*% x
		xtx <- s2*crossprod(x)
		inf3 <- rbind(zero, t(xtawxb), xtx)
		inf <- cbind(inf1, inf2, inf3)
		varb <- try(solve(inf, tol=tol.solve), silent=TRUE)
                if (inherits(varb, "try-error")) {
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                    con$fdHess <- TRUE
                    force_assign_eigen <- TRUE
                    warning(paste("inversion of asymptotic covariance",
                        "matrix failed for tol.solve =", tol.solve,
                        "\n", strsplit(attr(varb, "condition")$message,
                            ":")[[1]][2], "- using numerical Hessian."))
                } else {
                    varb <- (s2^2) * varb
		    rownames(varb) <- colnames(varb) <- 
			c("sigma", "rho", colnames(x))
		    rest.se <- sqrt(diag(varb))[-c(1:2)]
		    rho.se <- sqrt(varb[2,2])
		    TW <- (W %*% W) + crossprod(W)
		    T22 <- sum(diag(TW))
		    T21A <- sum(diag(TW %*% A))
		    LMtest <- ((t(r) %*% W %*% r)/s2)^2
		    LMtest <- LMtest/(T22 - ((T21A^2)*(rho.se^2)))
		    ase <- TRUE
                    timings[["eigen_se"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
	}

        if (con$fdHess) {
            coefs <- c(rho, coef.rho)
            if (con$compiled_sse) {
               ptr <- .Call("hess_lag_init", PACKAGE="spatialreg")
               assign("ptr", ptr, envir=env)
            }
            fdHess <- getVmatl(coefs, env,
               s2, trs, tol.solve=tol.solve, optim=con$optimHess,
               optimM=con$optimHessMethod)
            if (con$compiled_sse) {
                .Call("hess_lag_free", get("ptr", envir=env),
                     PACKAGE="spatialreg")
            }
            if (is.null(trs)) {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("rho", colnames(x))
            } else {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("sigma2", "rho", colnames(x))
            }
            if (is.null(trs)) {
 	        if (method != "eigen" || force_assign_eigen) {
                    rest.se <- sqrt(diag(fdHess)[-1])
		    rho.se <- sqrt(fdHess[1,1])
                }
            } else {
 	        if (method != "eigen" || force_assign_eigen) {
 	            rest.se <- sqrt(diag(fdHess)[-c(1,2)])
	            rho.se <- sqrt(fdHess[2,2])
                }
            }

            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
        } else fdHess <- FALSE
	call <- match.call()
        if (method=="SE_classic") {
            iC <- get("intern_classic", envir=env)
        } else iC <- NULL
#FIXME
	ret <- structure(list(type=type, dvars=dvars, rho=rho, 
		coefficients=coef.rho, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), #lm.model=lm.null,
                logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
		method=method, call=call, residuals=r, opt=optres,
                tarX=tarX, tary=tary, y=y, X=x,
		#lm.target=lm.lag, 
                fitted.values=fit,
		se.fit=NULL, #formula=formula,
                similar=similar,
		ase=ase, rho.se=rho.se, LMtest=LMtest, 
		resvar=varb, zero.policy=zero.policy, aliased=aliased,
                listw_style=listw$style, interval=interval, fdHess=fdHess,
                optimHess=con$optimHess, insert=!is.null(trs), trs=trs,
                LLNullLlm=LL_null_lm,
                timings=do.call("rbind", timings)[, c(1, 3)], 
                f_calls=get("f_calls", envir=env),
                hf_calls=get("hf_calls", envir=env), intern_classic=iC),
                class=c("Sarlm"))
        rm(env)
        GC <- gc()
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0L)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
        if (is.formula(Durbin)) attr(ret, "Durbin") <- deparse(Durbin)
	ret
}

sar.lag.mixed.f <- function(rho, env) {
        e.a <- get("e.a", envir=env)
        e.b <- get("e.b", envir=env)
        e.c <- get("e.c", envir=env)
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
        n <- get("n", envir=env)
	s2 <- SSE/n
	ldet <- do_ldet(rho, env)
	ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)
		- (1/(2*s2))*SSE)
	if (get("verbose", envir=env)) cat("rho:\t", rho, "\tfunction value:\t", ret, "\n")
        assign("f_calls", get("f_calls", envir=env)+1L, envir=env)

	ret
}



# Copyright 2010-12 by Roger Bivand
sacsarlm <- function(formula, data = list(), listw, listw2=NULL, na.action, 
	Durbin, type, method="eigen", quiet=NULL, zero.policy=NULL, 
	tol.solve=.Machine$double.eps, llprof=NULL, interval1=NULL, interval2=NULL,
        trs1=NULL, trs2=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(fdHess=NULL, #LAPACK=FALSE,
           Imult=2L, cheb_q=5L, MC_p=16L, MC_m=30L, spamPivot="MMD",
           in_coef=0.1, super=NULL, opt_method="nlminb", opt_control=list(),
           pars=NULL, npars=4L, pre_eig1=NULL, pre_eig2=NULL)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spatialregOptions)
#FIXME
        if (missing(type)) type <- "sac"
        if (type == "Durbin") type <- "sacmixed"
        if (missing(Durbin)) Durbin <- ifelse(type == "sac", FALSE, TRUE)
        if (listw$style != "W" && is.formula(Durbin)) {
            Durbin <- TRUE
            warning("formula Durbin requires row-standardised weights; set TRUE")
        }
        if (is.logical(Durbin) && isTRUE(Durbin)) type <- "sacmixed"
        if (is.formula(Durbin)) type <- "sacmixed"
        if (is.logical(Durbin) && !isTRUE(Durbin)) type <- "sac"
	switch(type, sac = if (!quiet) cat("\nSpatial ARAR model\n"),
	    sacmixed = if (!quiet) cat("\nSpatial ARAR mixed model (Manski)\n"),
	    stop("\nUnknown model type\n"))
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
        stopifnot(is.logical(zero.policy))
        if (!inherits(formula, "formula")) formula <- as.formula(formula)
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
        have_factor_preds <- have_factor_preds_mf(mf)
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(listw2)) listw2 <- listw
        if (!is.null(con$pre_eig1) && is.null(con$pre_eig2))
            con$pre_eig2 <- con$pre_eig1
        else if (!inherits(listw2, "listw")) stop("No 2nd neighbourhood list")
        if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        if (!is.null(con$pars)) {
            stopifnot(is.numeric(con$pars))
        }
#        stopifnot(is.integer(con$npars))
#        stopifnot(is.logical(con$fdHess))
#        stopifnot(is.logical(con$LAPACK))
#        stopifnot(is.logical(con$super))
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
            if (!is.null(con$pre_eig1)) {
                warning("NAs found, precomputed eigenvalues ignored")
                con$pre_eig1 <- NULL
            }
	}
	can.sim2 <- FALSE
	if (listw2$style %in% c("W", "S")) 
		can.sim2 <- can.be.simmed(listw2)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw2$neighbours) %in% na.act)
	    listw2 <- subset(listw2, subset, zero.policy=zero.policy)
            if (!is.null(con$pre_eig2)) {
                warning("NAs found, precomputed eigenvalues ignored")
                con$pre_eig2 <- NULL
            }
	}
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	n <- NROW(x)
	m <- NCOL(x)
        dvars <- c(m, 0L)
#	if (type != "sac") {
	if (is.formula(Durbin) || isTRUE(Durbin)) {
                prefix <- "lag"
                if (isTRUE(Durbin)) {
                    if (have_factor_preds) warn_factor_preds(have_factor_preds)
                    WX <- create_WX(x, listw, zero.policy=zero.policy,
                        prefix=prefix)
                } else {
                    data1 <- data
                    if (!is.null(na.act) && (inherits(na.act, "omit") ||
                        inherits(na.act, "exclude"))) {
                        data1 <- data1[-c(na.act),]
                    }
	            dmf <- lm(Durbin, data1, na.action=na.fail, 
		        method="model.frame")
	            formula_durbin_factors <- have_factor_preds_mf(dmf)
                    if (formula_durbin_factors) 
                        warn_factor_preds(formula_durbin_factors)
#	            dmf <- lm(Durbin, data, na.action=na.action, 
#		        method="model.frame")
                    fx <- try(model.matrix(Durbin, dmf), silent=TRUE)
                    if (inherits(fx, "try-error")) 
                        stop("Durbin variable mis-match")
                    WX <- create_WX(fx, listw, zero.policy=zero.policy,
                        prefix=prefix)
                    inds <- match(substring(colnames(WX), 5,
	                nchar(colnames(WX))), colnames(x))
                    if (anyNA(inds)) stop("WX variables not in X: ",
                        paste(substring(colnames(WX), 5,
                        nchar(colnames(WX)))[is.na(inds)], collapse=" "))
                    icept <- grep("(Intercept)", colnames(x))
                    iicept <- length(icept) > 0L
                    if (iicept) {
                        xn <- colnames(x)[-1]
                    } else {
                        xn <- colnames(x)
                    }
                    wxn <- substring(colnames(WX), nchar(prefix)+2,
                        nchar(colnames(WX)))
                    zero_fill <- NULL
                    if (length((which(!(xn %in% wxn)))) > 0L)
                        zero_fill <- length(xn) + (which(!(xn %in% wxn)))
                }
                dvars <- c(NCOL(x), NCOL(WX))
                if (is.formula(Durbin)) {
                    attr(dvars, "f") <- Durbin
                    attr(dvars, "inds") <- inds
                    attr(dvars, "zero_fill") <- zero_fill
                }
		x <- cbind(x, WX)
		m <- NCOL(x)
		rm(WX)
	}
	if (NROW(x) != length(listw2$neighbours))
	    stop("Input data and neighbourhood list2 have different dimensions")
	w2y <- lag.listw(listw2, y, zero.policy=zero.policy)
	w2wy <- lag.listw(listw2, wy, zero.policy=zero.policy)
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
          if (is.formula(Durbin)) {
	    stop("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
          } else {
	    warning("Aliased variables found: ",
                paste(names(aliased)[aliased], collapse=" "))
	    nacoef <- which(aliased)
		x <- x[,-nacoef]
          }
	}
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	if (m > 1 || (m == 1 && K == 1)) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw2, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw2, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)

#        env <- new.env(parent=globalenv())
        env <- new.env()
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("w2y", w2y, envir=env)
        assign("w2wy", w2wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("family", "SAR", envir=env)
        assign("first_time", TRUE, envir=env)
#        assign("LAPACK", con$LAPACK, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("can.sim2", can.sim2, envir=env)
        assign("listw", listw, envir=env)
        assign("listw2", listw2, envir=env)
        assign("similar", FALSE, envir=env)
        assign("similar2", FALSE, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (!quiet) cat(paste("\nSpatial autoregressive error model\n", 
		"Jacobian calculated using "))

        interval1 <- jacobianSetup(method, env, con, pre_eig=con$pre_eig1,
            trs=trs1, interval=interval1, which=1)
        assign("interval1", interval1, envir=env)
        interval2 <- jacobianSetup(method, env, con, pre_eig=con$pre_eig2,
            trs=trs2, interval=interval2, which=2)
        assign("interval2", interval2, envir=env)

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

        pars <- con$pars
        lower <- c(interval1[1], interval2[1])
        upper <- c(interval1[2], interval2[2])

        if (!is.null(llprof)) {
            llrho <- NULL
            lllambda <- NULL
            if (length(llprof) == 1L) {
                llrho <- seq(lower[1], upper[1], length.out=llprof)
                lllambda <- seq(lower[2], upper[2], length.out=llprof)
                llprof <- as.matrix(expand.grid(llrho, lllambda))
            }
            ll_prof <- numeric(nrow(llprof))
            for (i in 1:nrow(llprof)) 
                ll_prof[i] <- sacsar.f(llprof[i,], env=env)
            nm <- paste(method, "profile", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
        }

        if (is.null(pars)) {
          if (con$npars == 4L) {
             xseq <- c(lower[1], 0, upper[1], upper[1])*0.8
             yseq <- c(upper[2], 0, upper[2], lower[2])*0.8
             mpars <- cbind(xseq, yseq)
          } else {
             xseq <- seq(lower[1], upper[1], (upper[1]-lower[1])/2)*0.8
             yseq <- seq(lower[2], upper[2], (upper[2]-lower[2])/2)*0.8
             mpars <- as.matrix(expand.grid(xseq, yseq))
          }
        } else {
            mxs <- NULL
        }
        if (con$opt_method == "nlminb") {
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -nlminb(pp, sacsar.f, 
                    env=env, control=con$opt_control, lower=lower, 
                    upper=upper)$objective)
                pars <- mpars[which.max(mxs),]
                optres <- nlminb(pars, sacsar.f, env=env,
                    control=con$opt_control, lower=lower, upper=upper)
            } else {
                optres <- nlminb(pars, sacsar.f, env=env,
                    control=con$opt_control, lower=lower, upper=upper)
            }
        } else if (con$opt_method == "L-BFGS-B"){
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -optim(pars, sacsar.f, 
                    env=env, method="L-BFGS-B", control=con$opt_control, 
                    lower=lower, upper=upper)$objective)
                pars <- mpars[which.max(mxs),]
	        optres <- optim(pars, sacsar.f, env=env,
                    method="L-BFGS-B", control=con$opt_control,
                    lower=lower, upper=upper)
            } else {
	        optres <- optim(pars, sacsar.f, env=env,
                    method="L-BFGS-B", control=con$opt_control,
                    lower=lower, upper=upper)
            }
        } else {
            if (is.null(pars)) {
                mxs <- apply(mpars, 1, function(pp) -optim(pars, sacsar.f, 
                    env=env, method=con$opt_method, 
                    control=con$opt_control)$objective)
                pars <- mpars[which.max(mxs),]
	        optres <- optim(pars, sacsar.f, env=env,
                    method=con$opt_method, control=con$opt_control)
            } else {
	        optres <- optim(pars, sacsar.f, env=env,
                    method=con$opt_method, control=con$opt_control)
            }
        }
	LL <- -optres$objective
        if (optres$convergence != 0)
            warning(paste("convergence failure:", optres$message))
	rho <- optres$par[1]
	names(rho) <- "rho"
	lambda <- optres$par[2]
	names(lambda) <- "lambda"
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	lm.target <- lm(I(y - rho*wy - lambda*w2y + rho*lambda*w2wy) ~ 
            I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	coef.sac <- coefficients(lm.target)
        tarX <- model.matrix(lm.target)
        tary <- model.response(model.frame(lm.target))
	names(coef.sac) <- xcolnames
	lm.model <- lm(formula, data)
        logLik_lm.model <- logLik(lm.model)
        AIC_lm.model <- AIC(lm.model)
        ase <- FALSE
	asyvar1 <- FALSE
        force_assign_eigen <- FALSE
        if (method == "eigen") {
# taken from spatial/sac_models/sac.m
	    tr <- function(A) sum(diag(A))
            W1 <- listw2mat(listw)
            W2 <- listw2mat(listw2)
            A <- diag(n) - rho * W1
            AI <- solve(A)
            WA <- W1 %*% AI
            B <- diag(n) - lambda * W2
            BI <- solve(B)
            WB <- W2 %*% BI
            omeg <- s2*diag(n)
            omegi <- (1/s2)*diag(n)
            bhat <- coef.sac
            p3 <- p+3
            asyvar <- matrix(0.0, ncol=p3, nrow=p3)
            Bx <- B %*% x
            asyvar[4:p3, 4:p3] <- (1/s2)*crossprod(Bx)
            asyvar[1, 1] <- n/(2*s2*s2)
            term1 <- tr(WA %*% WA)
            BWABI <- B %*% WA %*% BI
            term2 = tr(omeg %*% t(BWABI) %*% omegi %*% BWABI)
            BWAxbhat <- B %*% WA %*% x %*% bhat
            term3 = t(BWAxbhat) %*% omegi %*% BWAxbhat
            asyvar[2, 2] <- term1 + term2 + term3
            term1 <- tr(WB %*% WB)
            term2 <- tr(omeg %*% t(WB) %*% omegi %*% WB)
            asyvar[3, 3] <- term1 + term2
            asyvar[2, 4:p3] <- (1/s2)*(t(Bx) %*% BWAxbhat)
            asyvar[4:p3, 2] <- asyvar[2, 4:p3]
            asyvar[2, 1] <- (1/s2)*tr(W1 %*% AI)
            asyvar[1, 2] <- asyvar[2, 1]
            asyvar[3, 1] <- (1/s2)*tr(W2 %*% BI)
            asyvar[1, 3] <- asyvar[3, 1]
            term1 = tr(t(WB) %*% omegi %*% BWABI %*% omeg)
            term2 = tr(W2 %*% WA %*% BI)
            asyvar[3, 2] <- term1 + term2
            asyvar[2, 3] <- asyvar[3, 2]
            asyvar1 <- try(solve(asyvar, tol.solve=tol.solve), silent=TRUE)
            if (inherits(asyvar1, "try-error")) {
                timings[["eigen_se"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                con$fdHess <- TRUE
                force_assign_eigen <- TRUE
                warning(paste("inversion of asymptotic covariance",
                    "matrix failed for tol.solve =", tol.solve,
                    "\n", strsplit(attr(asyvar1, "condition")$message,
                        ":")[[1]][2], "- using numerical Hessian."))
                ase=FALSE
            } else {
                rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "rho", "lambda", xcolnames)
                ase=TRUE
                rho.se <- sqrt(asyvar1[2, 2])
                lambda.se <- sqrt(asyvar1[3, 3])
                rest.se <- sqrt(diag(asyvar1)[-c(1:3)])
                nm <- "asymptotic vcov"
                timings[[nm]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
            }
        }
        if (con$fdHess) {
            coefs <- c(rho, lambda, coef.sac)
            fdHess <- getVmatsac(coefs, env, tol.solve=tol.solve)
            rownames(fdHess) <- colnames(fdHess) <- c("rho", "lambda",
                xcolnames)
            if (!ase) {
                rho.se <- sqrt(fdHess[1, 1])
                lambda.se <- sqrt(fdHess[2, 2])
                rest.se <- sqrt(diag(fdHess)[-c(1,2)])
            }
            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
        } else fdHess <- FALSE
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type=type, dvars=dvars, rho=rho, lambda=lambda,
	    coefficients=coef.sac, rest.se=rest.se, ase=ase,
	    LL=LL, s2=s2, SSE=SSE, parameters=(p+3), 
            logLik_lm.model=logLik_lm.model, AIC_lm.model=AIC_lm.model,
            #lm.model=lm.model, 
	    method=method, call=call, residuals=r, #lm.target=lm.target,
            tarX=tarX, tary=tary, y=y, X=x, W2X=WX, trs1=trs1, trs2=trs2,
	    opt=optres, pars=pars, mxs=mxs, fitted.values=fit, #formula=formula,
	    similar=get("similar", envir=env), rho.se=rho.se,
	    lambda.se=lambda.se, zero.policy=zero.policy, 
	    aliased=aliased, LLNullLlm=LL_null_lm,
            fdHess=fdHess, resvar=asyvar1, listw_style=listw$style,
            optimHess=FALSE, insert=FALSE, interval1=interval1,
            interval2=interval2, timings=do.call("rbind", timings)[, c(1, 3)]),
            class=c("Sarlm"))
        rm(env)
        GC <- gc()
        if (is.null(llprof)) ret$llprof <- llprof
        else {
            ret$llprof <- list(grd=llprof, ll=ll_prof, xseq=llrho,
                yseq=lllambda)
        }
        if (is.formula(Durbin)) attr(ret, "Durbin") <- deparse(Durbin)
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret

}


sacsar_sse <- function(coefs, env) {
    rho <- coefs[1]
    lambda <- coefs[2]
    yl <- get("y", envir=env) - rho * get("wy", envir=env) - 
        lambda * get("w2y", envir=env) + rho * lambda * get("w2wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    xl.q <- qr.Q(qr(xl, #LAPACK=get("LAPACK", envir=env)
))
    xl.q.yl <- crossprod(xl.q, yl)
    SSE <- crossprod(yl) - crossprod(xl.q.yl)
    SSE
}


sacsar.f <- function(coefs, env) {
    SSE <- sacsar_sse(coefs, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet1 <- do_ldet(coefs[1], env, which=1)
    ldet2 <- do_ldet(coefs[2], env, which=2)
    ret <- (ldet1 + ldet2 - ((n/2)*log(2*pi)) - (n/2)*log(s2) - 
        (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("rho:", coefs[1], " lambda:", coefs[2],
        " function:", ret, " Jacobian1:", ldet1, " Jacobian2:", ldet2,
        " SSE:", SSE, "\n")
    -ret
}


