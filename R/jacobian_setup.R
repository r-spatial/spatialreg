# Copyright 2012 by Roger Bivand

jacobianSetup <- function(method, env, con, pre_eig=NULL, trs=NULL, interval=NULL, which=1) {
    switch(method,
        eigen = {
            if (get("verbose", envir=env))
                cat("neighbourhood matrix eigenvalues\n")
            if (is.null(pre_eig)) {
                eigen_setup(env, which=which)
            } else {
                eigen_pre_setup(env, pre_eig=pre_eig, which=which)
            }
            er <- get("eig.range", envir=env)
            if (is.null(interval)) 
                interval <- c(er[1]+.Machine$double.eps,
                              er[2]-.Machine$double.eps)
        },
        Matrix = {
            if (get("listw", envir=env)$style %in% c("W", "S") &&
                !get("can.sim", envir=env))
                stop("Matrix method requires symmetric weights")
            if (get("listw", envir=env)$style %in% c("B", "C", "U") && 
                !(is.symmetric.glist(get("listw", envir=env)$neighbours,
                get("listw", envir=env)$weights)))
                stop("Matrix method requires symmetric weights")
            if (get("verbose", envir=env))
                cat("sparse matrix Cholesky decomposition\n")
            Imult <- con$Imult
            if (is.null(interval)) {
                if (get("listw", envir=env)$style == "B") {
                    interval <- c(-0.5, +0.25)
                } else interval <- c(-1, 0.999)
            }
            if (get("listw", envir=env)$style == "B") {
                Imult <- ceiling((2/3) * max(sapply(get("listw",
                    envir=env)$weights, sum)))
            }

            if (is.null(con$super)) con$super <- as.logical(NA)
            Matrix_setup(env, Imult, con$super, which=which)
        },
        Matrix_J = {
            if (get("listw", envir=env)$style %in% c("W", "S") &&
                !get("can.sim", envir=env))
                stop("Matrix method requires symmetric weights")
            if (get("listw", envir=env)$style %in% c("B", "C", "U") && 
                !(is.symmetric.glist(get("listw", envir=env)$neighbours,
                get("listw", envir=env)$weights)))
                stop("Matrix method requires symmetric weights")
            if (get("verbose", envir=env))
                cat("sparse matrix Cholesky decomposition\n")
            if (is.null(interval)) {
                if (get("listw", envir=env)$style == "B") {
                    interval <- c(-0.5, +0.25)
                } else interval <- c(-1, 0.999)
            }
            if (is.null(con$super)) con$super <- FALSE
            Matrix_J_setup(env, super=con$super, which=which)
        },
        spam = {
#            if (!require(spam)) stop("spam not available")
          if (requireNamespace("spam", quietly = TRUE)) {
            if (get("listw", envir=env)$style %in% c("W", "S") &&
                !get("can.sim", envir=env))
                stop("spam method requires symmetric weights")
            if (get("listw", envir=env)$style %in% c("B", "C", "U") && 
                !(is.symmetric.glist(get("listw", envir=env)$neighbours,
                get("listw", envir=env)$weights)))
                stop("spam method requires symmetric weights")
            if (get("verbose", envir=env))
                cat("sparse matrix Cholesky decomposition\n")
            spam_setup(env, pivot=con$spamPivot, which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
          } else {
            stop("spam not available")
          }
        },
        spam_update = {
#            if (!require(spam)) stop("spam not available")
          if (requireNamespace("spam", quietly = TRUE)) {
            if (get("listw", envir=env)$style %in% c("W", "S") &&
                !get("can.sim", envir=env))
                stop("spam method requires symmetric weights")
            if (get("listw", envir=env)$style %in% c("B", "C", "U") && 
                !(is.symmetric.glist(get("listw", envir=env)$neighbours,
                get("listw", envir=env)$weights)))
                stop("spam method requires symmetric weights")
            if (get("verbose", envir=env)) 
                cat("sparse matrix Cholesky decomposition\n")
            spam_update_setup(env, in_coef=con$in_coef,
                 pivot=con$spamPivot, which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
          } else {
            stop("spam not available")
          }
        },
        Chebyshev = {
            if (get("listw", envir=env)$style %in% c("W", "S") &&
                !get("can.sim", envir=env))
                stop("Chebyshev method requires symmetric weights")
            if (get("listw", envir=env)$style %in% c("B", "C", "U") && 
                !(is.symmetric.glist(get("listw", envir=env)$neighbours,
                get("listw", envir=env)$weights)))
                stop("Chebyshev method requires symmetric weights")
            if (get("verbose", envir=env)) 
                cat("sparse matrix Chebyshev approximation\n")
            cheb_setup(env, q=con$cheb_q, which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
        },
        MC = {
            if (!get("listw", envir=env)$style %in% c("W"))
                stop("MC method requires row-standardised weights")
            if (get("verbose", envir=env)) 
                cat("sparse matrix Monte Carlo approximation\n")
            mcdet_setup(env, p=con$MC_p, m=con$MC_m, which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
        },
        LU = {
            if (get("verbose", envir=env))
                cat("sparse matrix LU decomposition\n")
            LU_setup(env, which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
        },
        LU_prepermutate = {
            if (get("verbose", envir=env))
                cat("sparse matrix LU decomposition\n")
            LU_prepermutate_setup(env, coef=con$in_coef, order=con$LU_order,
                which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
        },
        moments = {
            if (get("verbose", envir=env))
                cat("Smirnov/Anselin (2009) trace approximation\n")
            moments_setup(env, trs=trs, m=con$MC_m, p=con$MC_p,
                type=con$type, correct=con$correct, trunc=con$trunc,
                which=which)
            if (is.null(interval)) interval <- c(-1,0.999)
       },
       SE_classic = {
            if (get("verbose", envir=env)) 
                cat("SE toolbox classic grid\n")
            if (is.null(interval)) interval <- c(-1,0.999)
            if (con$SE_method == "MC" && 
                !get("listw", envir=env)$style %in% c("W"))
                stop("MC method requires row-standardised weights")
            SE_classic_setup(env, SE_method=con$SE_method, p=con$MC_p,
                m=con$MC_m, nrho=con$nrho, interpn=con$interpn,
                interval=interval, SElndet=con$SElndet, which=which)
       },
       SE_whichMin = {
            if (get("verbose", envir=env))
                cat("SE toolbox which.min grid\n")
            if (is.null(interval)) interval <- c(-1,0.999)
            if (con$SE_method == "MC" &&
                !get("listw", envir=env)$style %in% c("W"))
                stop("MC method requires row-standardised weights")
            SE_whichMin_setup(env, SE_method=con$SE_method, p=con$MC_p,
                m=con$MC_m, nrho=con$nrho, interpn=con$interpn,
                interval=interval, SElndet=con$SElndet, which=which)
        },
        SE_interp = {
            if (get("verbose", envir=env))
                cat("SE toolbox which.min grid\n")
            if (is.null(interval)) interval <- c(-1,0.999)
            if (con$SE_method == "MC" &&
                !get("listw", envir=env)$style %in% c("W"))
                stop("MC method requires row-standardised weights")
            SE_interp_setup(env, SE_method=con$SE_method, p=con$MC_p,
                m=con$MC_m, nrho=con$nrho, interval=interval,
                which=which)
        },
        stop("...\n\nUnknown method\n"))
    interval
}

can.be.simmed <- function(listw) {
	res <- is.symmetric.nb(listw$neighbours, FALSE)
        if (is.null(attr(listw$weights, "comp")))
            stop("non-compliant listw object, re-build using spdep::nb2listw()")
	if (res) {
		if (attr(listw$weights, "mode") == "general")
			res <- attr(listw$weights, "glistsym")
	} else return(res)
	res
}


similar.listw_Matrix <- function(listw) {
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	ww <- as(listw, "CsparseMatrix")
	if (listw$style == "W") {
		d <- attr(listw$weights, "comp")$d
		d1 <- 1/(sqrt(d))
		dd <- as(as(Diagonal(x=d), "symmetricMatrix"), "CsparseMatrix")
		dd1 <- as(as(Diagonal(x=d1), "symmetricMatrix"),
		    "CsparseMatrix")
		ww1 <- dd %*% ww
		res <- dd1 %*% ww1 %*% dd1
	} else if (listw$style == "S") {
		q <- attr(listw$weights, "comp")$q
		Q <- attr(listw$weights, "comp")$Q
		eff.n <- attr(listw$weights, "comp")$eff.n
		q1 <- 1/(sqrt(q))
		qq <- as(as(Diagonal(x=q), "symmetricMatrix"), "CsparseMatrix")
		qq1 <- as(as(Diagonal(x=q1), "symmetricMatrix"),
		    "CsparseMatrix")
		ww0 <- (Q/eff.n) * ww
		ww1 <- qq %*% ww0
		sim0 <- qq1 %*% ww1 %*% qq1
		res <- (eff.n/Q) * sim0
	} else stop("Conversion not suitable for this weights style")
	res
}


similar.listw_spam <- function(listw) {
    if (requireNamespace("spam", quietly = TRUE)) {
#        if (!require(spam)) stop("spam not available")
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	sww <- as.spam.listw(listw)
	if (listw$style == "W") {
		sd <- attr(listw$weights, "comp")$d
		sd1 <- 1/(sqrt(sd))
                if (any(!is.finite(sd1))) {
                    sd1[!is.finite(sd1)] <- 0
                    warning("non-finite inverse diagonal values set to zero")
                }
		sdd <- spam::diag.spam(sd, n, n)
		sdd1 <- spam::diag.spam(sd1, n, n)
		sww1 <- sdd %*% sww
		res <- sdd1 %*% sww1 %*% sdd1
	} else if (listw$style == "S") {
		q <- attr(listw$weights, "comp")$q
		Q <- attr(listw$weights, "comp")$Q
		eff.n <- attr(listw$weights, "comp")$eff.n
		q1 <- 1/(sqrt(q))
                if (any(!is.finite(q1))) {
                    sd1[!is.finite(q1)] <- 0
                    warning("non-finite inverse diagonal values set to zero")
                }
		qq <- spam::diag.spam(q, n, n)
		qq1 <- spam::diag.spam(q1, n, n)
		ww0 <- (Q/eff.n) * sww
		ww1 <- qq %*% ww0
		sim0 <- qq1 %*% ww1 %*% qq1
		res <- (eff.n/Q) * sim0
	} else stop("Conversion not suitable for this weights style")
	return(res)
    } else {
        stop("spam not available")
    }
}

similar.listw <- function(listw) {
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(listw$neighbours)
	if (listw$style == "W") {
		d <- attr(listw$weights, "comp")$d
		glist <- vector(mode="list", length=n)
		for (i in 1:n) glist[[i]] <- d[i] * listw$weights[[i]]
		sd1 <- 1/sqrt(d)
		for (i in 1:n) {
			inb <- listw$neighbours[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
				for (j in 1:icd) {
					glist[[i]][j] <- sd1[i] * 
						glist[[i]][j] * sd1[inb[j]]
				}
			}
		}
		res <- listw
		res$weights <- glist
		attr(res$weights, "mode") <- "sim"
		attr(res$weights, "W") <- TRUE
		attr(res$weights, "comp") <- attr(listw$weights, "comp")
		res$style <- "W:sim"
	} else if (listw$style == "S") {
		q <- attr(listw$weights, "comp")$q
		Q <- attr(listw$weights, "comp")$Q
		eff.n <- attr(listw$weights, "comp")$eff.n
		glist <- vector(mode="list", length=n)
		for (i in 1:n) {
			glist[[i]] <- (Q/eff.n) * listw$weights[[i]]
			glist[[i]] <- q[i] * glist[[i]]
		}
		sq1 <- 1/sqrt(q)
		for (i in 1:n) {
			inb <- listw$neighbours[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
				for (j in 1:icd) {
					glist[[i]][j] <- sq1[i] * 
						glist[[i]][j] * sq1[inb[j]]
				}
				glist[[i]] <- (eff.n/Q) * glist[[i]]
			}
		}
		res <- listw
		res$weights <- glist
		attr(res$weights, "mode") <- "sim"
		attr(res$weights, "S") <- TRUE
		attr(res$weights, "comp") <- attr(listw$weights, "comp")
		res$style <- "S:sim"
	} else stop("Conversion not suitable for this weights style")
	sym_out <- is.symmetric.glist(res$neighbours, res$weights)
	if (!sym_out) {
	    if (attr(sym_out, "d") < .Machine$double.eps ^ 0.5)
		res <- listw2U(res)
	    else stop("defective similarity")
	}
	res
}

