

lmSLX <- function(formula, data = list(), listw, na.action, weights=NULL, Durbin=TRUE, zero.policy=NULL) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
        stopifnot(is.logical(zero.policy))
        if (!inherits(formula, "formula")) formula <- as.formula(formula)
#	mt <- terms(formula, data = data)
#	mf <- lm(formula, data, na.action=na.action, weights=weights,
#            method="model.frame")
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")

	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
        
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
        n <- nrow(x)
        if (n != length(listw$neighbours))
            stop("listw and data of different lengths")
        nclt <- colnames(x)

        weights <- as.vector(model.extract(mf, "weights"))
        if (!is.null(weights) && !is.numeric(weights)) 
            stop("'weights' must be a numeric vector")
        if (is.null(weights)) weights <- rep(as.numeric(1), n)
        if (any(is.na(weights))) stop("NAs in weights")
        if (any(weights < 0)) stop("negative weights")
        dvars <- c(NCOL(x), 0L)
        prefix <- "lag"
        if (isTRUE(Durbin)) {
            WX <- create_WX(x, listw, zero.policy=zero.policy,
               prefix=prefix)
        } else if (is.formula(Durbin)) {
            data1 <- data
            if (!is.null(na.act) && (inherits(na.act, "omit") ||
                inherits(na.act, "exclude"))) {
                data1 <- data1[-c(na.act),]
            }
            dmf <- lm(Durbin, data1, na.action=na.fail, 
	        method="model.frame")
#	    dmf <- lm(Durbin, data, na.action=na.action, 
#	         method="model.frame")
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
            zero_fill <- length(xn) + (which(!(xn %in% wxn)))
        } else stop("Durbin argument neither TRUE nor formula")
        dvars <- c(NCOL(x), NCOL(WX))
        if (is.formula(Durbin)) {
            attr(dvars, "f") <- Durbin
            attr(dvars, "inds") <- inds
            attr(dvars, "zero_fill") <- zero_fill
        }
	x <- cbind(x, WX)
	rm(WX)
#        WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
#        x <- cbind(x, WX)
# 180128 Mark L. Burkey summary.lm error for SLX object
        colnames(x) <- make.names(colnames(x))
        if (attr(mt, "intercept") == 1L) {
            lm.model <- lm(formula(paste("y ~ ", paste(colnames(x)[-1], collapse="+"))), data=as.data.frame(x), weights=weights)
        } else {
            lm.model <- lm(formula(paste("y ~ 0 + ", paste(colnames(x), collapse="+"))), data=as.data.frame(x), weights=weights)
        }
        sum_lm_model <- summary.lm(lm.model, correlation = FALSE)
        mixedImps <- NULL
	K <- ifelse(isTRUE(grep("Intercept",
            names(coefficients(lm.model))[1]) == 1L), 2, 1)
        if (isTRUE(Durbin)) {
          m <- length(coefficients(lm.model))
          odd <- (m%/%2) > 0
          if (odd) {
              m2 <- (m-1)/2
          } else {
              m2 <- m/2
          }
          if (K == 1 && odd) {
            warning("model configuration issue: no total impacts")
          } else {
            cm <- matrix(0, ncol=m, nrow=m2)
            if (K == 2) {
                if (odd) {
                    rownames(cm) <- nclt[2:(m2+1)]
                } else {
                    rownames(cm) <- nclt[1:m2]
                }
                for (i in 1:m2) cm[i, c(i+1, i+(m2+1))] <- 1
# drop bug fix 2016-09-21 Philipp Hunziker
                dirImps <- sum_lm_model$coefficients[2:(m2+1), 1:2, drop=FALSE]
                rownames(dirImps) <- rownames(cm)
                indirImps <- sum_lm_model$coefficients[(m2+2):m, 1:2, drop=FALSE]
                rownames(indirImps) <- rownames(cm)
            } else {
                rownames(cm) <- nclt[1:m2] # FIXME
                for (i in 1:m2) cm[i, c(i, i+m2)] <- 1
                dirImps <- sum_lm_model$coefficients[1:m2, 1:2, drop=FALSE]
                rownames(dirImps) <- rownames(cm)
                indirImps <- sum_lm_model$coefficients[(m2+1):m, 1:2, drop=FALSE]
                rownames(indirImps) <- rownames(cm)
            }
            totImps <- as.matrix(estimable(lm.model, cm)[, 1:2, drop=FALSE])
            
            
          } 
      } else if (is.formula(Durbin)) {
#FIXME
            m <- sum(dvars)
            m2 <- dvars[2]
            cm <- matrix(0, ncol=m, nrow=m2)
            for (i in 1:m2) {
                cm[i, c(inds[i], i+dvars[1])] <- 1
            }
            rownames(cm) <- wxn
            dirImps <- sum_lm_model$coefficients[2:dvars[1], 1:2,
              drop=FALSE]
            rownames(dirImps) <- xn
            indirImps <- sum_lm_model$coefficients[(dvars[1]+1):m, 1:2,
              drop=FALSE]
            if (!is.null(zero_fill)) {
              if (length(zero_fill) > 0L) {
               lres <- vector(mode="list", length=2L)
               for (j in 1:2) {
                 jindirImps <- rep(as.numeric(NA), (dvars[1]-1))
                   for (i in seq(along=inds)) {
                     jindirImps[(inds[i]-1)] <- indirImps[i, j]
                   }
                     lres[[j]] <- jindirImps
                   }
                   indirImps <- do.call("cbind", lres)
                 }
               }
               rownames(indirImps) <- xn
               totImps <- as.matrix(estimable(lm.model, cm)[, 1:2,
                 drop=FALSE])
               if (!is.null(zero_fill)) {
                 if (length(zero_fill) > 0L) {
                   lres <- vector(mode="list", length=2L)
                   for (j in 1:2) {
                     jtotImps <- dirImps[, j]
                     for (i in seq(along=inds)) {
                       jtotImps[(inds[i]-1)] <- totImps[i, j]
                     }
                     lres[[j]] <- jtotImps
                   }
                   totImps <- do.call("cbind", lres)
                  }
                }
                rownames(totImps) <- xn
        } else stop("undefined Durbin state")
        mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
            totImps=totImps)
        
        attr(lm.model, "mixedImps") <- mixedImps
        attr(lm.model, "dvars") <- dvars
        class(lm.model) <- c("SLX", class(lm.model))
        lm.model
}


predict.SLX <- function(object, newdata, listw, zero.policy=NULL, ...) {
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
    stopifnot(is.logical(zero.policy))
    if (missing(newdata)) {
        return(fitted(object))
    }
    if (!inherits(listw, "listw")) stop("No neighbourhood list")
    if (is(newdata, "Spatial")) newdata <- as(newdata, "data.frame")
    if (!inherits(newdata, "data.frame"))
        stop("newdata must be a Spatial*DataFrame or a data.frame")
    vars <- rownames(attr(object, "mixedImps")$dirImps)
    f <- formula(paste("~", paste(vars, collapse=" + ")))
    mf <- lm(f, newdata, method="model.frame")
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    if (any(is.na(x))) stop("NAs in independent variable")
    n <- nrow(x)
    if (n != length(listw$neighbours))
        stop("listw and data of different lengths")
    WX <- create_WX(x, listw, zero.policy=zero.policy, prefix="lag")
    x <- cbind(x, WX)
    res <- as.vector(x %*% coef(object))
    names(res) <- row.names(newdata)
    res
}


impacts.SLX <- function(obj, ...) {
    stopifnot(!is.null(attr(obj, "mixedImps")))
    n <- nrow(obj$model)
    k <- obj$qr$rank
    impactsWX(attr(obj, "mixedImps"), n, k, type="SLX", method="estimable")
}


impactsWX <- function(obj, n, k, type="SLX", method="estimable") {
    imps <- lapply(obj, function(x) x[, 1])
    names(imps) <- c("direct", "indirect", "total")
    attr(imps, "bnames") <- rownames(obj[[1]])
    ses <- lapply(obj, function(x) x[, 2])
    names(ses) <- c("direct", "indirect", "total")
    attr(ses, "bnames") <- rownames(obj[[1]])
    res <- list(impacts=imps, se=ses)
    attr(res, "n") <- n
    attr(res, "k") <- k
    attr(res, "type") <- type
    attr(res, "method") <- method
    attr(res, "bnames") <- rownames(obj[[1]])
    class(res) <- "WXImpact"
    res
}


print.WXImpact <- function(x, ...) {
    mat <- lagImpactMat(x$impacts)
    cat("Impact measures (", attr(x, "type"), ", ",
        attr(x, "method"), "):\n", sep="")
    print(mat, ...)
    cat("\n")
    invisible(x)
}


print.summary.WXImpact <- function(x, ...) {
    mat <- x$mat
    cat("Impact measures (", attr(x, "type"), ", ",
        attr(x, "method"), "):\n", sep="")
    print(mat, ...)
    cat("========================================================\n")
    mat <- x$semat
    cat("Standard errors:\n", sep="")
    print(mat, ...)
    cat("========================================================\n")
    cat("Z-values:\n")
    mat <- x$zmat
    rownames(mat) <- attr(x, "bnames")
    print(mat, ...)
    cat("\np-values:\n")
    xx <- apply(x$pzmat, 2, format.pval)
# 100928 Eelke Folmer
    if (length(attr(x, "bnames")) == 1L) {
        xx <- matrix(xx, ncol=3)
        colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- attr(x, "bnames")
    print(xx, ..., quote=FALSE)
    cat("\n")
    invisible(x)
}

summary.WXImpact <- function(object, ...,
 adjust_k=(attr(object, "type") == "SDEM")) {
    stopifnot(is.logical(adjust_k))
    stopifnot(length(adjust_k) == 1L)
    object$mat <- lagImpactMat(object$impacts)
    object$semat <- lagImpactMat(object$se)
    if (adjust_k) {
        object$semat <- sqrt((object$semat^2) * ((attr(object, "n") - 
            attr(object, "k"))/attr(object, "n")))
        attr(object, "method") <- paste(attr(object, "method"),
            ", n", sep="")
    } else {
        attr(object, "method") <- paste(attr(object, "method"),
            ", n-k", sep="")
    }
    object$zmat <- object$mat/object$semat
    object$pzmat <- 2*(1-pnorm(abs(object$zmat)))
    class(object) <- c("summary.WXImpact", class(object))
    object
}


create_WX <- function(x, listw, zero.policy=NULL, prefix="") {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
        stopifnot(is.logical(zero.policy))
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
	# check if there are enough regressors
	xcolnames <- colnames(x)
        stopifnot(!is.null(xcolnames))
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
        Wvars <- NULL
        wxI <- NULL
        WX <- NULL
	if (K == 2) {
        # unnormalized weight matrices
               	if (!(listw$style == "W")) {
 			intercept <- as.double(rep(1, n))
       	       		wxI <- lag.listw(listw, intercept, 
				zero.policy = zero.policy)
                        Wvars <- paste(prefix, ".(Intercept)", sep="")
               	} 
        }   
	if (m > 1 || (m == 1 && K == 1)) {
                WX <- matrix(as.numeric(NA), nrow=n,
                    ncol=ifelse(m==1, 1, (m-(K-1))))
		for (k in K:m) {
                        j <- ifelse(k==1, 1, k-(K-1))
			WX[,j] <- lag.listw(listw, x[,xcolnames[k]], 
			    zero.policy=zero.policy)
			if (any(is.na(WX[,j]))) 
			    stop("NAs in lagged independent variable")
                        Wvars <- c(Wvars, paste(prefix, ".",
                            xcolnames[k], sep=""))
		}
	}
        if (!is.null(wxI)) WX <- cbind(wxI, WX)
        colnames(WX) <- Wvars
        WX
}

