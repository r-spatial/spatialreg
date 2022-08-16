# Copyright 2006-14 by Roger Bivand
#

as.spam.listw <- function(listw) {
    if (requireNamespace("spam", quietly = TRUE)) {
#if (!require(spam)) stop("spam not available")
        N <- length(listw$neighbours)
        W_sn <- listw2sn(listw)
        rpts <- as.integer(cumsum(c(1, card(listw$neighbours))))
        W <- new("spam", entries=W_sn$weights, colindices=W_sn$to,
            rowpointers=rpts, dimension=as.integer(c(N, N)))
        stopifnot(spam::validate_spam(W))
        return(W)
    } else stop("spam not available")
}

listw2U_spam <- function(lw) 0.5 * (lw + t(lw))


setAs("listw", "CsparseMatrix", function(from) {as(as_dgRMatrix_listw(from), "CsparseMatrix")})
setAs("listw", "RsparseMatrix", function(from) {as_dgRMatrix_listw(from)})
setAs("listw", "symmetricMatrix", function(from) {as_dsTMatrix_listw(from)})


as_dgRMatrix_listw <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	p0 <- as.integer(c(0, cumsum(cardw)))
	scard <- sum(cardw)
	z <- .Call("listw2dgR", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard), PACKAGE="spatialreg")
	res <- as(new("dgRMatrix", j=z[[1]], p=p0, Dim=as.integer(c(n, n)),
		x=z[[2]]), "generalMatrix")
        colnames(res) <- attr(listw$neighbours, "region.id")
        rownames(res) <- colnames(res)
	res
}

as_dsTMatrix_listw <- function(listw) {
	if (!inherits(listw, "listw")) stop("not a listw object")
	if (!is.symmetric.glist(listw$neighbours, listw$weights))
		stop("not a symmetric matrix")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	scard <- sum(cardw)
	if (scard %% 2 != 0) stop("odd neighbours sum")
	z <- .Call("listw2dsT", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard/2), PACKAGE="spatialreg")

	res <- as(new("dsTMatrix", i=z[[1]], j=z[[2]], Dim=as.integer(c(n, n)),
		x=z[[3]]), "generalMatrix")
        colnames(res) <- attr(listw$neighbours, "region.id")
        rownames(res) <- colnames(res)
	res
}

as_dsCMatrix_I <- function(n) {
	if (n < 1) stop("matrix must have positive dimensions")
	as(as(Diagonal(n), "symmetricMatrix"), "generalMatrix")
}

as_dsCMatrix_IrW <- function(W, rho) {
	stopifnot(is(W, "symmetricMatrix"))
	n <- dim(W)[1]
	as_dsCMatrix_I(n) - rho * W
}

Jacobian_W <- function(W, rho) {
	sum(2*log(diag(chol(as_dsCMatrix_IrW(W, rho)))))
}


listw2U_Matrix <- function(lw) 	
	as(as(0.5 * (lw + t(lw)), "symmetricMatrix"), "CsparseMatrix")


powerWeights <- function(W, rho, order=250, X, tol=.Machine$double.eps^(3/5)) {
    timings <- list()
    .ptime_start <- proc.time()
    n <- dim(W)[1]
    dX <- dim(X)
    if (dX[1] == n) side <- "R"
    else if (dX[2] == n) side <- "L"
    else stop("W and X non-conformant")
    aW <- rho*W
    if (side == "R") last <- aW %*% X
    else last <- X %*% aW
    acc <- X + last
    conv <- FALSE
    iter <- 1
    series <- numeric(order)
    while (iter < order) {
        if (side == "R") {
            last <- aW %*% last
            acc <- acc + last
        } else {
            last <- last %*% aW
            acc <- acc + last
        }
# abs() added 2017-02-15, bug spotted by Yongwan Chun
        series[iter] <- mean(abs(as(last, "matrix")))
        if (series[iter] < tol) {
            conv <- TRUE
            break
        }
        iter <- iter+1
    }
    if (!conv) warning("not converged within order iterations")
    timings[["make_power_sum"]] <- proc.time() - .ptime_start
    attr(acc, "internal") <- list(series=series, order=order,
        tol=tol, iter=iter, conv=conv)
    attr(acc, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    acc
}

invIrM <- function(neighbours, rho, glist=NULL, style="W", method="solve", 
	feasible=NULL) {
	if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
	invIrW(nb2listw(neighbours, glist=glist, style=style), rho=rho, 
		method=method, feasible=feasible)
}

invIrW <- function(x, rho, method="solve", feasible=NULL) {
	if(inherits(x, "listw")) {
	  n <- length(x$neighbours)
	  V <- listw2mat(x)
	} else if (inherits(x, "Matrix") || inherits(x, "matrix")) {
	  if (method == "chol" && all(t(x) == x))
            stop("No Cholesky method for matrix or sparse matrix object")
          n <- dim(x)[1]
          V <- x
	} else stop("Not a weights list or a Sparse Matrix")
	if (is.null(feasible) || (is.logical(feasible) && !feasible)) {
		e <- eigen(V, only.values = TRUE)$values
		if (is.complex(e)) feasible <- 1/(range(Re(e)))
		else feasible <- 1/(range(e))
		if (rho <= feasible[1] || rho >= feasible[2])
			stop(paste("Rho", rho, "outside feasible range:",
                        paste(feasible, collapse=":")))
	}
        if (method == "chol"){
	    if (x$style %in% c("W", "S") && !(spatialreg::can.be.simmed(x)))
			stop("Cholesky method requires symmetric weights")
		if (x$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(x$neighbours, x$weights)))
			stop("Cholesky method requires symmetric weights")
		if (x$style %in% c("W", "S")) {
			V <- listw2mat(listw2U(spatialreg::similar.listw(x)))
		}
		mat <- diag(n) - rho * V
		res <- chol2inv(chol(mat))
	} else if (method == "solve") {
		mat <- diag(n) - rho * V
		res <- solve(mat)
	} else stop("unknown method")
	attr(res, "call") <- match.call()
	res
}


