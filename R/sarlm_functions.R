# Copyright 2009-2013 by Roger Bivand

sar_error_hess_sse <- function(lambda, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml1_sse_env", env, lambda, beta, PACKAGE="spreg")
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


sar_lag_hess_sse <- function(rho, beta, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spreg")
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


