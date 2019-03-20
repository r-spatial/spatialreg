# Copyright 2005-2012 by Roger Bivand

# Simultaneous autoregressive
SAR <- function(IlW, weights) {
    t(IlW) %*% weights %*% IlW
}

# Conditional  autoregressive
CAR <- function(IlW, weights) {
    IlW %*% weights
}

# Spatial moving average
SMA <- function(IlW, weights) {
    IlW <- solve(IlW)
    t(IlW) %*% weights %*% IlW
}


f_spautolm_hess_nlm <- function(coefs, env) {
    ret <- f_spautolm_hess(coefs, env)
    -ret
}

f_spautolm_hess <- function(coefs, env) {
    lambda <- coefs[1]
    int <- get("interval", envir=env)
    if (lambda <= int[1] || lambda >= int[2]) return(-Inf)
    beta <- coefs[-1]
    X <- get("X", envir=env)
    Y <- get("Y", envir=env)
    fitted <- X %*% beta
    residuals <- Y - fitted
    dmmf <- eval(parse(text=get("family", envir=env)))
    if (get("family", envir=env) == "SMA") IlW <- dmmf((get("I", envir=env) + 
        lambda * get("W", envir=env)), get("Sweights", envir=env))
    else IlW <- dmmf((get("I", envir=env) - lambda * get("W", envir=env)), 
        get("Sweights", envir=env))
    SSE <- c(crossprod(residuals, as.matrix(IlW %*% residuals)))
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    det <- ifelse(get("family", envir=env) == "CAR", 0.5*ldet, ldet)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env))
        cat("lambda:", lambda, "function:", ret, "Jacobian", ldet, "SSE",
            SSE, "\n")
    if (!is.finite(ret)) return(-Inf)
    ret
}


