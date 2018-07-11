listw2dgCMatrix<-function (listw, zero.policy=NULL) 
{
    if (!inherits(listw, "listw")) 
        stop("not a listw object")
    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    n <- length(listw$neighbours)
    cardw <- card(listw$neighbours)
    p0 <- as.integer(c(0, cumsum(cardw)))
    scard <- sum(cardw)
    t<-unlist(listw$neighbours)
    if (zero.policy) t <- t[t > 0]
    t<-t-1
    res <- new("dgCMatrix", i = as.integer(t), p = p0,  Dim = as.integer(c(n,n)), x = unlist(listw$weights))
    res<-t(res)
}
