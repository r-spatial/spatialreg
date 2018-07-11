tsls <- function(y,yend,X,Zinst,end, reg, yor, modified=FALSE, HAC=FALSE, distance=distance,  type=c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth=bandwidth, zero.policy=zero.policy) {
### estimation engine that deals with the various cases
#silent<-set.VerboseOption(FALSE)
if(modified){
	H <- cbind(reg, Zinst)
	Z <- cbind(yend, X)
	HH <- crossprod(H,H)
	Hye <- crossprod(H,Z)
	bz <- solve(HH,Hye)
	Zp <- H %*% bz
	ZpZp <- crossprod(Zp,Z)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	delta <- crossprod(ZpZpi,Zpy)
	names(delta)<-c(colnames(yend), colnames(X))
	Z1<-cbind(end, reg)
	yp <- Z1 %*% delta
    e <- yor - yp
     result <- list(coefficients=delta, yhat=yp, residuals=e)
	}
else{
		K<-ifelse(colnames(X)[1] == "(Intercept)" || all(X[,1]==1), 2, 1)
	#print(K)
	#if x does not have an intercept, one should add it to the instruments anyway
if (K==1)	H <- cbind(1,X,Zinst)
else	H <- cbind(X,Zinst)
#print(H[1,])
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
	HH <- crossprod(H,H)
	Hye <- crossprod(H,yend)
	bz <- solve(HH,Hye)
	yendp <- H %*% bz
	Zp <- cbind(yendp,X)
#	print(Zp[1,])
	ZpZp <- crossprod(Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	delta <- crossprod(ZpZpi,Zpy)
	names(delta)<-c(colnames(yend), colnames(X))
	yp <- Z %*% delta
	e <- y - yp
if(HAC){
	n<-nrow(X)
	#print(n)
		Ker<-vector(mode='list',length=n)
	#	print(Ker)
	#	print(match.arg(type))
		ker.fun<-switch(match.arg(type), Triangular={
			triangular
			}, Epanechnikov={
				epan
				}, Bisquare = {
					bisq
					}, TH={
						th
						}, QS = {
							qs
							}, Parzen = {
								parzen
								},
								Rectangular = {
									rectangular
								})
#					print(bandwith)
# print(ker.fun)
#print(is.numeric(bandwith))
if(is.null(attributes(distance)$GeoDa$dist)){
	Ker<-lapply(distance$weights,ker.fun, bandwidth=bandwidth)
	# print(Ker)
	Kern<-nb2listw(distance$neighbours,style="B", glist=Ker, zero.policy=zero.policy)
	} 
else{
	Ker<-lapply(attributes(distance)$GeoDa$dist,ker.fun, bandwidth=bandwidth)
	Kern<-nb2listw(distance,style="B", glist=Ker, zero.policy=zero.policy)
	# print(Kern[[1]])
	} 
# print(Ker[[1]])

# print(Kern$weights[[1]])
He<-matrix(,dim(H)[1],dim(H)[2])
#KHpe<-matrix(,dim(H)[1],dim(H)[2])
for (i in 1:dim(H)[2]) He[,i]<- H[,i] * e
#for(j in 1:dim(H)[2]){
#	for (i in 1:n) KHpe[i,j]<- sum(Ker[[i]]*He[distance$neigh[[i]],j])  + He[i,j]
#	 }

KHpe<-lag.listw(Kern,He, zero.policy=zero.policy) +He
KHeHe<-(t(He) %*% KHpe)
#KHeHe<-crossprod(He, KHpe)
HHp<-solve(HH)
ZpH<-crossprod(Z,H)
HpZ<-crossprod(H,Z)
fp<-ZpZpi%*%ZpH%*%HHp
vardelta<- (fp%*%KHeHe%*%t(fp))
s2 <- crossprod(e,e) / df
#sedelta<-sqrt(diag(vardelta))
#tdelta<-delta/sedelta
#pdelta<-pnorm(abs(tdelta), lower.tail=FALSE )*2
	}

else{
    	s2 <- crossprod(e,e) / df
	    vardelta <- ZpZpi * as.numeric(s2)
#	    sedelta <- sqrt(diag(vardelta))
#	    tdelta <- delta / sedelta
#	    pdelta <- pnorm(abs(tdelta),lower.tail=FALSE) * 2
}

	    result <- list(coefficients=delta,var=vardelta,s2=s2,
	          residuals=as.numeric(e),yhat=yp)
}	          
	result
}
