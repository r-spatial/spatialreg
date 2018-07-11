spreg<-function(formula, data=list(), listw, listw2=NULL, endog = NULL, instruments= NULL, lag.instr = FALSE, initial.value=0.2, model = c("sarar", "lag", "error", "ivhac", "ols"), het = FALSE, verbose=FALSE, na.action = na.fail,  HAC = FALSE, distance = NULL, type = c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth="variable" ,step1.c = FALSE, control = list()){

         		
#extract model objects	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action, method="model.frame")
	na.act<-attr(mf,'na.action')

# record call
cl<- match.call()


#generates x and y 
	y<-c(model.extract(mf,"response"))
	x<-model.matrix(mt,mf)

#checks on teh dimensions of x and y 	
if (length(y)!=nrow(x)) 
	stop("x and y have different length")

#check that X and y does not have missing values	
if (any(is.na(y))) 
        stop("NAs in dependent variable")
if (any(is.na(x))) 
        stop("NAs in independent variable")

if(HAC && model %in% c("lag","error","sarar")) stop("Model should be one of 'ivhac', or 'ols' when HAC is true ")

if(HAC){
 if(model != "ivhac" && !is.null(endog)) model <- 'ols.end'	
	if(is.null(distance)) stop("No distance measure specified")
	if(!inherits(distance,"distance")) 
	stop("The distance measure is not a distance object")
	
if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
}	
	

	
#fix the dimensions of the problem
	n<-nrow(x)
	k<-ncol(x)	
	xcolnames<-colnames(x)

	K<-ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
	
	
#check that W is an object of class listw or a Matrix 

if(!(model %in% c("ols", "ols.end"))){
	
if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
if(inherits(listw,"Matrix"))  Ws <- listw	


#check on the dimensions of x and W	
if (nrow(x) != nrow(Ws))
	stop("Input data and weights have different dimension")

if (k > 1) {
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
				         }
            wwx <- Ws %*% wx                    					         
    }


if(!is.null(listw2) && model != "sarar") stop("listw2 can be specified only with sarar")

if(is.null(listw2)) {
	twow <- FALSE		
	Ws2 <- Ws
	}
	
else{ 

twow <- TRUE	

if(!inherits(listw2,c("listw", "Matrix", "matrix"))) stop("listw2 format unknown")
if(inherits(listw2,"listw"))  Ws2<-listw2dgCMatrix(listw2)	
if(inherits(listw2,"matrix"))  Ws2<-Matrix(listw2)	

	 # # feeedback from user
if(identical(listw, listw2)){ 
	twow <- FALSE		
	Ws2 <- Ws
	}

	
if (k > 1) {
        w2x <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          W2x <- Ws2 %*% x[, i]          
            w2x[, (i - (K - 1))] <- as.matrix(W2x)
				         }
		  w2wx <- Ws2 %*% wx                   	
          w2wwx <- Ws2 %*% wwx                    	          

    }
	}

 wy<-Ws %*% y	
 colnames(wy)<-"lambda"

}


### Definition of the instruments for all cases: if there is an endogenous variable the instruments have to be specified in all models. 
if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")

if(model %in% c("sarar","lag", "ivhac")){
 	if (!is.null(endog)) {
		endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
			if(!is.null(instruments)){
					instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
					if(lag.instr) {
				         winst <- Ws %*% instruments
            			 wwinst<- Ws %*% winst	
								if(twow){
										w2i <- Ws2 %*% instruments 
										 w2wi <- Ws2 %*% winst 
										 w2wwi <- Ws2 %*% wwinst 	
AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst), as.matrix(w2i), as.matrix(w2wi),as.matrix(w2wwi))        
										}
								else  AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
									}
					else  AddH <- instruments        
if (K==2) {
	if(twow) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
	else  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), AddH)
		}
else {
	if(twow) Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
	else  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), AddH)
	}
									}
		 	# else stop("Instruments should be specified if there is an endogenous variable")
	Zmat<- cbind(x, endog, as.matrix(wy))            
	colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
}
	else {
	Zmat<- cbind(x, as.matrix(wy))                    
if (K==2){
	if(twow) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
	else  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx)) 
}
else {
	if(twow) Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
	else  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx)) 
	}
    }
    }

if(model == "error" ){
	if (!is.null(endog)) {
          endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
		  instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
	 if(lag.instr) { 
	        winst <- Ws %*% instruments           
	        AddH<- cbind(instruments, as.matrix(winst))
}
else AddH<- cbind(instruments)

Zmat<- cbind(x, endog)            
colnames(Zmat) <- c(colnames(x), colnames(endog)) 

if (K==2) Hmat<-cbind(x, wx, AddH) 
else Hmat<-cbind(1, x, wx, AddH)

 }
 
else {

if (K==2) Hmat<-cbind(x,wx)
else Hmat<-cbind(1, x,wx)	
	Zmat<- x
	}
}

if(model == "ols.end" ){
	if (!is.null(endog)) {
			endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))			
			instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))	
AddH<- cbind(instruments)
Zmat<- cbind(x, endog)            
colnames(Zmat) <- c(colnames(x), colnames(endog)) 
if (K==2) Hmat<-cbind(x, AddH) 
else Hmat<-cbind(1, x, AddH)
 }
}
		
	
if(model %in% c("sarar","error")){


firststep<-spatial.ivreg(y = y , Zmat = Zmat, Hmat = Hmat, het = het, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)
ubase<-residuals(firststep)


if (initial.value=="SAR"){
		Wubase<-Ws2 %*% ubase
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value


if(het){
	
Ggmat<-gg_het(Ws2, ubase, n)

optres <-nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control= control, v = Ggmat, verbose = verbose)

#list(abs.tol = abs.tol, rel.tol = rel.tol)

rhotilde<-optres$par
 # print(rhotilde)

if(step1.c){
 gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
 # print(gmm.weghts1.c$Phiinv)
# gmm.weghts1.c<-psirhorho_het_mod(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)

optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts1.c$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control)	

rhotilde<-optres$par
gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)

# gmm.weghts1.c<-psirhorho_het_mod(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)

vcmat_2sls <- Omega_het(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1, gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1, gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG, gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon, gmm.weghts1.c$Zstar, Ws2, step1.c = TRUE)


coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
rownames(coeff_2sls)<-c(colnames(Zmat), 'rho')
s2_2sls<-crossprod(ubase)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gm spatial"

k<-nrow(coeff_2sls)
R<-matrix(0,1,k)
R[,((k-1):k)]<-1
Rbeta<-R%*%coeff_2sls
Rvar<-R%*% vcmat_2sls$Omega %*%t(R)
stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
pval <- pchisq(stat,df=1,lower.tail=FALSE)
W<-list(stat=stat,pval=pval)



 results_2sls <- list(coefficients=coeff_2sls,var=vcmat_2sls$Omega, s2=s2_2sls, call=cl, residuals=as.numeric(ubase), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde)
 
 class(results_2sls)<-c("sphet", "gstsls")
	
}


}

else{
	
Ggmat<-gg_hom(Ws2, ubase, n)
optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control = control, v = Ggmat, verbose = verbose)
rhotilde<-optres$par
# print(rhotilde)	
	}
	




yt  <- y - rhotilde * Ws2 %*% y
wZmat <- Ws2 %*% Zmat
Zt <- Zmat - rhotilde * wZmat

# if(!sarar && is.matrix(endog)) Hmat <- cbind(x, wx, instruments)	
# else Hmat<- cbind(x,wx)	


# print(Hmat)
    secondstep<-spatial.ivreg(y =yt , Zmat = Zt, Hmat = Hmat, het = het, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)
delta <- coefficients(secondstep)
utildeb <- y - Zmat %*% delta
  # print(delta)

if(het){


Ggmat<-gg_het(Ws2, utildeb, n)

  gmm.weghts<-psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)

 # gmm.weghts<-psirhorho_het_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)

# print(gmm.weghts$Phiinv)
optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control)	

rhofin<-optres$par
	 gmm.weghts<-psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)

 # gmm.weghts<-psirhorho_het_mod(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)

 vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar, Ws2, step1.c = FALSE)
# vcmat <- Omega_het_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar, Ws2, step1.c = FALSE)
}

else{
	
Ggmat<-gg_hom(Ws2, utildeb, n)

   gmm.weghts<-psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )

   # gmm.weghts<-psirhorho_hom_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )


optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control )	

rhofin<-optres$par
 # print(rhofin)
     gmm.weghts<-psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
     # gmm.weghts<-psirhorho_hom_mod(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )

      vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)

     # vcmat <- Omega_hom_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2,gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)

	}





coeff <- as.matrix(c(as.numeric(delta), rhofin))
rownames(coeff)<-c(colnames(Zmat), 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gm spatial"

k<-nrow(coeff)
R<-matrix(0,1,k)
R[,((k-1):k)]<-1
Rbeta<-R%*%coeff
Rvar<-R%*% vcmat$Omega %*%t(R)
stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
pval <- pchisq(stat,df=1,lower.tail=FALSE)
W<-list(stat=stat,pval=pval)



if(het && step1.c) results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde,  twosls = results_2sls)

else  results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde)

 
 class(results)<-c("sphet", "gstsls")
 
 }

 
 if(model %in% c("lag", "ivhac", "ols.end")){
 	
results <-spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat, het = het, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)	
     model.data <- data.frame(cbind(y, x[, -1]))
    results$call <- cl
    results$model <- model.data
    results$type <- type
    results$bandwidth <- bandwidth
    results$method <- "s2slshac"
    results$HAC <- HAC
    class(results) <- c("sphet", "stsls_sphet")
 
 	}
 	
if(model == "ols")	{
	results <-hac.ols(y =y , x = x, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)	
     model.data <- data.frame(cbind(y, x[, -1]))
    results$call <- cl
    results$model <- model.data
    results$type <- type
    results$bandwidth <- bandwidth
    results$method <- "olshac"
    results$HAC <- HAC
    class(results) <- c("sphet")

}


 return(results)
}


impacts.gstsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
    tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
coefs <- drop(obj$coefficients)

    p2 <- length(coefs)
    rho <- coefs[(p2-1)]
    beta <- coefs[1:(p2-2)]
    p <- length(beta)
	p1 <- p + 1
    lambda <- coefs[p2]
# rho is lag coef., lambda is error coef (reversed from function use)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigmawor <- obj$var[-p2,-p2]
	Sigma <- Sigmawor[c(p1, (1:(p1-1))), c(p1, (1:(p1-1)))]    
	irho <- 1
    drop2beta <- 1
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, tol=tol,
        empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}






impacts.stsls_sphet <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
    tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    rho <- coefs[(p2)]
    beta <- coefs[1:(p2-1)]
# rho is lag coef., lambda is error coef (reversed from function use)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigma <- obj$var[c(p2, (1:(p2-1))), c(p2, (1:(p2-1)))]
    irho <- 1
    drop2beta <- 1
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, tol=tol,
        empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}
