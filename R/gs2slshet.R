gstslshet<-function(formula, data=list(),listw, na.action=na.fail, zero.policy=NULL, initial.value=0.2, abs.tol=1e-20, rel.tol=1e-10, eps=1e-5, inverse=T, sarar=T){

##functions that need to be sourced
	#source("twostagels.R")
	#source("utilities.R")
	#source("listw2dgCMatrix.R")
	#source("Omega.R")
	
	
	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
	
#extract model objects	
	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action, method="model.frame")
	na.act<-attr(mf,'na.action')

# record call
cl<- match.call()


#if(!inverse) warning("Approximated inverse could be inaccurate")

#preferences on missings values
if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }

#check that W ius an object of class listw
if(!inherits(listw,"listw")) 
	stop("The weights matrix is not a listw object")

#generates x and y 
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)

#checks on teh dimensions of x and y 	
if (length(y)!=nrow(x)) 
	stop("x and y have different length")

#check on the dimensions of x and W	
if (nrow(x) != length(listw$neighbours))
	stop("Input data and weights have different dimension")

#check that X and y does not have missing values	
if (any(is.na(y))) 
        stop("NAs in dependent variable")
if (any(is.na(x))) 
        stop("NAs in independent variable")

	
#fix the dimensions of the problem
	n<-nrow(x)
	k<-ncol(x)	
	xcolnames<-colnames(x)
	K<-ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)

if(sarar){	
	
	wy<-lag.listw(listw,y, zero.policy=zero.policy)
	wy<-array(wy,c(length(y),1))
	colnames(wy)<-("Wy")
	
	
if (any(is.na(wy)))
	stop("NAs in spatially lagged dependent variable")
	
if (k > 1) {
        WX <- matrix(nrow = n, ncol = (k  - (K - 1)))
        WWX <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
            wx <- lag.listw(listw, x[, i], zero.policy = zero.policy)
            wwx<- lag.listw(listw, wx, zero.policy = zero.policy)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WX[, (i - (K - 1))] <- wx
				 WWX[, (i - (K - 1))] <- wwx
				         }
    }

instr<-cbind(WX,WWX) 




##spatial two stage least square of the initial model
firststep<-tsls(y=y,yend=wy, X=x, Zinst = instr)

ubase<-residuals(firststep)


#GM step 1b
		int1<-Ggfastfast(listw=listw,ubase,n, zero.policy = zero.policy)

##One could start from this initial value but if uses optimize, the initial values are not an issue
if (initial.value=="SAR"){
		Wubase<-lag.listw(listw,ubase, zero.policy = zero.policy)
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value
		
		
optres1<-nlminb(pars, arg, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1)

		param<-optres1$par
		# print(param)
		reg<-x
		u<-ubase
		toinst<-wy
		
		fi<-fifour(x,listw,instr,ubase,wy,param,n,inverse,eps, zero.policy= zero.policy)
		# print(fi$nlw)

GMMfeas3<-nlminb(param, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1, VC=fi$nlw)


rhotilde <- GMMfeas3$par
# print(rhotilde)

fi1<-fifour(x,listw,instr,ubase,wy,rhotilde,n,inverse=inverse,eps=eps, zero.policy = zero.policy)
Om<-Omega(n,gamma=fi1$res1,H=fi1$instr,param=rhotilde,G=int1$bigG,FIinv=fi1$nlw,a=fi1$a,FI=fi1$nl,P=fi1$P,gammas=fi1$gammas,Ws=fi1$Ws)

yt  <- y - rhotilde * wy
xt <- x - rhotilde * lag.listw(listw,x, zero.policy = zero.policy)
wyt<-wy - rhotilde * lag.listw(listw,wy, zero.policy = zero.policy)
colnames(xt)<-xcolnames
colnames(wyt)<-c('Wyt')

secstepb<-tsls(y=yt, yend=wyt, X=xt, Zinst = instr, reg=x, end=wy, yor=y,modified=TRUE)

utildeb<-secstepb$residuals
int1b<-Ggfastfast(listw, utildeb,n, zero.policy = zero.policy)
psippb<-fistslsfast(reg=xt, Ws=fi$Ws, instr=instr, resid=utildeb, toinst=wyt, param=rhotilde,solo=x,end=wy)
pars1<-rhotilde


GMMstsls1b<-nlminb(pars1, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1b, VC=psippb$nlw)

rhohat<-GMMstsls1b$par
#xt <- x - rhohat * lag.listw(listw,x)
psifinal<-fistslsfast(reg=xt, Ws=fi$Ws, instr=instr, resid=utildeb, toinst=wyt, param=rhohat,solo=x,end=wy)
Omfinal<- Omegabis(gammas=psifinal$res1, Hs=psifinal$instr, param=rhohat, G=int1b$bigG, FIinv=psifinal$nlw, as=psifinal$a, n, FI=psifinal$nl, P=psifinal$P) 
coef<-(secstepb$coefficients[2:(ncol(xt)+1)])
coeff<-as.matrix(c(coef, secstepb$coefficients[1], rhohat ))
rownames(coeff)<-c(xcolnames, 'lambda', 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gs2slshac"

	k<-nrow(coeff)
	R<-matrix(0,1,k)
	R[,((k-1):k)]<-1
	Rbeta<-R%*%coeff
	Rvar<-R%*%Omfinal$VarCov%*%t(R)
	stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
	pval <- pchisq(stat,df=1,lower.tail=FALSE)
W<-list(stat=stat,pval=pval)



results<-list(coefficients=coeff,var=Omfinal$VarCov, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W,yhat=secstepb$yhat )
}
else{
	##this is an error model estimated by ols
firststep<-lm(y~x-1)
ubase<-residuals(firststep)

int1<-Ggfastfast(listw=listw,ubase,n, zero.policy = zero.policy)

##One could start from this initial value but if uses optimize, the initial values are not an issue
if (initial.value=="SAR"){
		Wubase<-lag.listw(listw,ubase, zero.policy = zero.policy)
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value
		
		
optres1<-nlminb(pars, arg, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1)

		param<-optres1$par
		reg<-x
		u<-ubase
				
		fi<-fierror(x,listw, ubase, param,n,inverse=inverse,eps=eps)

GMMfeas3<-nlminb(param, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1, VC=fi$nlw)

rhotilde<- GMMfeas3$par
# print(rhotilde)

fi1<-fierror(x,listw,ubase,rhotilde,n,inverse=inverse,eps=eps)
Om<-Omega(n,gamma=fi1$res1, H=x, param=rhotilde,G=int1$bigG,FIinv=fi1$nlw,a=fi1$a,FI=fi1$nl,P=fi1$P,gammas=fi1$gammas,Ws=fi1$Ws) ## in an error model H=X
#####################


yt  <- y - rhotilde * lag.listw(listw,y, zero.policy=zero.policy)
xt <- x - rhotilde * lag.listw(listw,x, zero.policy=zero.policy)
colnames(xt)<-xcolnames

secstepb<-lm(yt~xt-1)

utildeb<-y - x %*% as.matrix(coefficients(secstepb)) 

int1b<-Ggfastfast(listw, utildeb,n, zero.policy = zero.policy)

psippb<-fistslserror(reg=x, Ws=fi$Ws, resid=utildeb, param=rhotilde, solo = cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)))
pars1<-rhotilde


GMMstsls1b<-nlminb(pars1, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1b, VC=psippb$nlw)

rhohat<-GMMstsls1b$par
psifinal<-fistslserror(reg=x, Ws=fi$Ws, resid=utildeb, param=rhohat,solo=cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)))
Omfinal<- Omegabis(gammas=psifinal$res1, Hs = cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)) , param=rhohat, G=int1b$bigG, FIinv=psifinal$nlw, as=psifinal$a, n, FI=psifinal$nl, P=psifinal$P) 
coef<-coefficients(secstepb)
coeff<-as.matrix(c(coef, rhohat ))
rownames(coeff)<-c(xcolnames, 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gs2slshac"


results<-list(coefficients=coeff,var=Omfinal$VarCov, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=NULL,yhat=fitted(secstepb) )

	
	}

class(results)<-c("sphet", "gstsls")

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
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(beta, rho, lambda)
    Sigma <- obj$var
    irho <- p2-1
    drop2beta <- c(p2-1, p2)
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, tol=tol,
        empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}


