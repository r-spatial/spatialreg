
arg<-function(rhopar, v) 
{
    sys <- v$litg -v$bigG %*% c(rhopar, rhopar^2)
    value <- sum(sys^2)

value
}


arg1<-function(rhopar, v, VC)
{
    sys <- v$litg - v$bigG %*% c(rhopar, rhopar^2) 
    value <- t(sys) %*%  VC  %*% sys
    value
}


fifour<-function(reg, listw, instr, u,  toinst, param, n, inverse, eps, zero.policy = NULL) 
{
			
	H<-cbind(reg, instr) #elements to generate P
	Z<-cbind(reg, toinst)

	k<-dim(Z)[2]
	HH<-crossprod(H)
	HZ<-crossprod(H,Z)
	ZH<-crossprod(Z,H)
	HHn<-HH/n
	HZn<-HZ/n
	ZHn<-ZH/n
	HHninv<-solve(HHn)
	sec<-ZHn %*% HHninv %*% HZn
	secinv<-solve(sec)
	P<-HHninv %*% HZn %*% secinv
   Ws<-listw2dgCMatrix(listw, zero.policy = zero.policy)
    Wst<-t(Ws)
    WspWs<-crossprod(Ws)
    diag(WspWs)<-0
    Zt<-t(Z)
    ZtWst<-Zt%*%Wst
    ZprWsp<-Zt-param*ZtWst
    Wu<-Ws%*%u
    urWu<-u-param*Wu
    A1<-WspWs
    A2A2t<-Ws+Wst
    A1t<-t(A1)
    A1A1t<-A1+A1t
    alpha1<-(-1)*(ZprWsp%*%A1A1t%*%urWu)/n
    alpha2<-(-1)*(ZprWsp%*%A2A2t%*%urWu)/n
    gamma<-urWu^2 ##traces of the elements of psi
#	gammas<-as(Diagonal(,as.vector(gamma)),"sparseMatrix")
	gammas<-.symDiagonal(n, x = as.numeric(gamma), uplo = "U")
#	gammas<-Diagonal(x = as.numeric(gamma))
   tr1<-sum(diag(A1A1t%*%gammas%*%A1A1t%*%gammas))/2
   tr2<-sum(diag(A2A2t%*%gammas%*%A2A2t%*%gammas))/2
   tr3<-sum(diag(A1A1t%*%gammas%*%A2A2t%*%gammas))/2
   a1p<-H%*%P%*%alpha1
	a2p<-H%*%P%*%alpha2
	
	#a1<-matrix(,nrow=n,ncol=1)
	#a2<-matrix(,nrow=n,ncol=1)
		
		if(inverse){	
			IrWp<- -param*Wst
			diag(IrWp)<-1
			IrWpi<-solve(IrWp)
			}
			else{
				eps<-eps
				i<-2
				pWst<- param*Wst
				IrWpi<- pWst
				
				while(abs(sum(pWst))>eps){
					pWst<- param^i*Wst^i
					IrWpi<- IrWpi + pWst
					i<-i+1
					pWst<- param^i*Wst^i
								}
								
			diag(IrWpi)<-1
				}
				
			a1<-IrWpi%*%a1p
			a2<-IrWpi%*%a2p

			a<-cbind(as.matrix(a1),as.matrix(a2))
			a1gammaa1v<-sum(a1^2*gamma)
			a2gammaa2v<-sum(a2*gamma*a2)
			a2gammaa1v<-sum(a2*gamma*a1)
	      fi11<- (as.matrix(tr1) + a1gammaa1v)/n
          fi21<- (as.matrix(tr3) + a2gammaa1v)/n
		   fi22<- (as.matrix(tr2) + a2gammaa2v)/n
		   FIr1<- cbind(fi11, fi21)
		   FIr2<- cbind(fi21, fi22)
		    FI<-rbind(FIr1, FIr2)
		    FIinv<-solve(FI)		    
list(nlw=FIinv, a=a, res1=gamma,instr=H, dim=n, nl=FI, Ws=Ws, P=P,gammas=gammas)
}

fistslsfast<-function(reg, Ws, instr, resid, toinst, param, solo, end){
	n<-length(resid)
	H <- cbind(solo, instr)
	Z1<-cbind(solo,end)
	Z<-cbind(reg,toinst)
	HH<-crossprod(H)
	HZ<-crossprod(H,Z)
	ZH<-crossprod(Z,H)
	HZn<-HZ/n
	ZHn<-ZH/n
	HHninv<-solve(HH/n)
	sec<-ZHn %*% HHninv %*% HZn
	secinv<-solve(sec)
	P<-HHninv %*% HZn %*% secinv
	Wst<-t(Ws)
	A1<- Wst %*% Ws
	diag(A1)<- 0
	A1t<- t(A1)
	A1A1<-A1+A1t
	A2A2 <- Ws + Wst
	Z1pWp<-t(Z1)%*%Wst
	Z1prZ1pWp<-t(Z1)-param*Z1pWp
	Wu<-Ws%*%resid
	IrWu<-resid-param*Wu
	alpha1<-((-1)* Z1prZ1pWp %*%A1A1%*%IrWu)/n
	alpha2<-((-1)* Z1prZ1pWp %*%A2A2%*%IrWu)/n
	a1<-as.matrix(H%*%P%*%alpha1)
	a2<-as.matrix(H%*%P%*%alpha2)
	a<-as.matrix(cbind(as.matrix(a1),as.matrix(a2)))
	gamma<-IrWu^2
	gammas<-as(Diagonal(,as.vector(gamma)),"sparseMatrix")
	a1Ga1<-t(a1) %*% gammas%*% a1
	a2Ga2<-t(a2) %*% gammas%*% a2
	a1Ga2<-t(a1) %*% gammas%*% a2
   tr1<-sum(diag(A1A1%*%gammas%*%A1A1%*%gammas))/2
   tr2<-sum(diag(A2A2%*%gammas%*%A2A2%*%gammas))/2
   tr3<-sum(diag(A1A1%*%gammas%*%A2A2%*%gammas))/2
   fi11<- as.matrix((tr1 + a1Ga1)/n)
	fi22<- as.matrix((tr2 + a2Ga2)/n)
	fi12<- as.matrix((tr3 + a1Ga2)/n)
    FIr1<- cbind(fi11, fi12)
	FIr2<- cbind(fi12, fi22)
	FI<-rbind(as.matrix(FIr1), as.matrix(FIr2))
	FIinv<-solve(as.matrix(FI))
	list(nlw=as.matrix(FIinv), a=a, res1=gammas, instr=H, dim=n, nl=FI, P=P)
	}


Ggfastfast<-function(listw, u, n, zero.policy = NULL) 
{
     ub<-lag.listw(listw,u, zero.policy = zero.policy)
     ubb<-lag.listw(listw,ub, zero.policy = zero.policy)
     diag1 <- ub * u
     diag2<-ub^2
     Ws<-listw2dgCMatrix(listw, zero.policy = zero.policy)
     Wst<-t(Ws)
     Diag1<-as(Diagonal(,as.vector(diag1)),"sparseMatrix")
	  tra<-sum(diag(Ws%*%Diag1%*%Wst))
	  ubbub<-crossprod(ubb,ub)
 	  first<- (ubbub - tra)
 	  Diag2<-as(Diagonal(,as.vector(diag2)),"sparseMatrix")
	  tttt<-sum(diag(Ws%*%Diag2%*%Wst))
	  ubbubb<-crossprod(ubb,ubb)
     second<- (ubbubb - tttt)
#    second<- (ubbubb + tttt)
     uubb<-crossprod(u,ubb)
     ubub<-crossprod(ub,ub)
     third<-uubb + ubub  
     ububb<-crossprod(ub,ubb)
     bigG <- matrix(0, 2, 2)
     bigG[, 1] <- c(2*first, third)/n
     bigG[, 2] <- -c(second, ububb)/n
     diag3 <- u^2 
 	  Diag3<-as(Diagonal(,as.vector(diag3)),"sparseMatrix")
	  aa<-sum(diag(Ws%*%Diag3%*%Wst))
	  gamma1<-(ubub - aa)
	  uub<-crossprod(u,ub)
	  litg <- c(gamma1, uub)/n
     list(bigG = bigG, litg = litg)
}


fierror<-function(reg, listw, u, param, n,inverse,eps, zero.policy = NULL) 
{
			
	P<-solve(crossprod(reg)/n)
	
   Ws<-listw2dgCMatrix(listw, zero.policy = zero.policy)
    Wst<-t(Ws)
    WspWs<-crossprod(Ws)
    diag(WspWs)<-0
    Xt<-t(reg)
    XtWst<-Xt%*%Wst
    XprWsp<-Xt-param*XtWst
    Wu<-Ws%*%u
    urWu<-u-param*Wu
    A1<-WspWs
    A2A2<-Wst+Ws
    A2A2t<-t(A2A2)
    A1t<-t(A1)
    A1A1t<-A1+A1t
    alpha1<-(-1)*(XprWsp%*%A1A1t%*%urWu)/n
    alpha2<-(-1)*(XprWsp%*%A2A2t%*%urWu)/n
    gamma<-urWu^2 ##traces of the elements of psi
	gammas<-as(Diagonal(,as.vector(gamma)),"sparseMatrix")
   tr1<-sum(diag(A1A1t%*%gammas%*%A1A1t%*%gammas))/2
   tr2<-sum(diag(A2A2t%*%gammas%*%A2A2t%*%gammas))/2
   tr3<-sum(diag(A1A1t%*%gammas%*%A2A2t%*%gammas))/2
   a1p<-reg%*%P%*%alpha1
	a2p<-reg%*%P%*%alpha2
	a1<-matrix(,nrow=n,ncol=1)
	a2<-matrix(,nrow=n,ncol=1)
		if(inverse){	
	IrWp<- -param*Wst
			diag(IrWp)<-1
			IrWpi<-solve(IrWp)
			}
			else{
				eps<-eps
				i<-2
				pWst<- param*Wst
				IrWpi<- pWst
				
				while(abs(sum(pWst))>eps){
					pWst<- param^i*Wst^i
					IrWpi<- IrWpi + pWst
					i<-i+1
					pWst<- param^i*Wst^i
			}
			diag(IrWpi	)<-1
				}
			a1<-IrWpi%*%a1p
			a2<-IrWpi%*%a2p

#    	    print(a1)
			a<-cbind(as.matrix(a1),as.matrix(a2))
			a1gammaa1v<-sum(a1^2*gamma)
			a2gammaa2v<-sum(a2*gamma*a2)
			a2gammaa1v<-sum(a2*gamma*a1)
	        fi11<- (as.matrix(tr1) + a1gammaa1v)/n
           fi21<- (as.matrix(tr3) + a2gammaa1v)/n
		    fi22<- (as.matrix(tr2) + a2gammaa2v)/n
		    FIr1<- cbind(fi11, fi21)
		    FIr2<- cbind(fi21, fi22)
		    FI<-rbind(FIr1, FIr2)
		    FIinv<-solve(FI)
# print(FIinv)		    
list(nlw=FIinv, a=a, res1=gamma, dim=n, nl=FI, Ws=Ws, P=P,gammas=gammas)
}



fistslserror<-function(reg, Ws, resid, param, solo){
	 n<-length(resid)
		  H <- solo
		  Z1<-reg
 		 Z <- reg - param *Ws %*% Z1
 

		 HH<-crossprod(H)
		 HZ<-crossprod(H,Z)
		 ZH<-crossprod(Z,H)
		 HZn<-HZ/n
		 ZHn<-ZH/n
		 HHninv<-solve(HH/n)
		 sec<-ZHn %*% HHninv %*% HZn
		 secinv<-solve(sec)
	     P<-HHninv %*% HZn %*% secinv
 # P<-solve(crossprod(solo)/n)

	Wst<-t(Ws)
	A1<- Wst %*% Ws
	diag(A1)<- 0
	A1t<- t(A1)
	A1A1<-A1+A1t
	A2A2 <- Ws + Wst
	Z1pWp<-t(Z1)%*%Wst
	Z1prZ1pWp<-t(Z1)-param*Z1pWp
	Wu<-Ws%*%resid
	IrWu<-resid-param*Wu
	alpha1<-((-1)* Z1prZ1pWp %*%A1A1%*%IrWu)/n
	alpha2<-((-1)* Z1prZ1pWp %*%A2A2%*%IrWu)/n
	a1<-as.matrix(solo %*%P%*%alpha1)
	a2<-as.matrix(solo %*%P%*%alpha2)
	a<-as.matrix(cbind(as.matrix(a1),as.matrix(a2)))
	gamma<-IrWu^2
	gammas<-as(Diagonal(,as.vector(gamma)),"sparseMatrix")
	a1Ga1<-t(a1) %*% gammas%*% a1
	a2Ga2<-t(a2) %*% gammas%*% a2
	a1Ga2<-t(a1) %*% gammas%*% a2
   tr1<-sum(diag(A1A1%*%gammas%*%A1A1%*%gammas))/2
   tr2<-sum(diag(A2A2%*%gammas%*%A2A2%*%gammas))/2
   tr3<-sum(diag(A1A1%*%gammas%*%A2A2%*%gammas))/2
   fi11<- as.matrix((tr1 + a1Ga1)/n)
	fi22<- as.matrix((tr2 + a2Ga2)/n)
	fi12<- as.matrix((tr3 + a1Ga2)/n)
    FIr1<- cbind(fi11, fi12)
	FIr2<- cbind(fi12, fi22)
	FI<-rbind(as.matrix(FIr1), as.matrix(FIr2))
	FIinv<-solve(as.matrix(FI))
	 # print(FIinv)
	list(nlw=as.matrix(FIinv), a=a, res1=gammas, dim=n, nl=FI, P=P)
	}
