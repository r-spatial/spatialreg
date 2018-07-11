optimfunct <- function (rhopar,v,verbose=FALSE){
	 vv <- v$bigG %*% c(rhopar[1], rhopar[1]^2) - v$litg
    value <- sum(vv^2)
    if (verbose) cat("function:", value, "lambda:", rhopar[1], "\n")
    value
}


optimfunct_eff <- function (rhopar,v, vcmat, verbose=FALSE){
	 vv <- v$bigG %*% c(rhopar[1], rhopar[1]^2) - v$litg
    value <- t(vv) %*% vcmat %*% vv
    if (verbose) cat("function:", value, "lambda:", rhopar[1], "\n")
    value
}




gg_het<-function(Ws, u, n) 
{
     ub<- Ws %*% u
     ubb<-Ws %*% ub
     diag1 <- ub * u
     diag2<-ub^2
     Wst<-t(Ws)
     # Diag1<-as(Diagonal(,as.vector(diag1)),"sparseMatrix")
     Diag1<-Matrix(diag(n), sparse = TRUE)
     diag(Diag1)<-diag1
	  tra<-sum(diag(Ws%*%Diag1%*%Wst))
	  ubbub<-crossprod(ubb,ub)
 	  first<- (ubbub - tra)
 	  # Diag2<-as(Diagonal(,as.vector(diag2)),"sparseMatrix")
     Diag2 <- Matrix(diag(n), sparse = TRUE)
     diag(Diag2) <- diag2
	  tttt<-sum(diag(Ws%*%Diag2%*%Wst))
	  ubbubb<-crossprod(ubb,ubb)
     second<- (ubbubb - tttt)
     uubb<-crossprod(u,ubb)
     ubub<-crossprod(ub,ub)
     third<-uubb + ubub  
     ububb<-crossprod(ub,ubb)
     bigG <- matrix(0, 2, 2)
     bigG[, 1] <- c(2*as.numeric(first), as.numeric(third))/n
     bigG[, 2] <- -c(as.numeric(second), as.numeric(ububb))/n
     diag3 <- u^2 
 	  # Diag3<-as(Diagonal(,as.vector(diag3)),"sparseMatrix")
     Diag3 <- Matrix(diag(n), sparse = TRUE)
     diag(Diag3) <- diag3
	  aa<-sum(diag(Ws%*%Diag3%*%Wst))
	  gamma1<-(ubub - aa)
	  uub<-crossprod(u,ub)
	  litg <- c(as.numeric(gamma1), as.numeric(uub))/n
     list(bigG = bigG, litg = litg)
}




gg_hom<-function (Ws, u, n) 
{
    wu <- Ws %*% u
    wwu<- Ws %*% wu

    trwpw <- sum(Ws^2) 
    v.vec <- 1/ (1 + (trwpw/n)^2)
    d <- (1/n)* trwpw
    du <- d*u
    ubub<- crossprod(wu)
    udu<- crossprod(u,du)
    uA1u<- v.vec * (ubub - udu)
	uA2u<- crossprod(u, wu)
		    
	ubbubb <- crossprod(wwu)    
    dwu <- d*wu
    ubdub <- crossprod(wu,dwu)
 	wuA1wu <-  v.vec * (ubbubb - ubdub)   
    wuA2wu <- crossprod(wu, wwu)

    udub <- crossprod(u,dwu)
    ububb <-crossprod(wu,wwu)	
    upA1A1wu   <- 2 * v.vec *(ububb - udub)
    upA2A2wu <- crossprod(u,wwu) +  ubub
    
    
    bigG <- matrix(0, 2, 2)
        bigG[, 1] <- c(as.numeric(upA1A1wu), as.numeric(upA2A2wu))/n
        bigG[, 2] <- -c(as.numeric(wuA1wu), as.numeric(wuA2wu))/n

    litg <- c(as.numeric(uA1u), as.numeric(uA2u))/n
    
list(bigG = bigG, litg = litg, trwpw = trwpw, wu = wu, wwu = wwu, d = d, v.vec = v.vec)
}

psirhorho_hom <- function(rho, residuals, Hmat, Zmat, Ws, d, v.vec){
	n <- length(residuals)  
	epsilon<- residuals - rho * Ws %*% residuals 
	mu3<- sum(epsilon^3) /n 
	mu4<- sum(epsilon^4) /n 	
	sigma2<-crossprod(epsilon)/n
	
	Zstar <- Zmat - rho * Ws %*% Zmat
	Qhz<-crossprod(Hmat,Zstar)/n
	Qhh<-crossprod(Hmat)/n
	Qhhi<-solve(Qhh)
	sec.part<- t(Qhz) %*% Qhhi %*% Qhz
	Pmat<-Qhhi %*% Qhz %*% solve(sec.part)
	Tmat<-Hmat %*% Pmat	
	
    trA2A2 <- sum(Ws^2) 
    trA2A2I<- Matrix(diag(n), sparse=TRUE)* (trA2A2/n)
	A1 <- v.vec * (crossprod(Ws) - trA2A2I)
	A2 <- Ws
	
	
	A1pluA1 <- A1 + t(A1)	
	A2pluA2 <- A2 + t(A2)
	A1pluA1eps<- A1pluA1 %*% epsilon
	A2pluA2eps<- A2pluA2 %*% epsilon

	ZA1pluA1eps<- crossprod(Zstar, A1pluA1eps)
	ZA2pluA2eps<- crossprod(Zstar, A2pluA2eps)
		
	alpha1 <- -1*ZA1pluA1eps/n
	alpha2 <- -1*ZA2pluA2eps/n
	
	a.vec1 <- Tmat %*% alpha1 
	a.vec2 <- Tmat %*% alpha2 
	

    tr22<-sum(diag(A2pluA2 %*% A2pluA2))/(2*n)
#   tr22<-(2*trA2A2 + 2* t(as.vector(t(Ws))) %*% as.vector(Ws)) /(2*n)
    tr11<-sum(diag(A1pluA1 %*% A1pluA1))/(2*n)
	tr12<-sum(diag(A1pluA1 %*% A2pluA2))/(2*n)

	a1.a1 <-crossprod(a.vec1) /n 
	a2.a2 <-crossprod(a.vec2) /n 
	a1.a2 <-crossprod(a.vec1, a.vec2) /n 

    vecA1<- diag(A1)
    vecA2<- diag(A2)

    vecA1vecA1<- crossprod(vecA1)/n
    vecA2vecA2<- crossprod(vecA2)/n
    vecA1vecA2<- crossprod(vecA1,vecA2)/n 

     a1.vecA1<- crossprod(a.vec1, vecA1)
     a1.vecA2<- crossprod(a.vec1, vecA2)
     a2.vecA2<- crossprod(a.vec2, vecA2)
     a2.vecA1<- crossprod(a.vec2, vecA1)    
    
    la.term11 <- (a1.vecA1 + a1.vecA1)/n
    la.term12 <- (a1.vecA2 + a2.vecA1)/n
    la.term22 <- (a2.vecA2 + a2.vecA2)/n
    
    
     effe<- mu4-3*sigma2^2
    
    phi11<- as.numeric(sigma2^2 * tr11 + sigma2 * a1.a1 + effe *    vecA1vecA1 + mu3 * la.term11)
    phi12<- as.numeric(sigma2^2 * tr12 + sigma2 * a1.a2 + effe *    vecA1vecA2 + mu3 * la.term12)
    phi21<- as.numeric(sigma2^2 * tr12 + sigma2 * a1.a2 + effe *    vecA1vecA2 + mu3 * la.term12)
    phi22<- as.numeric(sigma2^2 * tr22 + sigma2 * a2.a2 + effe *    vecA2vecA2 + mu3 * la.term22)
    

     Phi<-cbind(rbind(phi11,phi21),rbind(phi12,phi22))
    Phiinv<-solve(Phi)
    
list(Phi = Phi, Phiinv = Phiinv, Pmat= Pmat, A1 = A1, A2 =A2, a.vec1 = a.vec1, a.vec2 = a.vec2, epsilon = epsilon, Zstar = Zstar )	
}



psirhorho_het<-function(rho, residuals, Hmat, Zmat, Ws, step1.c){
	n<-length(residuals)	

if(step1.c)	Zstar <- Zmat
else Zstar<- Zmat - rho * Ws %*% Zmat

 # print(dim(Zstar))
	epsilon<- residuals - rho * Ws %*% residuals 
    
	Qhz<-crossprod(Hmat,Zstar)/n
	Qhh<-crossprod(Hmat)/n
	Qhhi<-solve(Qhh)
	sec.part<- t(Qhz) %*% Qhhi %*% Qhz
	Pmat<-Qhhi %*% Qhz %*% solve(sec.part)
	Tmat<-Hmat %*% Pmat	

	trA2A2 <- sum(Ws^2) 
	A1 <- crossprod(Ws) 
	diag(A1)<-0
	A2 <- Ws
	
	
	A1pluA1 <- A1 + t(A1)	
	A2pluA2 <- A2 + t(A2)
	A1pluA1eps<- A1pluA1 %*% epsilon
	A2pluA2eps<- A2pluA2 %*% epsilon


	IrWp <- t(Diagonal(n) - rho * Ws)
	ZpIrWp<-crossprod(Zmat,IrWp)
	
	ZA1pluA1eps<- ZpIrWp %*% A1pluA1eps
	ZA2pluA2eps<- ZpIrWp %*% A2pluA2eps
	
		
	alpha1 <- -1*ZA1pluA1eps/n
	alpha2 <- -1*ZA2pluA2eps/n

if(step1.c){
	
	IrWpsi <- solve(Diagonal(n) - rho * t(Ws))
	a.vec1 <- IrWpsi%*% Tmat %*% alpha1 
	a.vec2 <- IrWpsi%*% Tmat %*% alpha2 
}

else{	
	a.vec1 <- Tmat %*% alpha1 
	a.vec2 <- Tmat %*% alpha2 
}

	
	gamma<-epsilon^2
	gammas<-Matrix(diag(n), sparse = TRUE)
     diag(gammas) <- gamma

	a1Ga1<-t(a.vec1) %*% gammas%*% a.vec1
	a2Ga2<-t(a.vec2) %*% gammas%*% a.vec2
	a1Ga2<-t(a.vec1) %*% gammas%*% a.vec2
	
   tr1<-sum(diag(A1pluA1%*%gammas%*%A1pluA1%*%gammas))/2
   tr2<-sum(diag(A2pluA2%*%gammas%*%A2pluA2%*%gammas))/2
   tr3<-sum(diag(A1pluA1%*%gammas%*%A2pluA2%*%gammas))/2
   
   phi11<- as.matrix((tr1 + a1Ga1)/n)
   phi22<- as.matrix((tr2 + a2Ga2)/n)
   phi12<- as.matrix((tr3 + a1Ga2)/n)
   phir1<- cbind(phi11, phi12)
   phir2<- cbind(phi12, phi22)
   
	Phi<-rbind(as.matrix(phir1), as.matrix(phir2))
	Phiinv<-solve(as.matrix(Phi))
	
list(Phi = Phi, Phiinv = Phiinv, Pmat= Pmat, A1 = A1, A2 =A2, a.vec1 = a.vec1, a.vec2 = a.vec2, epsilon = epsilon, Zstar = Zstar )		
	}




Omega_hom<-function(rho, Pmat, A1, A2, a.vec1, a.vec2, Hmat, bigG, Phirri, epsilon, Zstar){
	n<-length(epsilon)
	Jota<- bigG %*% c(1,2*rho)
	mu3<-sum(epsilon^3) / n
	sigma2<-crossprod(epsilon)/n
	Qhh <- crossprod(Hmat)/n
	a.vec <- cbind(as.numeric(a.vec1), as.numeric(a.vec2))
	a.vec.d<-cbind(diag(A1), diag(A2))
	
	Psidd <- as.numeric(sigma2)	* Qhh
	Psidr1<- as.numeric(sigma2)/n * crossprod(Hmat, a.vec)
	Psidr2<- mu3/n * crossprod(Hmat, a.vec.d)
	Psidr<- Psidr1 + Psidr2
	
	Omegadd<- as.matrix(t(Pmat) %*% Psidd %*% Pmat)
	Omegarr<- as.matrix(solve(t(Jota) %*% Phirri %*% Jota))
	Omegadr<- as.matrix(t(Pmat) %*% Psidr %*% Phirri %*% Jota %*% Omegarr)
	Omega<- cbind(rbind(Omegadd, t(Omegadr)), rbind(Omegadr, Omegarr))/n
	list(Omega = Omega)
}


Omega_het<-function(rho, Pmat, A1, A2, a.vec1, a.vec2, Hmat, bigG, Phirri, epsilon, Zstar, Ws, step1.c){
	n<-length(epsilon)
	gamma<-epsilon^2
	gammas<-Matrix(diag(n), sparse = TRUE)
    diag(gammas) <- gamma
     
	Jota<- bigG %*% c(1,2*rho)
	a.vec <- cbind(as.numeric(a.vec1), as.numeric(a.vec2))

if(step1.c){
	IrW <- Diagonal(n) - rho * Ws
	IrWp <-	t(IrW)
	Psidd <-  (t(Hmat) %*% IrW %*% gammas %*% IrWp %*% Hmat)/n
	Psidr<- (t(Hmat) %*% IrW %*% gammas %*% a.vec)/n
}

else{	
	Psidd <-  (t(Hmat) %*% gammas %*% Hmat)/n
	Psidr<- (t(Hmat) %*% gammas %*% a.vec)/n
	}
	Omegadd<- as.matrix(t(Pmat) %*% Psidd %*% Pmat)
	Omegarr<- as.matrix(solve(t(Jota) %*% Phirri %*% Jota))
	Omegadr<- as.matrix(t(Pmat) %*% Psidr %*% Phirri %*% Jota %*% Omegarr)
	Omega<- cbind(rbind(Omegadd, t(Omegadr)), rbind(Omegadr, Omegarr))/n
	list(Omega = Omega)
}







psirhorho_het_mod<-function(rho, residuals, Hmat, Zmat, Ws, step1.c){
	n<-length(residuals)	
if(step1.c)	Zstar <- Zmat
else Zstar<- Zmat - rho * Ws %*% Zmat

	epsilon<- residuals - rho * Ws %*% residuals 
    
	Qhz<-crossprod(Hmat,Zstar)/n
	Qhh<-crossprod(Hmat)/n
	Qhhi<-solve(Qhh)
	sec.part<- t(Qhz) %*% Qhhi %*% Qhz
	Pmat<-Qhhi %*% Qhz %*% solve(sec.part)
	Tmat<-Hmat %*% Pmat	

	trA2A2 <- sum(Ws^2) 
	A1 <- crossprod(Ws) 
	diag(A1)<-0
	A2 <- Ws
	
	
	A1pluA1 <- A1 + t(A1)	
	A2pluA2 <- A2 + t(A2)
	A1pluA1eps<- A1pluA1 %*% epsilon
	A2pluA2eps<- A2pluA2 %*% epsilon


	IrWp <- t(Diagonal(n) - rho * Ws)
	ZpIrWp<-crossprod(Zmat,IrWp)
	
	ZA1pluA1eps<- ZpIrWp %*% A1pluA1eps
	ZA2pluA2eps<- ZpIrWp %*% A2pluA2eps
	
		
	alpha1 <- -1*ZA1pluA1eps/n
	alpha2 <- -1*ZA2pluA2eps/n

if(step1.c){
	
	IrWpsi <- solve(Diagonal(n) - rho * t(Ws))
	a.vec1 <- IrWpsi%*% Tmat %*% alpha1 
	a.vec2 <- IrWpsi%*% Tmat %*% alpha2 
}

else{	
	a.vec1 <- Tmat %*% alpha1 
	a.vec2 <- Tmat %*% alpha2 
}

	
	gamma<-epsilon^2
	gammas<-Matrix(diag(n), sparse = TRUE)
     diag(gammas) <- gamma

	a1Ga1<-t(a.vec1) %*% gammas%*% a.vec1
	a2Ga2<-t(a.vec2) %*% gammas%*% a.vec2
	a1Ga2<-t(a.vec1) %*% gammas%*% a.vec2
	
   tr1<-sum(diag(A1pluA1%*%gammas%*%A1pluA1%*%gammas))/2
   tr2<-sum(diag(A2pluA2%*%gammas%*%A2pluA2%*%gammas))/2
   tr3<-sum(diag(A1pluA1%*%gammas%*%A2pluA2%*%gammas))/2
   
   phi11<- as.matrix((tr1 )/n)
   phi22<- as.matrix((tr2 )/n)
   phi12<- as.matrix((tr3 )/n)
   phir1<- cbind(phi11, phi12)
   phir2<- cbind(phi12, phi22)
   
	Phi<-rbind(as.matrix(phir1), as.matrix(phir2))
	Phiinv<-solve(as.matrix(Phi))
	
list(Phi = Phi, Phiinv = Phiinv, Pmat= Pmat, A1 = A1, A2 =A2, a.vec1 = a.vec1, a.vec2 = a.vec2, epsilon = epsilon, Zstar = Zstar )		
	}




psirhorho_hom_mod <- function(rho, residuals, Hmat, Zmat, Ws, d, v.vec){
	n <- length(residuals)  
	epsilon<- residuals - rho * Ws %*% residuals 
	
	mu3<- sum(epsilon^3) /n 
	mu4<- sum(epsilon^4) /n 	
	sigma2<-crossprod(epsilon)/n
	
	Zstar <- Zmat - rho * Ws %*% Zmat
	Qhz<-crossprod(Hmat,Zstar)/n
	Qhh<-crossprod(Hmat)/n
	Qhhi<-solve(Qhh)
	sec.part<- t(Qhz) %*% Qhhi %*% Qhz
	Pmat<-Qhhi %*% Qhz %*% solve(sec.part)
	Tmat<-Hmat %*% Pmat	
	
    trA2A2 <- sum(Ws^2) 
    trA2A2I<- Matrix(diag(n), sparse=TRUE)* (trA2A2/n)
	A1 <- v.vec * (crossprod(Ws) - trA2A2I)
	A2 <- Ws
	
	
	A1pluA1 <- A1 + t(A1)	
	A2pluA2 <- A2 + t(A2)
	A1pluA1eps<- A1pluA1 %*% epsilon
	A2pluA2eps<- A2pluA2 %*% epsilon

	ZA1pluA1eps<- crossprod(Zstar, A1pluA1eps)
	ZA2pluA2eps<- crossprod(Zstar, A2pluA2eps)
		
	alpha1 <- -1*ZA1pluA1eps/n
	alpha2 <- -1*ZA2pluA2eps/n
	
	a.vec1 <- Tmat %*% alpha1 
	a.vec2 <- Tmat %*% alpha2 
	

    tr22<-sum(diag(A2pluA2 %*% A2pluA2))/(2*n)
#   tr22<-(2*trA2A2 + 2* t(as.vector(t(Ws))) %*% as.vector(Ws)) /(2*n)
    tr11<-sum(diag(A1pluA1 %*% A1pluA1))/(2*n)
	tr12<-sum(diag(A1pluA1 %*% A2pluA2))/(2*n)

	a1.a1 <-crossprod(a.vec1) /n 
	a2.a2 <-crossprod(a.vec2) /n 
	a1.a2 <-crossprod(a.vec1, a.vec2) /n 

    vecA1<- diag(A1)
    vecA2<- diag(A2)

    vecA1vecA1<- crossprod(vecA1)/n
    vecA2vecA2<- crossprod(vecA2)/n
    vecA1vecA2<- crossprod(vecA1,vecA2)/n 

     a1.vecA1<- crossprod(a.vec1, vecA1)
     a1.vecA2<- crossprod(a.vec1, vecA2)
     a2.vecA2<- crossprod(a.vec2, vecA2)
     a2.vecA1<- crossprod(a.vec2, vecA1)    
    
    la.term11 <- (a1.vecA1 + a1.vecA1)/n
    la.term12 <- (a1.vecA2 + a2.vecA1)/n
    la.term22 <- (a2.vecA2 + a2.vecA2)/n
    
    
     effe<- mu4-3*sigma2^2
    
    # phi11<- as.numeric(sigma2^2 * tr11 )
    # phi12<- as.numeric(sigma2^2 * tr12 )
    # phi21<- as.numeric(sigma2^2 * tr12 )
    # phi22<- as.numeric(sigma2^2 * tr22 )

    phi11<- as.numeric(sigma2^2 * tr11 + effe *    vecA1vecA1 )
    phi12<- as.numeric(sigma2^2 * tr12 + effe *    vecA1vecA2 )
    phi21<- as.numeric(sigma2^2 * tr12 + effe *    vecA1vecA2 )
    phi22<- as.numeric(sigma2^2 * tr22 + effe *    vecA2vecA2 )
  

     Phi<-cbind(rbind(phi11,phi21),rbind(phi12,phi22))
    Phiinv<-solve(Phi)
    
list(Phi = Phi, Phiinv = Phiinv, Pmat= Pmat, A1 = A1, A2 =A2, a.vec1 = a.vec1, a.vec2 = a.vec2, epsilon = epsilon, Zstar = Zstar )	
}



Omega_hom_mod<-function(rho, Pmat, A1, A2, a.vec1, a.vec2, Hmat, bigG, Phirri, epsilon, Zstar){

	n<-length(epsilon)
	Jota<- bigG %*% c(1,2*rho)
	mu3<-sum(epsilon^3) / n
	sigma2<-crossprod(epsilon)/n
	Qhh <- crossprod(Hmat)/n
	 a.vec <- cbind(as.numeric(a.vec1), as.numeric(a.vec2))
	 a.vec.d<-cbind(diag(A1), diag(A2))
	
	Psidd <- as.numeric(sigma2)	* Qhh
	 Psidr1<- as.numeric(sigma2)/n * crossprod(Hmat, a.vec)
	Psidr2<- mu3/n * crossprod(Hmat, a.vec.d)
	Psidr<-  Psidr2
	
	Omegadd<- as.matrix(t(Pmat) %*% Psidd %*% Pmat)
	Omegarr<- as.matrix(solve(t(Jota) %*% Phirri %*% Jota))
	Omegadr<- as.matrix(t(Pmat) %*% Psidr %*% Phirri %*% Jota %*% Omegarr)
	Omega<- cbind(rbind(Omegadd, t(Omegadr)), rbind(Omegadr, Omegarr))/n
	list(Omega = Omega)
}


Omega_het_mod<-function(rho, Pmat, A1, A2, a.vec1, a.vec2, Hmat, bigG, Phirri, epsilon, Zstar, Ws, step1.c){
	n<-length(epsilon)
	gamma<-epsilon^2
	gammas<-Matrix(diag(n), sparse = TRUE)
    diag(gammas) <- gamma
     
	Jota<- bigG %*% c(1,2*rho)
	 a.vec <- cbind(as.numeric(a.vec1), as.numeric(a.vec2))

if(step1.c){
	IrW <- Diagonal(n) - rho * Ws
	IrWp <-	t(IrW)
	Psidd <-  (t(Hmat) %*% IrW %*% gammas %*% IrWp %*% Hmat)/n
	Psidr<- (t(Hmat) %*% IrW %*% gammas %*% a.vec)/n
}

else{	
	Psidd <-  (t(Hmat) %*% gammas %*% Hmat)/n
	 Psidr<- (t(Hmat) %*% gammas %*% a.vec)/n
	}
	Omegadd<- as.matrix(t(Pmat) %*% Psidd %*% Pmat)
	Omegarr<- as.matrix(solve(t(Jota) %*% Phirri %*% Jota))
	 Omegadr<- as.matrix(t(Pmat) %*% Psidr %*% Phirri %*% Jota %*% Omegarr)
	 # print(dim(Omegadr))
	Omega<- cbind(rbind(Omegadd, matrix(0,1,2)), rbind(matrix(0,2,1), Omegarr))/n
	list(Omega = Omega)
}
