Omega<-function(n, gamma, H, param, G, FIinv, a, FI, P, gammas,Ws) {
	HpW<-t(H)%*%Ws
	rHpW<- -param*HpW
	HprHpW<-t(H)+rHpW
	HprHpWG<-HprHpW%*%gammas
	WpH<-t(Ws)%*%H
	rWpH<- -param*WpH
	HrWpH<-H+rWpH
	FIdd<-(HprHpWG%*%HrWpH)/n
#	print(a)
	FIdr<-as.matrix((HprHpWG%*%a)/n)
	FIzero1<-cbind(as.matrix(FIdd), FIdr)
	FIzero2<-cbind(t(FIdr), as.matrix(FI))
	Fizero<-rbind(as.matrix(FIzero1),as.matrix(FIzero2))
	w<-rbind(1,2*param)
	J<-G%*%w
	JFJ1<- t(J) %*% FIinv %*% J
     JFJ<- JFJ1^(-1)
     JFIinv<-crossprod(J,FIinv)
     OM1<- JFJ %*% JFIinv 
     OM2<- FIinv %*% J %*% JFJ
     P<-as.matrix(P)
      OMEGA1r1<-cbind(t(P),  matrix(0, nrow = nrow(t(P)), ncol= ncol(OM1)) )
     OMEGA1r2<-cbind(matrix(0, nrow= nrow(OM1), ncol = ncol(t(P)) ), OM1)

     OMEGA1<- rbind(OMEGA1r1, OMEGA1r2)
     
     OMEGA2r1<-cbind(P, matrix(0, nrow= nrow(P), ncol= ncol(OM2)))
	 OMEGA2r2<-cbind(matrix(0, nrow= nrow(OM2), ncol= ncol(P)), OM2)
     OMEGA2<-rbind(OMEGA2r1, OMEGA2r2)
    OMEGAA<-OMEGA1 %*% Fizero %*% OMEGA2
     OMEGA<-OMEGAA/n
	 list(VarCov=OMEGA)
	}
	
	Omegabis<-function(gammas,Hs,param,G,FIinv,as,n,FI,P){
	FIdd<-as.matrix((t(Hs)%*%gammas%*%Hs)/n)
	FIdr<-as.matrix((t(Hs)%*%gammas%*%as)/n)
     Fizeror1<- cbind(FIdd, FIdr)
     Fizeror2<- cbind(t(FIdr), FI)
     FIzeros<- rbind(as.matrix(Fizeror1), as.matrix(Fizeror2))
     w<-rbind(1, 2*param)
     J<-G %*% w
     JFJ1<- t(J) %*% FIinv %*% J
     JFJ<- solve(JFJ1)
     JFIinv<-crossprod(J,FIinv)
     OM1<- as.matrix(JFJ %*% JFIinv )
     OM2<- as.matrix(FIinv %*% J %*% JFJ)
     # print(dim(P))
     # print(dim(OM1))
      # print(matrix(0, nrow= nrow(t(P)), ncol= ncol(OM1) ))
      # print(t(P))
      P<-as.matrix(P)

     OMEGA1r1<-cbind(t(P),  matrix(0, nrow = nrow(t(P)), ncol= ncol(OM1)) )
     OMEGA1r2<-cbind(matrix(0, nrow= nrow(OM1), ncol = ncol(t(P)) ), OM1)
     # print(OMEGA1r1)
     # print(dim(OMEGA1r2))

     OMEGA1 <- rbind(OMEGA1r1, OMEGA1r2)
# print(OMEGA1)
     OMEGA2r1<-cbind(P, matrix(0, nrow= nrow(P), ncol= ncol(OM2)))
	 OMEGA2r2<-cbind(matrix(0, nrow= nrow(OM2), ncol= ncol(P)), OM2)
	 OMEGA2<-rbind(OMEGA2r1, OMEGA2r2)
	OMEGAA<-OMEGA1 %*% FIzeros %*% OMEGA2
	OMEGA<-OMEGAA/n
	 list(VarCov=OMEGA)
		}
		
		
