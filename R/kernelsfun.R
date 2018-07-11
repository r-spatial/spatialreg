triangular <- function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z <- q/ban
		k <- rep(0, length(z))
        k[which(z<1)  ] <- 1-z[which(z<1)  ]
		k
		}

rectangular<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
	    z <-q/ban
		k <- rep(0, length(z))
        k[which(z<1)  ] <- 1
		k
		k	
}



epan<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
				z<-q/ban			
		k <- rep(0, length(z))
        k[which(z<1)  ] <- 1 - z[which(z<1)]^2
		k
				}
				
bisq<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k <- rep(0, length(z))		
		k[which(z<1)] <- (1-z[which(z<1)]^2)^2
		}


parzen<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z <- q/ban
		k <- rep(0, length(z))		
		tmp1<-which(z<=0.5)
		tmp2 <-which(z > 0.5 & z<1)
k[tmp1]<-1-6*z[tmp1]^2+6*abs(z[tmp1])^3
k[tmp2]<-2*(1-abs(z[tmp2]))^3
		k
		}

th <- function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k <- rep(0, length(z))		
k[which(z<1)] <- (1+cos(pi*z[which(z<1)]))/2
		k
		}


qs<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k <- rep(0, length(z))		
k[which(z<1)] <-(25/(12*pi^2*z[which(z<1)]^2)) *((sin((6*pi*z[which(z<1)])/5)/((6*pi*z[which(z<1)])/5))-cos((6*pi*z[which(z<1)])/5))
		k
		}
