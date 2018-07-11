stslshac<-function(formula, data=list(),listw,na.action=na.fail,zero.policy=NULL,HAC=TRUE, distance=NULL,type=c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth="variable",W2X=TRUE){

##functions that need to be sourced
	#source("twostagels.R")
	#source("kernelsfun.R")	
#extract model objects	
	
	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))

	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action,method="model.frame")
	na.act<-attr(mf,'na.action')


#preferences on missings values
if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    

#check that W ius an object of class listw
if(!inherits(listw,"listw")) 
	stop("The weights matrix is not a listw object")
	
##check that an exiting kernel is specified
if(HAC){
if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
#check that the distance measure is specificed
if(is.null(distance) ) stop("No distance measure specified")

#check that dist is an object of class sphet distance
if(!inherits(distance,"distance")) 
	stop("The distance measure is not a distance object")

}
#generates x and y 
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)

cl<-match.call()
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
if(W2X)      wwx<- lag.listw(listw, wx, zero.policy = zero.policy)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WX[, (i - (K - 1))] <- wx
				 WWX[, (i - (K - 1))] <- wwx
				         }
    }

#instr<-cbind(WX[,-c(1:8)],WWX[,-c(1:8)]) 
instr<-cbind(WX,WWX) 
#print(cbind(x,instr)[1,])

##spatial two stage least square of the initial model
#print(type)
results<-tsls(y=y,yend=wy, X=x, Zinst = instr, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance, zero.policy=zero.policy)
model.data<-data.frame(cbind(y,x[,-1]))

results$call<-cl
results$model<-model.data
results$type<-type
results$bandwidth<-bandwidth
results$method<-"s2slshac"
results$HAC<-HAC
results$zero.policy<-zero.policy
class(results)<-c("sphet", "stsls")

return(results)
}
