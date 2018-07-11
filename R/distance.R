distance <- function(coord,region.id=NULL,output=TRUE,type=c("NN","distance","inverse"),measure=c("euclidean","gcircle","chebyshev","braycur","canberra"),nn=6, cutoff=FALSE, miles=TRUE,R=NULL, shape.name=NULL,region.id.name=NULL,firstline=FALSE,file.name=NULL) {
	
	if ( !type %in% c("NN","distance","inverse") ) stop("unknown type")
	if ( !measure[1] %in% c("euclidean","gcircle","chebyshev","braycur","canberra") ) stop("unknown measure")
	if(!inherits(coord,c("data.frame","matrix"))){
		coord <-  as.matrix(coord)
		}
	if (is.null(region.id) && ncol(coord)!=3) {
		id<-seq(1,nrow(coord))
		warning("region.id variable not specified")
		}
	if (is.null(region.id) && ncol(coord)==3) id<-coord[,1]
	if (!is.null(region.id) && ncol(coord) == 3){ 
		check<- all.equal(region.id,coord[,1])
		if ( check!=TRUE ) stop ("region.id and coord[,1] are different")
		id <- region.id
		}
	if (!is.null(region.id) && ncol(coord) != 3) id<-region.id 
	
	if (length(unique(id)) != length(id)) 
            stop("non-unique region.id given")

		
    k<- ifelse(ncol(coord)==3,2,1)
	 x<- coord[,k]
	 y<- coord[,k+1]
    h <- length(x)
    hh2 <- h*(h-1)/2
    hh <- h*(h-1)
    dd <- matrix(,nrow=hh2,ncol=3) 
	 vec<-id
	 vec1<-x
	 vec2<-y
		coord<-cbind(vec,vec1,vec2)
Weights<-cbind(0,0,0)
		Sq2<-c(1,seq(h-1,2))
		#print(Sq2)
		SSq2<-cumsum(Sq2)
		#print(SSq2)
						for (i in vec[-length(vec)]) {
									ref<-	vec[i]
									w<-which(vec!=i & vec > i)

									coordi<-cbind(vec1[i],vec2[i])
									coordj<-cbind(vec1[w],vec2[w])
####

		dist.fun<-switch(match.arg(measure), euclidean={
			dist.euclidean
			}, gcircle={
				dist.gcircle
				}, chebyshev= {
					dist.chebyshev
					}, braycur={
						dist.braycur
						}, canberra = {
							dist.canberra
								})

weights<-dist.fun(coordi,coordj,miles=TRUE,R=NULL)
if(measure[1]=="gcircle") weights<-t(weights)            
									
#####									
									
if(type=="inverse") tmp <- rbind(cbind(rep(vec[i],length(w)),vec[w],1/weights),cbind(vec[w],rep(vec[i],length(w)),1/weights))
else tmp <- rbind(cbind(rep(vec[i],length(w)),vec[w],weights),cbind(vec[w],rep(vec[i],length(w)),weights))

#print(tmp)
						if (!is.numeric(weights)) stop
									Weights<-rbind(Weights,tmp)
			}
			Weights<-Weights[-1,]

	if (type=="NN"){
		or<-order(Weights[,1],Weights[,3])
Weights<- Weights[or,]
final<-matrix(,h*nn,3)
for (i in 1:h){
	tmp<-Weights[Weights[,1]==i,][1:nn,]
	final[(((i*nn)-(nn-1)):(i*nn)),]<-as.matrix(tmp)
	}
Weights<-final
		}	
		
if (type=="distance" && cutoff){
			tmp<-Weights[,3]
		quant<-quantile(tmp)[cutoff+1]
tokeep<-which(tmp <= quant)
#print(tokeep)
Weights<-Weights[tokeep,]
}

if (type=="inverse" && cutoff){
			tmp<-Weights[,3]
		quant<-quantile(tmp)[cutoff+1]
tokeep<-which(tmp <= quant)
Weights<-Weights[tokeep,]
}

	    worder <- order(Weights[,1],Weights[,2])
   		 Weights<-Weights[worder,]
    colnames(Weights)<-c("from","to","distance")
    class(Weights)<- c("matrix","distance.matrix")
   
    if (output) {
   		 firstl<-c("0",h,shape.name,region.id.name)
   		 #print(firstline)
if(firstline) write(firstl,file.name,ncolumns=4)
    	write.table(Weights,file.name,row.names=FALSE,quote=FALSE,append=TRUE,col.names=FALSE)
    	}
Weights
    }
