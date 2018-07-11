circular<-function(nrow,ncol, ab){
    nrow <- as.integer(nrow)
    if (nrow < 1) 
        stop("nrow nonpositive")
    ncol <- as.integer(ncol)
    if (ncol < 1) 
        stop("nrow nonpositive")
 ab <- as.integer(ab)
	if(ab< 1) 
		stop("ab nonpositive")


N<-nrow*ncol	
vec<-seq(1,N)
fp<-vec[(length(vec)-ab+1):length(vec)]
lp<-vec[1:ab]
vec<-c(fp,vec,lp)

res<-vector("list",length = N)
rownames <- character(N)

for(i in (ab+1): (N+ab)){
	
	res[[i-ab]]<-sort(vec[(i-ab):(i+ab)][-(ab+1)])
	rownames[i-ab] <- paste(vi2mrc(i-ab, nrow, ncol), collapse = ":")	
		}


class(res)<-"nb"
class(res) <- "nb"
attr(res, "call") <- match.call()
attr(res, "region.id") <- rownames
attr(res, "cell") <- TRUE
res <- sym.attr.nb(res)
res
	
	}


