summary.sphet<-function(object, width=getOption("width"),digits=getOption("digits"),obsinfo=FALSE,...){

	if(inherits(object,"distance")){
		weights<-	attributes(object)$GeoDa$dist
		neighbors<- object
		n<- length(weights)
		maxdist<-unlist(lapply(weights,max))
		nneigh<-unlist(lapply(weights,length))
		object$weights<-weights
		object$n<-n
		object$maxdist <- maxdist
		object$nneigh <- nneigh
		class(object)<- c('summary.sphet','sphet','distance')
		}
	else{
	coeff<-object$coefficients
	vcmat<-object$var
	se<-sqrt(diag(vcmat))
   t<-coeff/se
   pval<-pnorm(abs(t),lower.tail=FALSE) * 2
   CoefTable <- cbind(coeff,se,t,pval)
if(object$method=="s2slshac")  { 
	 if (object$HAC) colnames(CoefTable) <- c("Estimate","SHAC St.Er.","t-value","Pr(>|t|)")
	 else colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	}
else colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
   object$CoefTable <- CoefTable
	class(object)<- c('summary.sphet','sphet')
	}
	object
	}
	
	print.summary.sphet <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),obsinfo=FALSE,...){
if(inherits(x,"distance")){
	cat("\n Number of observations:\n")
	cat(" n:", x$n, "\n")
if (obsinfo){
	cat("\n Maximum distance for each observation:\n") 
	cat(x$maxdist,  fill=TRUE)
	}
	cat("\n Distance summary:\n")
	print(summary(x$maxdist), width=width)
if (obsinfo){	
	cat("\n Number of non-zero elements for each observation:\n")
	cat( x$nneigh, fill=TRUE)
	}
	cat("\n Neighbors summary:\n")	
	print(summary(x$nneigh), width=width)
	}
			
			else{
				
if (x$method=="s2slshac") {
if (x$HAC) cat("\n Stsls with Spatial HAC standard errors\n")
else cat("\n Stsls \n")
}
else cat("\n Generalized stsls\n")		  				
		  cat("\nCall:\n")
  print(x$call)

  cat("\nResiduals:\n")
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(sumres(x))


		  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable,digits=digits)

if(x$method=="gs2slshac") {
	
	if(!is.null(x$W)){
		  cat("\nWald test that rho and lambda are both zero:\n")
cat(" Statistics:", x$W$stat,"p-val:", x$W$pval, "\n")
}
}
}
  cat("\n")  
  invisible(x)
}

sumres <- function(x){
  sr <- summary(residuals(x))
  srm <- mean(residuals(x))
  if (abs(srm) < 1e-10){
    sr <- sr[c(1:3,5:6)]
  }
  sr
}


print.sphet <- function(x, digits = max(3, getOption("digits") - 3),...) 
{
	if(inherits(x,"distance")){ 
		tmp<-summary.sphet(x)
		print.summary.sphet(tmp)		
		}

else{	
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
}
    invisible(x)
}
