# The assumption is that listw1 and listw2 are different, otherwise there would be no reason to use the spatial J-test

kpjtest <- function(H0model, H1model, data = list(), listw0 = NULL, listw1 = NULL, endogH0 = NULL, endogH1 = NULL, instrumentsH0 = NULL, instrumentsH1 = NULL, lag.instr = FALSE, model = c("lag", "sarar"), het = FALSE, HAC = F, distance = NULL, type = c("Epanechnikov", "Triangular", "Bisquare", "Parzen", "QS", "TH", "Rectangular"), bandwidth = "variable",  na.action = na.fail){
	 cl <- match.call()
	if(is.null(listw0)) 
		stop("listw for the null model was not specified")
	if(is.null(listw1))
		stop("listw for the alternative model was not specified")
	
	
	 if (!inherits(listw0, c("listw", "Matrix", "matrix"))) 
            stop("listw format unknown")
        if (inherits(listw0, "listw"))  Ws0 <- listw2dgCMatrix(listw0)
        if (inherits(listw0, "matrix")) Ws0 <- Matrix(listw0)
        if (inherits(listw0, "Matrix")) Ws0 <- listw0

	 if (!inherits(listw1, c("listw", "Matrix", "matrix"))) 
            stop("listw format unknown")
        if (inherits(listw1, "listw"))  Ws1 <- listw2dgCMatrix(listw1)
        if (inherits(listw1, "matrix")) Ws1 <- Matrix(listw1)
        if (inherits(listw1, "Matrix")) Ws1 <- listw1
	
	if((dim(Ws0)[1] != dim(Ws1)[1]) &&  (dim(Ws0)[2] != dim(Ws1)[2])) 
			stop("listw0 and listw1 must have the same dimension")
			
	if( all(Ws0 == Ws1)) 
			stop("listw0 and listw1 cannot be the same")		
	
	
if(model == "lag") res <- jtestlag(H0model, H1model, data = data, listw0 = Ws0, listw1 = Ws1, endogH0 = endogH0, endogH1 = endogH1, instrumentsH0 = instrumentsH0, instrumentsH1 = instrumentsH1, lag.instr = lag.instr, model = model,  HAC = HAC, distance = distance, type = type, bandwidth = bandwidth,  na.action = na.fail, het = het, cl = cl)


if(model == "sarar") stop ("Method for sarar not yet implemented")
	
	
	return(res)
}



jtestlag <- function(H0model, H1model, data = list(), listw0, listw1, endogH0, endogH1, instrumentsH0, instrumentsH1, lag.instr, model,  HAC, distance, type, bandwidth,  na.action, het, cl){
	
Alt <- spreg(H1model, data = data, listw = listw1, endog = endogH1, instruments = instrumentsH1, lag.instr = lag.instr, model = "lag", het = F)	

mt <- terms(H1model, data = data)
mf <- lm(H1model, data, na.action = na.action, method = "model.frame")
y1 <- c(model.extract(mf, "response"))
x1 <- model.matrix(mt, mf)
wy1 <- as.matrix(listw1 %*% y1)
reg <- as.matrix(cbind(x1, wy1))
yp <- reg %*% coefficients(Alt)
data$yp <- yp
# H0fm <- as.formula(paste(names(mf)[1], paste(names(mf)[-1], collapse = " + "), sep = " ~" ))

H0 <- augmented(H0model, data = data, listw0 = listw0, listw1 = listw1, yp = yp, x1 = x1, endogH0 = endogH0, endogH1 = endogH1, instrumentsH0 = instrumentsH0, instrumentsH1 = instrumentsH1, lag.instr = lag.instr, model = model,  HAC = HAC, distance = distance, type = type, bandwidth = bandwidth,  na.action = na.fail, het = het, cl = cl)	
# print(H0)
return(H0)
	
}

# # jtestsarar <- function(){
	
	
# }



augmented <- function(H0model, data, listw0, listw1, yp, x1, endogH0, endogH1, instrumentsH0, instrumentsH1, lag.instr, model,  HAC, distance, type, bandwidth,  na.action, het, cl){
	
	
	    if (!is.null(endogH0) && is.null(instrumentsH0)) 
        stop("No instruments specified for the endogenous variable in the null model")

	
				if(colnames(x1)[1] == "(Intercept)")	 x1 <- x1[,-1]
				else x1 <- x1

						mt <- terms(H0model, data = data)
						mf <- lm(H0model, data, na.action = na.action, method = "model.frame")

						y0 <- c(model.extract(mf, "response"))
						x0 <- model.matrix(mt, mf)
						w0y0 <- listw0 %*% y0
						
				if(colnames(x0)[1] == "(Intercept)")	 x0 <- x0[,-1]
				else x0 <- x0
						
						xaug <- cbind(x0 , x1[,-which(colnames(x1) %in% colnames(x0))])
						w0x0 <- listw0 %*% x0
						w02x0 <- listw0 %*% w0x0
						w1x1 <- listw1 %*% x1
						w12x1 <- listw1 %*% w1x1

	

		if(is.null(instrumentsH0) && is.null(instrumentsH1) ){
		
						Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1)	
		}


		if(!is.null(instrumentsH0) && !is.null(instrumentsH1) ){
			
			instrumentsH0 <- as.matrix(lm(instrumentsH0, data, na.action = na.action, method = "model.frame"))
			instrumentsH1 <- as.matrix(lm(instrumentsH1, data, na.action = na.action, method = "model.frame"))
			instaug <- cbind(instrumentsH0 , instrumentsH1[,-which(colnames(instrumentsH1) %in% colnames(instrumentsH0))])
			
			if(lag.instr){
				w0inH0 <- listw0 %*% instrumentsH0
				ww0inH0 <- listw0 %*% w0inH0
				w1inH1 <- listw1 %*% instrumentsH1
				ww1inH1 <- listw1 %*% w1inH1
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instaug, w0inH0, ww0inH0, w1inH1, ww1inH1)			
			}	
							Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instaug )			
		}


		if(!is.null(instrumentsH0) && is.null(instrumentsH1) ){
						
						instrumentsH0 <- as.matrix(lm(instrumentsH0, data, na.action = na.action, method = "model.frame"))
						
						if(lag.instr){
				w0inH0 <- listw0 %*% instrumentsH0
				ww0inH0 <- listw0 %*% w0inH0
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH0, w0inH0, ww0inH0)			
			}
		Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH0)				
		}
	
	
			if(is.null(instrumentsH0) && !is.null(instrumentsH1) ){
						
						instrumentsH1 <- as.matrix(lm(instrumentsH1, data, na.action = na.action, method = "model.frame"))
						
						if(lag.instr){
				w1inH1 <- listw1 %*% instrumentsH1
				ww1inH1 <- listw1 %*% w1inH1
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH1, w1inH1, ww1inH1)			
			}
		Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH1)				
					
		}


Haug <- as.matrix(Haug)

if(!is.null(endogH0)){
	endogH0 <- as.matrix(lm(endogH0, data, na.action = na.action, method = "model.frame"))
	Zaug <- cbind(1, x0, as.matrix(w0y0), endogH0, as.matrix(yp))
}
else  Zaug <- cbind(1, x0, as.matrix(w0y0), as.matrix(yp))
colnames(Zaug) <- c("(Intercept)", colnames(x0), "w0y", "yp")

results <- spatial.ivreg(y0, Zaug, Haug, het = het, HAC = HAC, distance = distance, type = type, bandwidth = bandwidth)
    
 model.data <- data.frame(cbind(as.matrix(y0), as.matrix(Zaug)))
        results$call <- cl
        results$model <- model.data
        results$type <- type
        results$bandwidth <- bandwidth
        results$method <- "s2slshac"
        results$HAC <- HAC
         class(results) <- c("sphet", "stsls_sphet")
# print(summary(results))

# df <- nrow(Zaug) - ncol(Zaug)
# HHaug <- crossprod(Haug)
# HHaugi <- solve(HHaug)
# fp <- Haug %*% HHaugi
# sp <- crossprod(Haug,Zaug)
# Zaugp <- fp %*% sp
# daug <- solve(crossprod(Zaugp), crossprod(Zaugp, y0))

# raug <- y0 - Zaug %*% daug

# s2 <- as.numeric(crossprod(raug)/df)
# vc <- s2 * solve(crossprod(Zaugp))
# print(ret)

# class(results) <- c("KP_jtest", "sphet", "stsls_sphet")
return(results)
}



# summary.KP_jtest <- function(x, digits= max(3, getOption("digits") - 2),...){
# cp <- coefficients(x)
# vc <- x$var
# test <- cp[length(cp)]^2/vc[length(cp),length(cp)]

# STAT <- qchisq(0.05, 1, lower.tail=FALSE)

# # ret <- list(test, STAT)	
	
# cat("\n Kelejian and Piras J-test for spatial models:\n")
# cat("\n H0 vs H1:\n")
# cat("\n ---------------------------------------------\n")
# cat("\n Stat:\n")
# print(test)
# cat("\n Chi-sq (df = 1):\n")
# print(STAT)

# }