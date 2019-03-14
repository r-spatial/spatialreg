# Copyright 2003-2015 by Roger Bivand 

set.listw_is_CsparseMatrix_Option <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("listw_is_CsparseMatrix", envir = .spregOptions)
	assign("listw_is_CsparseMatrix", check, envir = .spregOptions)
	res
}

get.listw_is_CsparseMatrix_Option <- function() {
	get("listw_is_CsparseMatrix", envir = .spregOptions)
}

set.VerboseOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("verbose", envir = .spregOptions)
	assign("verbose", check, envir = .spregOptions)
	res
}

get.VerboseOption <- function() {
	get("verbose", envir = .spregOptions)
}

set.ZeroPolicyOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("zeroPolicy", envir = .spregOptions)
	assign("zeroPolicy", check, envir = .spregOptions)
	res
}

get.ZeroPolicyOption <- function() {
	get("zeroPolicy", envir = .spregOptions)
}

set.ClusterOption <- function(cl) {
	if (!is.null(cl)) {
            if (!inherits(cl, "cluster")) stop ("cluster required")
        }
	assign("cluster", cl, envir = .spregOptions)
        invisible(NULL)
}

get.ClusterOption  <- function() {
	get("cluster", envir = .spregOptions)
}

set.mcOption <- function(value) {
        stopifnot(is.logical(value))
        stopifnot(length(value) == 1)
	res <- get("mc", envir = .spregOptions)
        if (.Platform$OS.type == "windows") {
            if (value) warning("multicore not available on Windows")
        } else {
	    assign("mc", value, envir = .spregOptions)
        }
	res
}

get.mcOption  <- function() {
	get("mc", envir = .spregOptions)
}

set.coresOption <- function(value) {
	res <- get("cores", envir = .spregOptions)
        if (is.null(value)) {
            assign("cores", value, envir = .spregOptions)
        } else {
            stopifnot(is.integer(value))
            stopifnot(length(value) == 1)
            stopifnot(!is.na(value))
	    assign("cores", value, envir = .spregOptions)
        }
	res
}

get.coresOption  <- function() {
	get("cores", envir = .spregOptions)
}


