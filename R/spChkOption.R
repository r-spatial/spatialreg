# Copyright 2003-2015 by Roger Bivand 

set.listw_is_CsparseMatrix_Option <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("listw_is_CsparseMatrix", envir = .spatialregOptions)
	assign("listw_is_CsparseMatrix", check, envir = .spatialregOptions)
	res
}

get.listw_is_CsparseMatrix_Option <- function() {
	get("listw_is_CsparseMatrix", envir = .spatialregOptions)
}

set.VerboseOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("verbose", envir = .spatialregOptions)
	assign("verbose", check, envir = .spatialregOptions)
	res
}

get.VerboseOption <- function() {
	get("verbose", envir = .spatialregOptions)
}

set.ZeroPolicyOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("zeroPolicy", envir = .spatialregOptions)
	assign("zeroPolicy", check, envir = .spatialregOptions)
	res
}

get.ZeroPolicyOption <- function() {
	get("zeroPolicy", envir = .spatialregOptions)
}

set.ClusterOption <- function(cl) {
	if (!is.null(cl)) {
            if (!inherits(cl, "cluster")) stop ("cluster required")
        }
	assign("cluster", cl, envir = .spatialregOptions)
        invisible(NULL)
}

get.ClusterOption  <- function() {
	get("cluster", envir = .spatialregOptions)
}

set.mcOption <- function(value) {
        stopifnot(is.logical(value))
        stopifnot(length(value) == 1)
	res <- get("mc", envir = .spatialregOptions)
        if (.Platform$OS.type == "windows") {
            if (value) warning("multicore not available on Windows")
        } else {
	    assign("mc", value, envir = .spatialregOptions)
        }
	res
}

get.mcOption  <- function() {
	get("mc", envir = .spatialregOptions)
}

set.coresOption <- function(value) {
	res <- get("cores", envir = .spatialregOptions)
        if (is.null(value)) {
            assign("cores", value, envir = .spatialregOptions)
        } else {
            stopifnot(is.integer(value))
            stopifnot(length(value) == 1)
            stopifnot(!is.na(value))
	    assign("cores", value, envir = .spatialregOptions)
        }
	res
}

get.coresOption  <- function() {
	get("cores", envir = .spatialregOptions)
}


