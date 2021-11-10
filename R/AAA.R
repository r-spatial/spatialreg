# Copyright 2019 by Roger Bivand 
#

.spatialregOptions <- new.env(FALSE, globalenv())
#assign("spChkID", FALSE, envir = .spatialregOptions)
assign("zeroPolicy", spdep::get.ZeroPolicyOption(), envir = .spatialregOptions)
assign("verbose", spdep::get.VerboseOption(), envir = .spatialregOptions)
assign("mc", spdep::get.mcOption(), envir = .spatialregOptions)
assign("cores", spdep::get.coresOption(), envir = .spatialregOptions)
assign("cluster", spdep::get.ClusterOption(), envir = .spatialregOptions)
assign("rlecuyerSeed", rep(12345, 6), envir = .spatialregOptions)
#assign("listw_is_CsparseMatrix", FALSE, envir = .spatialregOptions)

setOldClass(c("listw"))

Hausman.test <- function(object, ...)
    UseMethod("Hausman.test")

