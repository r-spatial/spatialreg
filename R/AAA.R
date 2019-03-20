# Copyright 2019 by Roger Bivand 
#

.spregOptions <- new.env(FALSE, globalenv())
#assign("spChkID", FALSE, envir = .spregOptions)
assign("zeroPolicy", spdep::get.ZeroPolicyOption(), envir = .spregOptions)
assign("verbose", spdep::get.VerboseOption(), envir = .spregOptions)
assign("mc", spdep::get.mcOption(), envir = .spregOptions)
assign("cores", spdep::get.coresOption(), envir = .spregOptions)
assign("cluster", spdep::get.ClusterOption(), envir = .spregOptions)
#assign("rlecuyerSeed", rep(12345, 6), envir = .spregOptions)
#assign("listw_is_CsparseMatrix", FALSE, envir = .spregOptions)

setOldClass(c("listw"))

Hausman.test <- function(object, ...)
    UseMethod("Hausman.test", object)

