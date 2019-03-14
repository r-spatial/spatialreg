# Copyright 2019 by Roger Bivand 
#

.spregOptions <- new.env(FALSE, globalenv())
#assign("spChkID", FALSE, envir = .spregOptions)
assign("zeroPolicy", spdep::get.ZeroPolicyOption(), envir = .spregOptions)
assign("verbose", spdep::get.VerboseOption(), envir = .spregOptions)
#assign("mc", ifelse(.Platform$OS.type == "windows", FALSE, TRUE),
# envir = .spregOptions)
#assign("cores", NULL, envir = .spregOptions)
#assign("cluster", NULL, envir = .spregOptions)
#assign("rlecuyerSeed", rep(12345, 6), envir = .spregOptions)
#assign("listw_is_CsparseMatrix", FALSE, envir = .spregOptions)

#setOldClass(c("listw"))


