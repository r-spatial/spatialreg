library(spatialreg)
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb)
COL.OLD$fEW <- factor(COL.OLD$EW)
COL.OLD$fDISCBD <- ordered(cut(COL.OLD$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(CRIME ~ INC + HOVAL + fDISCBD*fEW)
expect_warning(COL.SLX0 <- lmSLX(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.SLX1 <- lmSLX(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.SLX2 <- lmSLX(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fEW))
expect_silent(COL.SLX3 <- lmSLX(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))
expect_warning(COL.err0 <- errorsarlm(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.err1 <- errorsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.err2 <- errorsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD))
expect_silent(COL.err3 <- errorsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))
expect_warning(COL.lag0 <- lagsarlm(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.lag1 <- lagsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.lag2 <- lagsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD))
expect_silent(COL.lag3 <- lagsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))
expect_warning(COL.sac0 <- sacsarlm(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.sac1 <- sacsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.sac2 <- sacsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD))
expect_silent(COL.sac3 <- sacsarlm(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))
expect_warning(COL.lag0 <- spBreg_lag(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.lag1 <- spBreg_lag(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.lag2 <- spBreg_lag(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD))
expect_silent(COL.lag3 <- spBreg_lag(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))
expect_warning(COL.err0 <- spBreg_err(f, data=COL.OLD, lw, Durbin=TRUE))
expect_warning(COL.err1 <- spBreg_err(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD*fEW))
expect_warning(COL.err2 <- spBreg_err(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL + fDISCBD))
expect_silent(COL.err3 <- spBreg_err(f, data=COL.OLD, lw, Durbin=~ INC + HOVAL))

