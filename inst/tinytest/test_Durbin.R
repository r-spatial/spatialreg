library(spatialreg)
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb)
f <- formula(CRIME ~ INC + HOVAL)
COL.mix.eig <- lagsarlm(f, data=COL.OLD, lw, Durbin=TRUE)
p1 <- predict(COL.mix.eig)
p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2, p3, check.attributes=FALSE)))
COL.mix.eig.f <- lagsarlm(f, data=COL.OLD, lw, Durbin=update(f, ~ . - HOVAL))
p2f <- predict(COL.mix.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3f <- predict(COL.mix.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2f, p3f, check.attributes=FALSE)))
COL.SDerr.eig <- errorsarlm(f, data=COL.OLD, lw, Durbin=TRUE)
p1 <- predict(COL.SDerr.eig)
p2 <- predict(COL.SDerr.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3 <- predict(COL.SDerr.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2, p3, check.attributes=FALSE)))
COL.SDerr.eig.f <- errorsarlm(f, data=COL.OLD, lw, Durbin=update(f, ~ . - HOVAL))
p2f <- predict(COL.SDerr.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3f <- predict(COL.SDerr.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2f, p3f, check.attributes=FALSE)))

