# Contributed by Tobias Rũttenauer,
# https://github.com/ruettenauer/mixed/blob/master/spatialreg_intercept-tests_lmslx_errorsarlm.R
# https://github.com/r-spatial/spatialreg/issues/13#issuecomment-1450272546

library(sf)
library(spatialreg)
nc <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
row.names(nc) <- as.character(nc$FIPSNO)
nc$ft.SID74 <- sqrt(1000)*(sqrt(nc$SID74/nc$BIR74) + sqrt((nc$SID74+1)/nc$BIR74))
nc$ft.NWBIR74 <- sqrt(1000)*(sqrt(nc$NWBIR74/nc$BIR74) + sqrt((nc$NWBIR74+1)/nc$BIR74))
gal_file <- system.file("weights/ncCR85.gal", package="spData")[1]
ncCR85 <- spdep::read.gal(gal_file, region.id=nc$FIPSNO)
nc_lw <- spdep::nb2listw(ncCR85, style="B")
nc_lw_w <- spdep::nb2listw(ncCR85, style="W") # row-standardized for no intercept

#### --------- Tests for SLX--------- ####

### SLX with unstandardized W and intercept in lagX
m2b <- spatialreg::lmSLX(ft.SID74 ~ ft.NWBIR74, 
                              weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
# summary(m2b)
# summary(spatialreg::impacts(m2b))

# Impact equals coef
expect_true(isTRUE(all.equal(spatialreg::impacts(m2b)[[1]]$indirect, 
          m2b$coefficients["lag.ft.NWBIR74"], 
          tolerance = 1e-5, check.attributes = FALSE)))



### SLX with unstandardized W, intercept in lagX, and multiple covariates
m2c <- spatialreg::lmSLX(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                              weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
# summary(m2c)
# summary(spatialreg::impacts(m2c))

# Impacts equal coefs
m2c.imps <- spatialreg::impacts(m2c)
expect_true(isTRUE(all.equal(c(m2c.imps[[1]]$indirect["ft.NWBIR74"],
            m2c.imps[[1]]$indirect["east"],
            m2c.imps[[1]]$indirect["AREA"]), 
          c(m2c$coefficients["lag.ft.NWBIR74"],
            m2c$coefficients["lag.east"],
            m2c$coefficients["lag.AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m2c, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m2c.imps[[1]]$total["east"],
            summary(m2c.imps)$semat["east", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### SLX with standardized W, without intercept in lagX, and multiple covariates
m2d <- spatialreg::lmSLX(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw_w, Durbin=TRUE)
# summary(m2d)
# summary(spatialreg::impacts(m2d))

# Impacts equal coefs
m2d.imps <- spatialreg::impacts(m2d)
expect_true(isTRUE(all.equal(c(m2d.imps[[1]]$indirect["ft.NWBIR74"],
            m2d.imps[[1]]$indirect["east"],
            m2d.imps[[1]]$indirect["AREA"]), 
          c(m2d$coefficients["lag.ft.NWBIR74"],
            m2d$coefficients["lag.east"],
            m2d$coefficients["lag.AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m2d, linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m2d, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m2d.imps[[1]]$total["ft.NWBIR74"],
            summary(m2d.imps)$semat["ft.NWBIR74", "Total"],
            m2d.imps[[1]]$total["east"],
            summary(m2d.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with standardized W
m2e <- spatialreg::lmSLX(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                              weights=BIR74, data=nc, listw=nc_lw_w, 
                              Durbin = ~ ft.NWBIR74 + east)
# summary(m2e)
# summary(spatialreg::impacts(m2e))

# Impacts equal coefs
m2e.imps <- spatialreg::impacts(m2e)
expect_true(isTRUE(all.equal(c(m2e.imps[[1]]$indirect["ft.NWBIR74"],
            m2e.imps[[1]]$indirect["east"]), 
          c(m2e$coefficients["lag.ft.NWBIR74"],
            m2e$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m2e, linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m2e, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m2e.imps[[1]]$total["ft.NWBIR74"],
            summary(m2e.imps)$semat["ft.NWBIR74", "Total"],
            m2e.imps[[1]]$total["east"],
            summary(m2e.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with unstandardized W and intercept
m2f <- spatialreg::lmSLX(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw, 
                         Durbin = ~ ft.NWBIR74 + east)
# summary(m2f)
# summary(spatialreg::impacts(m2f))

# Impacts equal coefs
m2f.imps <- spatialreg::impacts(m2f)
expect_true(isTRUE(all.equal(c(m2f.imps[[1]]$indirect["ft.NWBIR74"],
            m2f.imps[[1]]$indirect["east"]), 
          c(m2f$coefficients["lag.ft.NWBIR74"],
            m2f$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m2f, 
               linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m2f, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m2f.imps[[1]]$total["ft.NWBIR74"],
            summary(m2f.imps)$semat["ft.NWBIR74", "Total"],
            m2f.imps[[1]]$total["east"],
            summary(m2f.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))


#### --------- Tests intercept in SLX formula --------- ####


### SLX with unstandardized W and intercept in lagX
m3b <- spatialreg::lmSLX(ft.SID74 ~ -1 + ft.NWBIR74, 
                         weights=BIR74, data=nc, listw=nc_lw_w, Durbin=TRUE)
# summary(m3b)
# summary(spatialreg::impacts(m3b))

# Impact equals coef
m3b.imps <- spatialreg::impacts(m3b)
expect_true(isTRUE(all.equal(m3b.imps[[1]]$indirect, 
          m3b$coefficients["lag.ft.NWBIR74"], 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m3b, linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
expect_true(isTRUE(all.equal(c(m3b.imps[[1]]$total,
            summary(m3b.imps)$semat["ft.NWBIR74", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### SLX with unstandardized W, no intercept in X, and multiple covariates
m3c <- spatialreg::lmSLX(ft.SID74 ~ - 1 + ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
# summary(m3c)
# summary(spatialreg::impacts(m3c))

# Impacts equal coefs
m3c.imps <- spatialreg::impacts(m3c)
expect_true(isTRUE(all.equal(c(m3c.imps[[1]]$indirect["ft.NWBIR74"],
            m3c.imps[[1]]$indirect["east"],
            m3c.imps[[1]]$indirect["AREA"]), 
          c(m3c$coefficients["lag.ft.NWBIR74"],
            m3c$coefficients["lag.east"],
            m3c$coefficients["lag.AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m3c, 
               linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m3c, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m3c.imps[[1]]$total["ft.NWBIR74"],
            summary(m3c.imps)$semat["ft.NWBIR74", "Total"],
            m3c.imps[[1]]$total["east"],
            summary(m3c.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with standardized W
m3e <- spatialreg::lmSLX(ft.SID74 ~ - 1 + ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw_w, 
                         Durbin = ~ ft.NWBIR74 + east)
# summary(m3e)
# summary(spatialreg::impacts(m3e))

# Impacts equal coefs
m3e.imps <- spatialreg::impacts(m3e)
expect_true(isTRUE(all.equal(c(m3e.imps[[1]]$indirect["ft.NWBIR74"],
            m3e.imps[[1]]$indirect["east"]), 
          c(m3e$coefficients["lag.ft.NWBIR74"],
            m3e$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m3e, 
               linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m3e, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m3e.imps[[1]]$total["ft.NWBIR74"],
            summary(m3e.imps)$semat["ft.NWBIR74", "Total"],
            m3e.imps[[1]]$total["east"],
            summary(m3e.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))

# In this setup, the total impact equals the direct impact if omitted from Durbin=formula
expect_true(isTRUE(all.equal(c(m3e.imps[[1]]$total["AREA"]), 
          c(m3e.imps[[1]]$direct["AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with unstandardized W and no intercept, and -1 in Durbin formula
m3f <- spatialreg::lmSLX(ft.SID74 ~ -1 + ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw, 
                         Durbin = ~ - 1 + ft.NWBIR74 + east)
# summary(m3f)
# summary(spatialreg::impacts(m3f))

# Impacts equal coefs
m3f.imps <- spatialreg::impacts(m3f)
expect_true(isTRUE(all.equal(c(m3f.imps[[1]]$indirect["ft.NWBIR74"],
            m3f.imps[[1]]$indirect["east"]), 
          c(m3f$coefficients["lag.ft.NWBIR74"],
            m3f$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m3f, 
               linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m3f, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m3f.imps[[1]]$total["ft.NWBIR74"],
            summary(m3f.imps)$semat["ft.NWBIR74", "Total"],
            m3f.imps[[1]]$total["east"],
            summary(m3f.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with unstandardized W and no intercept, no intercept declaration in Durbin formula
m3fb <- spatialreg::lmSLX(ft.SID74 ~ -1 + ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw, 
                         Durbin = ~ ft.NWBIR74 + east)
# summary(m3fb)
# summary(spatialreg::impacts(m3fb))

# Impacts equal coefs
m3fb.imps <- spatialreg::impacts(m3fb)
expect_true(isTRUE(all.equal(c(m3fb.imps[[1]]$indirect["ft.NWBIR74"],
            m3fb.imps[[1]]$indirect["east"]), 
          c(m3fb$coefficients["lag.ft.NWBIR74"],
            m3fb$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc1 <- summary(multcomp::glht(m3fb, 
               linfct = c("ft.NWBIR74 + lag.ft.NWBIR74 = 0")))
lc2 <- summary(multcomp::glht(m3fb, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m3fb.imps[[1]]$total["ft.NWBIR74"],
            summary(m3fb.imps)$semat["ft.NWBIR74", "Total"],
            m3fb.imps[[1]]$total["east"],
            summary(m3fb.imps)$semat["east", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))


#### --------- Tests for SDEM --------- ####

### SDEM with unstandardized W and intercept in lagX
m1b <- spatialreg::errorsarlm(ft.SID74 ~ ft.NWBIR74, 
                              weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
# summary(m1b)
# summary(spatialreg::impacts(m1b))

# Impact equals coef
expect_true(isTRUE(all.equal(spatialreg::impacts(m1b)[[1]]$indirect, 
                             m1b$coefficients["lag.ft.NWBIR74"], 
          tolerance = 1e-5, check.attributes = FALSE)))



### SDEM with unstandardized W, intercept in lagX, and multiple covariates
m1c <- spatialreg::errorsarlm(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                              weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
# summary(m1c)
# summary(spatialreg::impacts(m1c))

# Impacts equal coefs
m1c.imps <- spatialreg::impacts(m1c)
expect_true(isTRUE(all.equal(c(m1c.imps[[1]]$indirect["ft.NWBIR74"],
            m1c.imps[[1]]$indirect["east"],
            m1c.imps[[1]]$indirect["AREA"]), 
          c(m1c$coefficients["lag.ft.NWBIR74"],
            m1c$coefficients["lag.east"],
            m1c$coefficients["lag.AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m1c, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m1c.imps[[1]]$total["east"],
            summary(m1c.imps)$semat["east", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### SDEM with standardized W, without intercept in lagX, and multiple covariates
m1d <- spatialreg::errorsarlm(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                   weights=BIR74, data=nc, listw=nc_lw_w, Durbin=TRUE)
# summary(m1d)
# summary(spatialreg::impacts(m1d))

# Impacts equal coefs
m1d.imps <- spatialreg::impacts(m1d)
expect_true(isTRUE(all.equal(c(m1d.imps[[1]]$indirect["ft.NWBIR74"],
            m1d.imps[[1]]$indirect["east"],
            m1d.imps[[1]]$indirect["AREA"]), 
          c(m1d$coefficients["lag.ft.NWBIR74"],
            m1d$coefficients["lag.east"],
            m1d$coefficients["lag.AREA"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m1d, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m1d.imps[[1]]$total["east"],
            summary(m1d.imps)$semat["east", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))



### Test Durbin is formula, with standardized W
m1e <- spatialreg::errorsarlm(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                              weights=BIR74, data=nc, listw=nc_lw_w, 
                              Durbin = ~ ft.NWBIR74 + east)
# summary(m1e)
# summary(spatialreg::impacts(m1e))

# Impacts equal coefs
m1e.imps <- spatialreg::impacts(m1e)
expect_true(isTRUE(all.equal(c(m1e.imps[[1]]$indirect["ft.NWBIR74"],
            m1e.imps[[1]]$indirect["east"]), 
          c(m1e$coefficients["lag.ft.NWBIR74"],
            m1e$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m1e, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m1e.imps[[1]]$total["east"],
            summary(m1e.imps)$semat["east", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))


### Test Durbin is formula, with unstandardized W and intercept
# Not allowed, switches to Durbin = TRUE
m1f <- spatialreg::errorsarlm(ft.SID74 ~ ft.NWBIR74 + east + AREA, 
                         weights=BIR74, data=nc, listw=nc_lw, 
                         Durbin = ~ ft.NWBIR74 + east)
# summary(m1f)
# summary(spatialreg::impacts(m1f))
# 
# Impacts equal coefs
m1f.imps <- spatialreg::impacts(m1f)
expect_true(isTRUE(all.equal(c(m1f.imps[[1]]$indirect["ft.NWBIR74"],
            m1f.imps[[1]]$indirect["east"]), 
          c(m1f$coefficients["lag.ft.NWBIR74"],
            m1f$coefficients["lag.east"]), 
          tolerance = 1e-5, check.attributes = FALSE)))

# Total impact is linear combination
lc <- summary(multcomp::glht(m1f, linfct = c("east + lag.east = 0")))
expect_true(isTRUE(all.equal(c(m1f.imps[[1]]$total["east"],
            summary(m1f.imps)$semat["east", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)))




