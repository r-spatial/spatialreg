# Version 1.2-9 (development)

* raised #39, no support for weights in SEM/SDEM/SLX #39

* address #37; #38 remains (no formula Durbin support for prediction using any Sarlm object)

* address #19 by not reporting `AIC` where case weights are used in `spautolm` or `errorsarlm`

* address bug in `predict()` for new data, SDEM. Others in #37, #38 need work.

* Further added checking for SLX/SDEM impacts and edge/corner cases; starting transition to use **multcomp** in place og **gmodels**

# Version 1.2-8 (2023-03-01)

* Attending to SLX/Durbin/non-W-style weights: #7, #36, #26, #35, #30 #24, #23, partly based on #13

# Version 1.2-6 (2022-10-07)

* make local copy of `gmodels::estimable()` for lm objects only, add authors to package contributors

* -Wstrict-prototypes fixes

# Version 1.2-5 (2022-08-16)

* updating coercion for **Matrix** 1.4-2

* updating dontrun/donttest for package split (previously unchecked, mostly in `aple`)

# Version 1.2-3 (2022-04-18)

* protect definition of USE_FC_LEN_T

* unescape underscores in help pages

# Version 1.2-1 (2021-11-11)

* Add Fortran character handling USE_FC_LEN_T WRE §6.6.1

* Add **spdep** split-out functionality

# Version 1.1-8 (2021-05-03)

* #18 standardize use of `coef()` methods for (some) fitted model summary objects

* https://github.com/tidymodels/broom/issues/1003#issuecomment-798694400 changing **spatialreg** model output class names: **spdep** `sarlm` -> **spatialreg** `Sarlm`, `spautolm` -> `Spautolm`, `stsls` -> `Stsls`, `gmsar` -> `Gmsar`, `lagmess` -> `Lagmess`, `SLX` -> , `SlX`, `MCMC_s*_g` -> `MCMC_s*_G`, `SFResult` -> `SfResult`, `ME_res` -> `Me_res`, `lagImpact` -> `LagImpact`, `WXImpact` -> `WXimpact`

* #16 merged coordination of impacts methods (Gianfranco Piras)

* #14 merged correction to SDEM and SLX impacts when a lagged intercept is present (Tobias Rüttenauer).

# Version 1.1-5 (2019-12-01)

* #6, #11 na.action and precomputed eigenvalue bug

* #9 Griddy Gibbs issue

* #8 Predict method for SLX

* #7, #13-14 Offset impacts for SDEM/SLX

* #5, #10 Panel Durbin= argument


# Version 1.1-3 (2019-04-01)

* #2 Split spatialreg from spdep; spdep functions still present there as deprecated

