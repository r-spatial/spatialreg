# spatialreg
spatialreg: spatial models estimation and testing

(name of package and repo subject to change)

This new project aims at collecting all the estimation functions for spatial cross-sectional models (on lattice) currently contained in spdep, sphet and spse.

For now only has functions from spdep, where they are shown as deprecated. This package only loads the namespace of spdep; if you attach it, the same functions in the other packages will be masked. Some feed through adequately, others do not (mostly where model.matrix facilities do not like the extra level of passing arguments).
