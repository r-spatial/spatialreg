# Options for parallel support

Provides support for the use of parallel computation in the parallel
package.

## Usage

``` r
set.mcOption(value)
get.mcOption()
set.coresOption(value)
get.coresOption()
set.ClusterOption(cl)
get.ClusterOption()
```

## Arguments

- value:

  valid replacement value

- cl:

  a cluster object created by `makeCluster` in parallel

## Details

Options in the spatialreg package are held in an environment local to
the package namespace and not exported. Option values are set and
retrieved with pairs of access functions, get and set. The `mc` option
is set by default to FALSE on Windows systems, as they cannot fork the R
session; by default it is TRUE on other systems, but may be set FALSE.
If `mc` is FALSE, the `Cluster` option is used: if `mc` is FALSE and the
`Cluster` option is NULL no parallel computing is done, or the `Cluster`
option is passed a “cluster” object created by the parallel or snow
package for access without being passed as an argument. The `cores`
option is set to NULL by default, and can be used to store the number of
cores to use as an integer. If `cores` is NULL, facilities from the
parallel package will not be used.

## Value

The option access functions return their current settings, the
assignment functions usually return the previous value of the option.

## Note

An extended example is shown in the documentation of
[`mom_calc`](https://r-spatial.github.io/spatialreg/reference/trW.md),
including treatment of seeding of RNG for multicore/cluster.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Examples

``` r
ls(envir=spatialreg:::.spatialregOptions)
#> [1] "cluster"      "cores"        "mc"           "rlecuyerSeed" "verbose"     
#> [6] "zeroPolicy"  
library(parallel)
nc <- max(2L, detectCores(logical=FALSE), na.rm = TRUE)-1L
nc
#> [1] 5
# set nc to 1L here
if (nc > 1L) nc <- 1L
#nc <- ifelse(nc > 2L, 2L, nc)
coresOpt <- get.coresOption()
coresOpt
#> NULL
if (!is.na(nc)) {
 invisible(set.coresOption(nc))
 print(exists("mom_calc"))
 if(.Platform$OS.type == "windows") {
# forking not permitted on Windows - start cluster
# removed for Github actions 210502
# \dontrun{
  print(get.mcOption())
  cl <- makeCluster(get.coresOption())
  print(clusterEvalQ(cl, exists("mom_calc")))
  set.ClusterOption(cl)
  clusterEvalQ(get.ClusterOption(), library(spatialreg))
  print(clusterEvalQ(cl, exists("mom_calc")))
  clusterEvalQ(get.ClusterOption(), detach(package:spatialreg))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("mom_calc")))
  stopCluster(cl)
# }
 } else {
  mcOpt <- get.mcOption()
  print(mcOpt)
  print(mclapply(1:get.coresOption(), function(i) exists("mom_calc"),
   mc.cores=get.coresOption()))
  invisible(set.mcOption(FALSE))
  cl <- makeCluster(nc)
  print(clusterEvalQ(cl, exists("mom_calc")))
  set.ClusterOption(cl)
  clusterEvalQ(get.ClusterOption(), library(spatialreg))
  print(clusterEvalQ(cl, exists("mom_calc")))
  clusterEvalQ(get.ClusterOption(), detach(package:spatialreg))
  set.ClusterOption(NULL)
  print(clusterEvalQ(cl, exists("mom_calc")))
  stopCluster(cl)
  invisible(set.mcOption(mcOpt))
 }
 invisible(set.coresOption(coresOpt))
}
#> [1] TRUE
#> [1] TRUE
#> [[1]]
#> [1] TRUE
#> 
#> [[1]]
#> [1] FALSE
#> 
#> [[1]]
#> [1] TRUE
#> 
#> [[1]]
#> [1] FALSE
#> 
```
