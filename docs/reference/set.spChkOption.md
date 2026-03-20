# Control checking of spatial object IDs

Provides support for checking the mutual integrity of spatial neighbour
weights and spatial data; similar mechanisms are used for passing global
verbose and zero.policy options, and for providing access to a running
cluster for embarrassingly parallel tasks.

## Usage

``` r
set.VerboseOption(check)
get.VerboseOption()
set.ZeroPolicyOption(check)
get.ZeroPolicyOption()
#set.listw_is_CsparseMatrix_Option(check)
#get.listw_is_CsparseMatrix_Option()
```

## Arguments

- check:

  a logical value, TRUE or FALSE

## Details

Analysis functions will have an spChk argument by default set to NULL,
and will call
[`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
to get the global spatial option for whether to check or not — this is
initialised to FALSE, and consequently should not break anything. It can
be changed to TRUE using `set.spChkOption(TRUE)`, or the spChk argument
can be assigned in analysis functions.
[`spNamedVec()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
is provided to ensure that rownames are passed on to single columns
taken from two-dimensional arrays and data frames.

## Value

[`set.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
returns the old logical value,
[`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
returns the current logical value, and
[`chkIDs()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
returns a logical value for the test lack of difference.
[`spNamedVec()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.html)
returns the selected column with the names set to the row names of the
object from which it has been extracted.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Examples

``` r
get.VerboseOption()
#> [1] FALSE
get.ZeroPolicyOption()
#> [1] FALSE
```
