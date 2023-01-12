# **bayesWatch** 
  
![](https://www.r-pkg.org/badges/version/bayesWatch) ![](https://www.r-pkg.org/badges/last-release/bayesWatch)

EXPLAIN IT BETTER SOON THANK YOU

This repository is organized as a stand-alone R package.  For questions, issues, or clarifications please reach out to Murph: <acmurph@unc.edu>.  Feel free to email any applications; we'd be happy to highlight them here.


## Installation

You can install the latest version from CRAN using:

``` r
install.packages( "bayesWatch" )
```

``` r
require( "bayesWatch" )
```

## Examples

```r
full_data = data("example_data.rds")

day_of_observations = readRDS("bayesWatch/data/day_of_observations.rds")

day_dts = readRDS("bayesWatch/data/day_dts.rds")
x       = fit_regime_vector(full_data, day_of_observations, day_dts, 
                            iterations = 500, g.prior = 1, linger_parameter = 20, n.cores=3,
                            wishart_df_inital = 3, hyperprior_b = 3, lambda = 5)
saveRDS(x,"bayesWatch_fit_object.rds")
print(x)

xx = readRDS("bayesWatch_fit_object.rds")
detect_faults(xx)

sawnuti(string1="a b c", string2="d b c", times1="1 2 3",times2="3 2 1", alpha = 1, 
        match_function = matchFunction, gap_penalty = 1)
# $ScoreingMatrix
#   [,1] [,2] [,3] [,4]
# [1,]    0   -3   -5   -6
# [2,]   -3   -1   -4   -5
# [3,]   -1   -1    0   -2
# [4,]   -4   -4   -3   -1
#
# $AlignmentScore
# [1] "-1"
#
# $Alignment
#   [,1] [,2] [,3]
# [1,] "d"  "b"  "c"
# [2,] "|"  "|"  "|"
# [3,] "a"  "b"  "c"
```

## Packages Required

None.

## Citation

A. Murph, A. Flynt, B. R. King (2021). Comparing finite sequences of discrete events with non-uniform time intervals, <em>Sequential Analysis</em>, 40(3), 291-313.