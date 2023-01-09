# **sawnuti** 
  
![](https://www.r-pkg.org/badges/version/sawnuti) ![](https://www.r-pkg.org/badges/last-release/sawnuti)

An algorithm for the comparison of two sequences in time.  This package implements the methods introduced in [Murph et al. 2021](https://www.tandfonline.com/doi/full/10.1080/07474946.2021.1940491).

This repository is organized as a stand-alone R package.  For questions, issues, or clarifications please reach out to Murph: <acmurph@unc.edu>.  Feel free to email any applications; we'd be happy to highlight them here.


## Installation

You can install the latest version from CRAN using:

``` r
install.packages( "sawnuti" )
```

``` r
require( "sawnuti" )
```

## Examples

```r
matchFunction = function(a,b){ifelse(a==b, 1, -1)}

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