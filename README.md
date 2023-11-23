
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TSL.logconcave

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the released version of TSL.logconcave from
[GitHub](https://github.com/) with:

``` r
#Uninstalling previous versions, if installed
devtools::install_github("nilanjanalaha/TSL.logconcave", force=TRUE)
#> Downloading GitHub repo nilanjanalaha/TSL.logconcave@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/ss/v9sfqr114lsfx9wvjl3397zc0000gr/T/RtmpxxL8TZ/remotes15bfc239201a7/nilanjanalaha-TSL.logconcave-f60c489/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/ss/v9sfqr114lsfx9wvjl3397zc0000gr/T/RtmpxxL8TZ/remotes15bfc239201a7/nilanjanalaha-TSL.logconcave-f60c489/DESCRIPTION’
#>   ─  preparing ‘TSL.logconcave’:
#>      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  building ‘TSL.logconcave_0.0.0.9000.tar.gz’
#>      
#> 
#> Installing package into '/Users/nlaha/Library/R/arm64/4.3/library'
#> (as 'lib' is unspecified)
```

This package calculates an estimator for estimating the location shift
in a two sample location model. Suppose $$X=\mu+Z_1$$ and
$$Y=\mu+\Delta+Z_2$$ where $Z_1$ and $Z_2$ are two i.i.d. random
variables with common density $g_0$. We also assume that $g_0$ is
log-concave.

Let us simulate some Toy data from the above model. We will take
$\mu=0$, $\Delta=0.50$, and we take $g_0$ to be the standard normal
density.

``` r
library(TSL.logconcave)
x <- rnorm(100, 0, 1)
y <- rnorm(50, 0.5 ,1) 
```

Note that the two samples can have different sample sizes. Then we
create a list called datam which represents the dataset.

``` r
datam <- list(x=x, y=y)
```

## The one step estimator

Our estimator is a one step estimator, which can be calculated using the
function pooled.smoothed. The output is a vector, whose first element is
the estimate of $\Delta$, and second element is an estimator of the
Fisher information for estimating $\Delta$. Type ?pooled_smoothed in the
R console for more details.

``` r
temp <- pooled_smoothed(datam, eta = 0.001)
# estimator
temp[1]
#> [1] 0.1248414
```

Here eta is the parameter which controls the truncation level. Our
suggestion is to use the default level of eta, which is 0. Simulations
show that this is the most efficient estimator in terms of asymptotic
variance.

``` r
temp <- pooled_smoothed(datam)
# estimator
temp[1]
#> [1] 0.1245371
```
