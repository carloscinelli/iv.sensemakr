
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iv.sensemakr: Sensitivity Analysis Tools for IV

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/iv.sensemakr)](https://CRAN.R-project.org/package=iv.sensemakr)
[![R-CMD-check](https://github.com/carloscinelli/iv.sensemakr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/carloscinelli/iv.sensemakr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/carloscinelli/iv.sensemakr/graph/badge.svg)](https://app.codecov.io/gh/carloscinelli/iv.sensemakr)
<!-- badges: end -->

`iv.sensemakr` implements a suite of sensitivity analysis tools for
instrumental variable estimates, as discussed in Cinelli, C. and
Hazlett, C. (2025) “An Omitted Variable Bias Framework for Sensitivity
Analysis of Instrumental Variables”, *Biometrika*
([doi:10.1093/biomet/asaf004](https://doi.org/10.1093/biomet/asaf004);
[PDF](https://carloscinelli.com/files/Cinelli%20and%20Hazlett%20-%20OVB%20for%20IV.pdf)).

# News

- Paper published in Biometrika: Cinelli, C. and Hazlett, C. (2025), “An
  Omitted Variable Bias Framework for Sensitivity Analysis of
  Instrumental Variables.”
  [doi:10.1093/biomet/asaf004](https://doi.org/10.1093/biomet/asaf004)You
  can find a free version here
  ([PDF](https://carloscinelli.com/files/Cinelli%20and%20Hazlett%20-%20OVB%20for%20IV.pdf))
- Package website is now online
  <https://carloscinelli.com/iv.sensemakr/>
- Video of the [talk at
  PCIC](https://www.youtube.com/watch?v=ykZqPB1FPNo)
- Video of the [talk at
  JSM](https://www.youtube.com/watch?v=0b5VYWh8fN4)

# Development version

To install the development version on GitHub make sure you have the
package `devtools` installed.

``` r
# install.packages("devtools") 
devtools::install_github("carloscinelli/iv.sensemakr")
```

Please also make sure you have the latest version of `sensemakr`
installed.

``` r
# install.packages("devtools") 
devtools::install_github("carloscinelli/sensemakr")
```

# Basic usage

``` r
# loads package
library(iv.sensemakr)

# loads dataset
data("card")

# prepares data
y <- card$lwage  # outcome
d <- card$educ   # treatment
z <- card$nearc4 # instrument
x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
                   data = card) # covariates
# fits IV model
card.fit <- iv_fit(y,d,z,x)

# see results
card.fit
#> 
#> Instrumental Variable Estimation
#> (Anderson-Rubin Approach)
#> =============================================
#> IV Estimates:
#>   Coef. Estimate: 0.132
#>   t-value: 2.33
#>   p-value: 0.02
#>   Conf. Interval: [0.025, 0.285]
#> Note: H0 = 0, alpha = 0.05, df = 2994.
#> =============================================
#> See summary for first stage and reduced form.

# runs sensitivity analysis
card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

# see results
card.sens
#> 
#> Sensitivity Analysis for Instrumental Variables
#> (Anderson-Rubin Approach)
#> =============================================================
#> IV Estimates:
#>   Coef. Estimate: 0.132
#>   t-value: 2.33
#>   p-value: 0.02
#>   Conf. Interval: [0.025, 0.285]
#> 
#> Sensitivity Statistics:
#>   Extreme Robustness Value: 0.000523
#>   Robustness Value: 0.00667
#> 
#> Bounds on Omitted Variable Bias:
#>  Bound Label  R2zw.x R2y0w.zx Lower CI Upper CI Crit. Thr.
#>     1x black 0.00221   0.0750  -0.0212    0.402       2.59
#>      1x smsa 0.00639   0.0202  -0.0192    0.396       2.57
#> 
#> Note: H0 = 0, q >= 1, alpha = 0.05, df = 2994.
#> =============================================================
#> See summary for first stage and reduced form.

# sensitivity contour plot
plot(card.sens, lim = 0.09)
```

<img src="man/figures/README-basic-usage-1.png" width="100%" style="display: block; margin: auto;" />
