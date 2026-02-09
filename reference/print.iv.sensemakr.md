# Sensitivity analysis print and summary methods for `iv.sensemakr`

The `print` and `summary` methods provide verbal descriptions of the
sensitivity analysis results obtained with the function
[`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md).

## Usage

``` r
# S3 method for class 'iv.sensemakr'
summary(object, ...)

# S3 method for class 'summary.iv.sensemakr'
print(x, digits = 3, ...)
```

## Arguments

- object:

  an object of class
  [`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md).

- ...:

  arguments passed to other methods.

- x:

  an object of class
  [`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md).

- digits:

  minimal number of *significant* digits.

## Value

`summary.iv.sensemakr` returns an object of class
`summary.iv.sensemakr`. The `print` methods return the object invisibly.

## Examples

``` r
data("card")
y <- card$lwage
d <- card$educ
z <- card$nearc4
x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
                   data = card)
card.fit <- iv_fit(y, d, z, x)
card.sens <- sensemakr(card.fit, benchmark_covariates = "black")
print(card.sens)
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
#>     1x black 0.00221    0.075  -0.0212    0.402       2.59
#> 
#> Note: H0 = 0, q >= 1, alpha = 0.05, df = 2994.
#> =============================================================
#> See summary for first stage and reduced form.
#> 
summary(card.sens)
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
#>  Bound Label  R2zw.x r2y0w.zx Lower CI Upper CI Crit. Thr.
#>     1x black 0.00221    0.075  -0.0212    0.402       2.59
#> 
#> Note: H0 = 0, q >= 1, alpha = 0.05, df = 2994.
#> -------------------------------------------------------------
#> FS Estimates:
#>   Coef. Estimate: 0.32
#>   Standard Error: 0.0879
#>   t-value: 3.64
#>   p-value: 0.000276
#>   Conf. Interval: [0.148, 0.492]
#> 
#> Sensitivity Statistics:
#>   Extreme Robustness Value: 0.00313
#>   Robustness Value: 0.0302
#> 
#> Bounds on Omitted Variable Bias:
#>  Bound Label  R2zw.x R2dw.zx Lower CI Upper CI Crit. Thr.
#>     1x black 0.00221  0.0334    0.109    0.531        2.4
#> 
#> Note: H0 = 0, q = 1, alpha = 0.05, df = 2994.
#> -------------------------------------------------------------
#> RF Estimates:
#>   Coef. Estimate: 0.0421
#>   Standard Error: 0.0181
#>   t-value: 2.33
#>   p-value: 0.02
#>   Conf. Interval: [0.007, 0.078]
#> 
#> Sensitivity Statistics:
#>   Extreme Robustness Value: 0.000523
#>   Robustness Value: 0.00667
#> 
#> Bounds on Omitted Variable Bias:
#>  Bound Label  R2zw.x R2yw.zx Lower CI Upper CI Crit. Thr.
#>     1x black 0.00221  0.0657 -0.00418   0.0883       2.56
#> 
#> Note: H0 = 0, q = 1, alpha = 0.05, df = 2994.
#> =============================================================
#> 
```
