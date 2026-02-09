# Extract estimates of an `iv.sensemakr` object

This function extracts the estimate, lower limit, upper limit, t-value,
and (extreme) robustness values of an `iv.sensemakr` object, created
with the function
[`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md).

## Usage

``` r
# S3 method for class 'iv.sensemakr'
coef(object, parm = "iv", ...)
```

## Arguments

- object:

  an object of class
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

- parm:

  which estimate to return. Options are `"iv"` for instrumental variable
  estimate, `"fs"` for the first-stage estimate and `"rf"` for the
  reduced-form estimate.

- ...:

  arguments passed to other methods.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with the
sensitivity statistics for the requested parameters.

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
coef(card.sens)
#>     estimate        lwr       upr  t.value       xrv_qa       rv_qa q min alpha
#> iv 0.1315038 0.02480484 0.2848236 2.327075 0.0005232443 0.006666407 1   1  0.05
#>     dof
#> iv 2994
coef(card.sens, parm = "fs")
#>     estimate       lwr       upr t.value      xrv_qa      rv_qa q min alpha
#> fs 0.3198989 0.1476194 0.4921785 3.64085 0.003129076 0.03023129 1   1  0.05
#>     dof
#> fs 2994
```
