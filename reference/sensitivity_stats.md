# Sensitivity statistics for instrumental variable estimates

Convenience function that computes robustness values for IV estimates as
well as auxiliary first stage and reduced form regressions.

## Usage

``` r
sensitivity_stats(...)

# S3 method for class 'iv_fit'
sensitivity_stats(model, parm = "iv", q = 1, alpha = 0.05, min = TRUE, ...)

# S3 method for class 'iv.sensemakr'
sensitivity_stats(model, parm = "iv", q = 1, alpha = 0.05, min = TRUE, ...)
```

## Arguments

- ...:

  further arguments passed to or from other methods.

- model:

  a model created with the function
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

- parm:

  contour plots of which estimate? Options are `iv` for instrumental
  variable estimates, `fs` for first-stage estimates, and `rf` for
  reduced-form estimates.

- q:

  percent change of the effect estimate that would be deemed
  problematic. Default is 1, which means a reduction of 100% of the
  current effect estimate (bring estimate to zero).

- alpha:

  significance level.

- min:

  should we consider biases as large or larger than a certain amount?
  Default is `TRUE`.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with columns
for the estimate, confidence interval bounds (lower and upper), t-value,
extreme robustness value (`xrv_qa`), robustness value (`rv_qa`), and the
parameters used (`q`, `min`, `alpha`, `dof`).

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

# sensitivity statistics for the IV estimate
sensitivity_stats(card.fit)
#>     estimate        lwr       upr  t.value       xrv_qa       rv_qa q min alpha
#> iv 0.1315038 0.02480484 0.2848236 2.327075 0.0005232443 0.006666407 1   1  0.05
#>     dof
#> iv 2994

# sensitivity statistics for the first-stage
sensitivity_stats(card.fit, parm = "fs")
#>     estimate       lwr       upr t.value      xrv_qa      rv_qa q min alpha
#> fs 0.3198989 0.1476194 0.4921785 3.64085 0.003129076 0.03023129 1   1  0.05
#>     dof
#> fs 2994
```
