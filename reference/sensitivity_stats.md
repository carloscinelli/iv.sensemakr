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
