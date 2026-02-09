# Computes the (extreme) robustness value for IV

Computes robustness values for
[`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md)
objects, adapting the robustness value definitions of sensemakr to
instrumental variables as described in Cinelli and Hazlett (2025). For
`parm = "iv"`, returns robustness values for the IV estimate; for
`parm = "fs"` or `parm = "rf"`, dispatches to sensemakr methods on the
corresponding [`lm`](https://rdrr.io/r/stats/lm.html) models.

## Usage

``` r
robustness_value(...)

extreme_robustness_value(...)

xrv(...)

rv(...)

# S3 method for class 'iv_fit'
extreme_robustness_value(
  model,
  parm = "iv",
  q = 1,
  alpha = 0.05,
  min = TRUE,
  ...
)

# S3 method for class 'iv_fit'
robustness_value(model, parm = "iv", q = 1, alpha = 0.05, min = TRUE, ...)
```

## Arguments

- ...:

  further arguments passed to or from other methods.

- model:

  an
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md)
  model

- parm:

  parameter for which the robustness value is computed. Default is `iv`,
  meaning that the robustness value of the IV estimate is computed.
  Other options are to compute the robustness value of auxiliary
  estimates, such as the first stage (`fs`) or the reduced form (`rf`).

- q:

  percent change of the effect estimate that would be deemed
  problematic. Default is 1, which means a reduction (increase) of 100%
  of the current effect estimate (bring estimate to zero). It has to be
  greater than zero.

- alpha:

  significance level.

- min:

  in many cases, researchers are interested in biases as large or larger
  than a certain amount (for instance, the strength of confounding to
  bring a positive estimate to zero or below). Setting `min = TRUE`
  (default) computes the robustness value for such cases. Setting
  `min = FALSE` computes the robustness value for a bias of exactly `q`.

## References

Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework
for Sensitivity Analysis of Instrumental Variables." Biometrika.
[doi:10.1093/biomet/asaf004](https://doi.org/10.1093/biomet/asaf004)
