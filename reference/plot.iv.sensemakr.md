# Sensitivity analysis plots for IV sensemakr

This function provides the contour plots of the sensitivity analysis
results obtained with the function
[`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md)
for IV. It is basically a dispatcher to the core plot function
[`ovb_contour_plot`](https://carloscinelli.com/iv.sensemakr/reference/ovb_contour_plot.md).

## Usage

``` r
# S3 method for class 'iv.sensemakr'
plot(x, sensitivity.of = c("ci", "lwr", "upr", "t-value"), parm = "iv", ...)
```

## Arguments

- x:

  an object of class `iv.sensemakr` created with the
  [`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md)
  function.

- sensitivity.of:

  should the contour plot show adjusted lower limits of confidence
  intervals (`"lwr"`), upper limit of confidence intervals (`"upr"`) or
  t-values (`"t-value"`)?

- parm:

  contour plots of which estimate? Options are `iv` for instrumental
  variable estimates, `fs` for first-stage estimates, and `rf` for
  reduced-form estimates.

- ...:

  further arguments and graphical parameters.
