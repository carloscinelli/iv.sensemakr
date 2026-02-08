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
