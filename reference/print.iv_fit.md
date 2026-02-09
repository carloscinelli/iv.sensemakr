# print and summary methods for `iv_fit`

The print and summary methods provide verbal descriptions of the results
obtained with the function
[`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

## Usage

``` r
# S3 method for class 'iv_fit'
summary(object, ...)

# S3 method for class 'iv_fit'
print(x, digits = 3, ...)
```

## Arguments

- object:

  an object of class
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

- ...:

  arguments passed to other methods.

- x:

  an object of class
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

- digits:

  minimal number of significant digits

## Value

`print.iv_fit` returns the object `x` invisibly. `summary.iv_fit`
returns an object of class `summary.iv_fit`. `print.summary.iv_fit`
returns its argument invisibly.

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
print(card.fit)
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
#> 
summary(card.fit)
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
#> ---------------------------------------------
#> FS Estimates:
#>   Coef. Estimate: 0.32
#>   Standard Error: 0.0879
#>   t-value: 3.64
#>   p-value: 0.000276
#>   Conf. Interval: [0.148, 0.492]
#> Note: H0 = 0, alpha = 0.05, df = 2994.
#> ---------------------------------------------
#> RF Estimates:
#>   Coef. Estimate: 0.0421
#>   Standard Error: 0.0181
#>   t-value: 2.33
#>   p-value: 0.02
#>   Conf. Interval: [0.007, 0.078]
#> Note: H0 = 0, alpha = 0.05, df = 2994.
#> =============================================
#> 
```
