# Extracts point estimates and confidence intervals of an `iv_fit` model.

The function `coef` extracts point estimates of an
[`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md)
model.

The function `confint` extracts confidence intervals of an
[`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md)
model.

## Usage

``` r
# S3 method for class 'iv_fit'
coef(object, parm = "iv", ...)

# S3 method for class 'iv_fit'
confint(object, parm = c("iv", "fs", "rf"), level, ...)
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

- level:

  coverage level (i.e, 1-alpha). If not provided, it uses the same level
  as the one provided in
  [`iv_fit`](https://carloscinelli.com/iv.sensemakr/reference/iv_fit.md).

## Value

`coef` returns a numeric vector with the estimates of interest.

`confint` returns a numeric vector with the confidence interval of
interest.

## Examples

``` r
# prepare data
data("card")
y <- card$lwage
d <- card$educ
z <- card$nearc4
x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
                   data = card)

# fit iv model
card.fit <- iv_fit(y, d, z, x)

# extract coefficients
coef(card.fit)
#>        iv 
#> 0.1315038 
coef(card.fit, parm = "fs")
#>        fs 
#> 0.3198989 
coef(card.fit, parm = "rf")
#>         rf 
#> 0.04206794 

# extract confidence intervals
confint(card.fit)
#> IV conf. interval:
#> [0.025, 0.285]
#> level = 0.95
confint(card.fit, parm = "fs")
#>        2.5 %    97.5 %
#> fs 0.1476194 0.4921785
confint(card.fit, parm = "rf")
#>          2.5 %     97.5 %
#> rf 0.006622162 0.07751371
```
