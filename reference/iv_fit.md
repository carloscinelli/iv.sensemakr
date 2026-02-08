# Instrumental Variable Estimation using the Anderson-Rubin approach

`iv_fit` computes instrumental variable estimates and confidence
intervals using the Anderson-Rubin (AR) approach (Anderson and Rubin,
1949). This approach is numerically identical to Fieller's theorem
(Fieller, 1954). See Cinelli and Hazlett (2025) for further discussion.

The AR point estimate is numerically identical to the point estimate of
two-stage least squares (2SLS) and it is given by the ratio of the
reduced-form to the first-stage regression coefficient. Confidence
intervals, however, are constructed differently. 2SLS is equivalent to
using the delta-method to obtain the variance of the ratio estimator,
and then proceeding by assuming the ratio is asymptotically normal. This
approximation can fail when instruments are "weak." The Anderson-Rubin
approach instead uses a test inversion procedure to construct confidence
intervals. This procedure has correct coverage regardless of instrument
strength, at the (inevitable) cost of eventually obtaining unbounded
confidence intervals.

## Usage

``` r
iv_fit(y, d, z, x = NULL, h0 = 0, alpha = 0.05)
```

## Arguments

- y:

  [`numeric`](https://rdrr.io/r/base/numeric.html)
  [`vector`](https://rdrr.io/r/base/vector.html) with the outcome.

- d:

  [`numeric`](https://rdrr.io/r/base/numeric.html)
  [`vector`](https://rdrr.io/r/base/vector.html) with the treatment.

- z:

  [`numeric`](https://rdrr.io/r/base/numeric.html)
  [`vector`](https://rdrr.io/r/base/vector.html) with the instrument.

- x:

  (optional) [`numeric`](https://rdrr.io/r/base/numeric.html)
  [`matrix`](https://rdrr.io/r/base/matrix.html) with observed
  covariates.

- h0:

  null hypothesis for the target parameter (the IV estimate).

- alpha:

  significance test for hypothesis tests and confidence intervals.

## Value

An object of class `iv_fit`, containing:

- `data` :

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with the data
  used for fitting the models.

- `models` :

  A [`list`](https://rdrr.io/r/base/list.html) with the
  [`lm`](https://rdrr.io/r/stats/lm.html) models used for obtaining the
  IV estimates. This includes the first-stage (FS), reduced-form (RF),
  and Anderson-Rubin (AR) regressions.

- `estimates` :

  A [`list`](https://rdrr.io/r/base/list.html) with the summary
  information of IV estimates, as well as summary information of the
  auxiliary estimates of the FS, RF, and AR regression.

- `pars`:

  A [`list`](https://rdrr.io/r/base/list.html) with the parameters of
  the call, such as the null hypothesis `h0` and the significance level
  `alpha`.

@returns The function returns an object of class `iv_fit`.

## References

Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a
single equation in a complete system of stochastic equations, Annals of
Mathematical Statistics, 20, 46-63.

Fieller, E. C. (1954). Some problems in interval estimation. Journal of
the Royal Statistical Society: Series B (Methodological), 16(2),
175-185.

Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework
for Sensitivity Analysis of Instrumental Variables." Biometrika.
[doi:10.1093/biomet/asaf004](https://doi.org/10.1093/biomet/asaf004)
