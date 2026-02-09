# Minimal sensitivity reporting for IV estimates

This function produces LaTeX or HTML code for a minimal sensitivity
analysis table for instrumental variable estimates, as suggested in
Cinelli and Hazlett (2025). For objects of class
[`sensemakr`](https://carloscinelli.com/iv.sensemakr/reference/sensemakr.md)
(from the sensemakr package), it dispatches to
`sensemakr::`[`ovb_minimal_reporting`](https://rdrr.io/pkg/sensemakr/man/print.sensemakr.html).

## Usage

``` r
ovb_minimal_reporting(
  x,
  digits = 3,
  verbose = TRUE,
  format = c("latex", "html", "pure_html"),
  ...
)
```

## Arguments

- x:

  an object of class `iv.sensemakr` or `sensemakr`.

- digits:

  minimal number of *significant* digits.

- verbose:

  if `TRUE`, the function prints the code with
  [`cat`](https://rdrr.io/r/base/cat.html).

- format:

  code format to print: `"latex"`, `"html"` (requires mathjax), or
  `"pure_html"`.

- ...:

  further arguments passed to the table-building functions. Optional
  overrides include `outcome_label` and `treatment_label`.

## Value

The function returns the LaTeX or HTML code invisibly as a character
string and also prints it with [`cat`](https://rdrr.io/r/base/cat.html)
when `verbose = TRUE`.

## References

Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework
for Sensitivity Analysis of Instrumental Variables." Biometrika.
[doi:10.1093/biomet/asaf004](https://doi.org/10.1093/biomet/asaf004)

## Examples

``` r
# loads package
library(iv.sensemakr)

# loads dataset
data("card")

# prepares data
y <- card$lwage
d <- card$educ
z <- card$nearc4
x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
                   data = card)

# fits IV model and runs sensitivity analysis
card.fit <- iv_fit(y, d, z, x)
card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

# latex code
ovb_minimal_reporting(card.sens)
#> \begin{table}[!h]
#> \centering
#> \begin{tabular}{lrrrrrr}
#> \multicolumn{7}{c}{Outcome: \textit{y}} \\
#> \hline \hline 
#> Treatment: & Est. & Lower CI & Upper CI & t-value & $XRV_{q = 1, \alpha = 0.05}$ & $RV_{q = 1, \alpha = 0.05}$  \\ 
#> \hline 
#> \textit{d} & 0.132 & 0.025 & 0.285 & 2.327 & 0.1\% & 0.7\% \\ 
#> \hline 
#> df = 2994 & & \multicolumn{5}{r}{ \small \textit{Bound (1x black)}: $R^2_{Z\sim W| {\bf X}}$ = 0.2\%, $R^2_{Y(0)\sim W| Z, {\bf X}}$ = 7.5\%} \\
#> \end{tabular}
#> \end{table}

# html code (pure html, no mathjax needed)
ovb_minimal_reporting(card.sens, format = "pure_html")
#> <table style='align:center'>
#> <thead>
#> <tr>
#>  <th style="text-align:left;border-bottom: 1px solid transparent;border-top: 1px solid black"> </th>
#>  <th colspan = 6 style="text-align:center;border-bottom: 1px solid black;border-top: 1px solid black"> Outcome: y</th>
#> </tr>
#> <tr>
#>  <th style="text-align:left;border-top: 1px solid black"> Treatment </th>
#>  <th style="text-align:right;border-top: 1px solid black"> Est. </th>
#>  <th style="text-align:right;border-top: 1px solid black"> Lower CI </th>
#>  <th style="text-align:right;border-top: 1px solid black"> Upper CI </th>
#>  <th style="text-align:right;border-top: 1px solid black"> t-value </th>
#>  <th style="text-align:right;border-top: 1px solid black"> XRV<sub>q = 1, &alpha; = 0.05</sub> </th>
#>  <th style="text-align:right;border-top: 1px solid black"> RV<sub>q = 1, &alpha; = 0.05</sub> </th>
#> </tr>
#> </thead>
#> <tbody>
#>  <tr>
#>  <td style="text-align:left; border-bottom: 1px solid black"><i>d</i></td>
#>  <td style="text-align:right;border-bottom: 1px solid black">0.132 </td>
#>  <td style="text-align:right;border-bottom: 1px solid black">0.025 </td>
#>  <td style="text-align:right;border-bottom: 1px solid black">0.285 </td>
#>  <td style="text-align:right;border-bottom: 1px solid black">2.327 </td>
#>  <td style="text-align:right;border-bottom: 1px solid black">0.1% </td>
#>  <td style="text-align:right;border-bottom: 1px solid black">0.7% </td>
#> </tr>
#> </tbody>
#> <tr>
#> <td colspan = 7 style='text-align:right;border-bottom: 1px solid transparent;font-size:11px'>Note: df = 2994; Bound ( 1x black ):  R<sup>2</sup><sub>Z~W|X</sub> = 0.2%, R<sup>2</sup><sub>Y(0)~W|Z,X</sub> = 7.5%</td>
#> </tr>
#> </table>
```
