# National Longitudinal Survey of Young Men (NLSYM)

Data used in Card (1995). Consists of a sample of 3,010 individuals from
the National Longitudinal Survey of Young Men (NLSYM).

The treatment is `educ`, the outcome is `lwage` and the instrument is
`nearc4`.

## Usage

``` r
data('card')
```

## Format

A data.frame with 3010 observations and 34 variables:

- **id:** person identifier

- **nearc2:** =1 if near 2 yr college, 1966

- **nearc4:** =1 if near 4 yr college, 1966

- **educ:** years of schooling, 1976

- **age:** in years

- **fatheduc:** father's schooling

- **motheduc:** mother's schooling

- **weight:** NLS sampling weight, 1976

- **momdad14:** =1 if live with mom, dad at 14

- **sinmom14:** =1 if with single mom at 14

- **step14:** =1 if with step parent at 14

- **reg661:** =1 for region 1, 1966

- **reg662:** =1 for region 2, 1966

- **reg663:** =1 for region 3, 1966

- **reg664:** =1 for region 4, 1966

- **reg665:** =1 for region 5, 1966

- **reg666:** =1 for region 6, 1966

- **reg667:** =1 for region 7, 1966

- **reg668:** =1 for region 8, 1966

- **reg669:** =1 for region 9, 1966

- **south66:** =1 if in south in 1966

- **black:** =1 if black

- **smsa:** =1 in in SMSA, 1976

- **south:** =1 if in south, 1976

- **smsa66:** =1 if in SMSA, 1966

- **wage:** hourly wage in cents, 1976

- **enroll:** =1 if enrolled in school, 1976

- **KWW:** knowledge world of work score

- **IQ:** IQ score

- **married:** =1 if married, 1976

- **libcrd14:** =1 if lib. card in home at 14

- **exper:** age - educ - 6

- **lwage:** log(wage)

- **expersq:** exper^2

## References

Card, D. "Using Geographic Variation in College Proximity to Estimate
the Return to Schooling". In L.N. Christofides, E.K. Grant, and R.
Swidinsky, editors, Aspects of Labor Market Behaviour: Essays in Honour
of John Vanderkamp. Toronto: University of Toronto Press, 1995.

## Examples

``` r
data('card')
head(card)
#>   id nearc2 nearc4 educ age fatheduc motheduc weight momdad14 sinmom14 step14
#> 1  2      0      0    7  29       NA       NA 158413        1        0      0
#> 2  3      0      0   12  27        8        8 380166        1        0      0
#> 3  4      0      0   12  34       14       12 367470        1        0      0
#> 4  5      1      1   11  27       11       12 380166        1        0      0
#> 5  6      1      1   12  34        8        7 367470        1        0      0
#> 6  7      1      1   12  26        9       12 380166        1        0      0
#>   reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669 south66 black
#> 1      1      0      0      0      0      0      0      0      0       0     1
#> 2      1      0      0      0      0      0      0      0      0       0     0
#> 3      1      0      0      0      0      0      0      0      0       0     0
#> 4      0      1      0      0      0      0      0      0      0       0     0
#> 5      0      1      0      0      0      0      0      0      0       0     0
#> 6      0      1      0      0      0      0      0      0      0       0     0
#>   smsa south smsa66 wage enroll KWW  IQ married libcrd14 exper    lwage expersq
#> 1    1     0      1  548      0  15  NA       1        0    16 6.306275     256
#> 2    1     0      1  481      0  35  93       1        1     9 6.175867      81
#> 3    1     0      1  721      0  42 103       1        1    16 6.580639     256
#> 4    1     0      1  250      0  25  88       1        1    10 5.521461     100
#> 5    1     0      1  729      0  34 108       1        0    16 6.591674     256
#> 6    1     0      1  500      0  38  85       1        1     8 6.214608      64
```
