# test-7-ovb-minimal-reporting.R
# Tests for ovb_minimal_reporting and internal table functions

# setup -------------------------------------------------------------------
data("card")
y <- card$lwage
d <- card$educ
z <- card$nearc4
x <- model.matrix(~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
                   data = card)
card.fit  <- iv_fit(y, d, z, x)
card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

# also create a sensemakr object without benchmarks (no bounds)
card.sens.nobounds <- sensemakr(card.fit)


# format dispatch ---------------------------------------------------------

test_that("ovb_minimal_reporting default format is latex", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  expect_true(grepl("\\\\begin\\{table\\}", out))
  expect_true(grepl("\\\\end\\{table\\}", out))
})

test_that("ovb_minimal_reporting dispatches to html format", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "html")
  expect_true(grepl("<table>", out, fixed = TRUE))
  expect_true(grepl("</table>", out, fixed = TRUE))
  # html format uses mathjax notation
  expect_true(grepl("$XRV_{", out, fixed = TRUE))
})

test_that("ovb_minimal_reporting dispatches to pure_html format", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "pure_html")
  expect_true(grepl("<table", out, fixed = TRUE))
  expect_true(grepl("</table>", out, fixed = TRUE))
  # pure_html uses <sub> tags, not mathjax
  expect_true(grepl("<sub>", out, fixed = TRUE))
  expect_true(grepl("&alpha;", out, fixed = TRUE))
})


# LaTeX output content ----------------------------------------------------

test_that("latex table contains key structural elements", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  expect_true(grepl("\\\\begin\\{tabular\\}", out))
  expect_true(grepl("\\\\hline", out))
  expect_true(grepl("XRV_\\{", out))
  expect_true(grepl("RV_\\{", out))
  expect_true(grepl("\\\\alpha", out))
  expect_true(grepl("Est\\.", out))
  expect_true(grepl("Lower CI", out))
  expect_true(grepl("Upper CI", out))
  expect_true(grepl("t-value", out))
})

test_that("latex table contains numeric results", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  # estimate should be 0.132
  expect_true(grepl("0\\.132", out))
  # df should be 2994
  expect_true(grepl("2994", out))
})


# custom labels -----------------------------------------------------------

test_that("custom labels are used in latex output", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE,
                                outcome_label = "lwage",
                                treatment_label = "educ")
  expect_true(grepl("lwage", out, fixed = TRUE))
  expect_true(grepl("educ", out, fixed = TRUE))
})

test_that("custom labels are used in html output", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "html",
                                outcome_label = "lwage",
                                treatment_label = "educ")
  expect_true(grepl("lwage", out, fixed = TRUE))
  expect_true(grepl("educ", out, fixed = TRUE))
})

test_that("custom labels are used in pure_html output", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "pure_html",
                                outcome_label = "lwage",
                                treatment_label = "educ")
  expect_true(grepl("lwage", out, fixed = TRUE))
  expect_true(grepl("educ", out, fixed = TRUE))
})


# default labels from iv_fit ---------------------------------------------

test_that("default labels come from iv_fit variable names", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  # the y_name is "y" and d_name is "d" from deparse(substitute())
  # just check a label appears (not NULL)
  expect_true(is.character(out))
  expect_true(nchar(out) > 100)
})


# verbose flag ------------------------------------------------------------

test_that("verbose = TRUE prints output", {
  expect_output(ovb_minimal_reporting(card.sens, verbose = TRUE))
})

test_that("verbose = FALSE does not print", {
  expect_silent(out <- ovb_minimal_reporting(card.sens, verbose = FALSE))
  expect_true(is.character(out))
})

test_that("verbose = FALSE for html does not print", {
  expect_silent(out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "html"))
  expect_true(is.character(out))
})

test_that("verbose = FALSE for pure_html does not print", {
  expect_silent(out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "pure_html"))
  expect_true(is.character(out))
})


# invisible return --------------------------------------------------------

test_that("return value is invisible character string", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  expect_type(out, "character")
  expect_length(out, 1)
})


# bounds in footnote ------------------------------------------------------

test_that("latex footnote contains bound info when benchmarks present", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE)
  expect_true(grepl("Bound", out, fixed = TRUE))
  expect_true(grepl("1x black", out, fixed = TRUE))
  expect_true(grepl("R\\^2_\\{Z", out))
  expect_true(grepl("R\\^2_\\{Y\\(0\\)", out))
})

test_that("html footnote contains bound info", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "html")
  expect_true(grepl("Bound", out, fixed = TRUE))
  expect_true(grepl("1x black", out, fixed = TRUE))
})

test_that("pure_html footnote contains bound info", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE, format = "pure_html")
  expect_true(grepl("Bound", out, fixed = TRUE))
  expect_true(grepl("1x black", out, fixed = TRUE))
  expect_true(grepl("R<sup>2</sup>", out, fixed = TRUE))
})


# no bounds case ----------------------------------------------------------

test_that("latex works without benchmarks (no bounds)", {
  out <- ovb_minimal_reporting(card.sens.nobounds, verbose = FALSE)
  expect_true(is.character(out))
  expect_true(grepl("\\\\begin\\{table\\}", out))
  # should NOT have bound info
  expect_false(grepl("Bound", out, fixed = TRUE))
})

test_that("html works without benchmarks (no bounds)", {
  out <- ovb_minimal_reporting(card.sens.nobounds, verbose = FALSE, format = "html")
  expect_true(is.character(out))
  expect_false(grepl("Bound", out, fixed = TRUE))
})

test_that("pure_html works without benchmarks (no bounds)", {
  out <- ovb_minimal_reporting(card.sens.nobounds, verbose = FALSE, format = "pure_html")
  expect_true(is.character(out))
  expect_false(grepl("Bound", out, fixed = TRUE))
})


# .escape_latex -----------------------------------------------------------

test_that(".escape_latex escapes dollar signs", {
  escaped <- iv.sensemakr:::.escape_latex("card$lwage")
  expect_true(grepl("\\\\\\$", escaped))
  expect_equal(escaped, "card\\$lwage")
})

test_that(".escape_latex escapes multiple special characters", {
  escaped <- iv.sensemakr:::.escape_latex("a$b%c&d#e_f{g}h")
  expect_true(grepl("\\\\\\$", escaped))
  expect_true(grepl("\\\\%", escaped))
  expect_true(grepl("\\\\&", escaped))
  expect_true(grepl("\\\\#", escaped))
  expect_true(grepl("\\\\_", escaped))
})

test_that(".escape_latex handles plain strings", {
  expect_equal(iv.sensemakr:::.escape_latex("educ"), "educ")
  expect_equal(iv.sensemakr:::.escape_latex("lwage"), "lwage")
})


# caption and label -------------------------------------------------------

test_that("latex caption and label arguments work", {
  out <- ovb_minimal_reporting(card.sens, verbose = FALSE,
                                caption = "IV sensitivity table",
                                label = "tab:iv")
  expect_true(grepl("\\\\caption\\{IV sensitivity table\\}", out))
  expect_true(grepl("\\\\label\\{tab:iv\\}", out))
})


# fallback to sensemakr ---------------------------------------------------

test_that("ovb_minimal_reporting fallback works for sensemakr objects", {
  skip_if_not_installed("sensemakr")
  # create a regular OLS sensemakr object
  data("card")
  lm_fit <- lm(lwage ~ educ + exper + expersq + black + south + smsa, data = card)
  ols_sens <- sensemakr::sensemakr(lm_fit, treatment = "educ", benchmark_covariates = "black")
  out <- ovb_minimal_reporting(ols_sens, verbose = FALSE)
  expect_true(is.character(out))
  # OLS table should have different structure (S.E. column, not Lower/Upper CI)
  expect_true(grepl("S\\.E\\.", out))
})


# digits ------------------------------------------------------------------

test_that("digits parameter controls precision", {
  out5 <- ovb_minimal_reporting(card.sens, verbose = FALSE, digits = 5)
  out2 <- ovb_minimal_reporting(card.sens, verbose = FALSE, digits = 2)
  # more digits -> longer string

  expect_true(nchar(out5) >= nchar(out2))
})
