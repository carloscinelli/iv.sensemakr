# Tests that all documentation examples run without error.
# This ensures @examples in Rd files are valid and functional.

# --- Helper: set up Card data ---
card_setup <- function() {
  data("card", package = "iv.sensemakr")
  y <- card$lwage
  d <- card$educ
  z <- card$nearc4
  x <- model.matrix(
    ~ exper + expersq + black + south + smsa + reg661 + reg662 +
      reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
    data = card
  )
  list(y = y, d = d, z = z, x = x)
}


test_that("iv_fit example from sensemakr.Rd runs correctly", {
  dat <- card_setup()

  # Fits IV model
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  expect_s3_class(card.fit, "iv_fit")

  # See results
  expect_output(print(card.fit), "Instrumental Variable Estimation")

  # Runs sensitivity analysis
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))
  expect_s3_class(card.sens, "iv.sensemakr")

  # See results
  expect_output(print(card.sens), "Sensitivity Analysis")

  # Sensitivity contour plot
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(card.sens, lim = 0.09))
})


test_that("coef.iv_fit examples run correctly", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Extract coefficients
  expect_true(is.numeric(coef(card.fit)))
  expect_true(is.numeric(coef(card.fit, parm = "fs")))
  expect_true(is.numeric(coef(card.fit, parm = "rf")))

  # Extract confidence intervals
  ci.iv <- confint(card.fit)
  expect_true(is.numeric(ci.iv))
  expect_true(length(ci.iv) >= 2)

  ci.fs <- confint(card.fit, parm = "fs")
  expect_true(is.numeric(ci.fs))

  ci.rf <- confint(card.fit, parm = "rf")
  expect_true(is.numeric(ci.rf))
})


test_that("summary methods produce output", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # iv_fit summary
  expect_output(print(summary(card.fit)), "Instrumental Variable Estimation")

  # sensemakr summary
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))
  expect_output(print(summary(card.sens)), "Sensitivity Analysis")
})


test_that("robustness_value examples run correctly", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # RV and XRV
  rv.val <- rv(card.fit)
  expect_true(is.numeric(rv.val))
  expect_true(rv.val > 0 && rv.val < 1)

  xrv.val <- xrv(card.fit)
  expect_true(is.numeric(xrv.val))
  expect_true(xrv.val > 0 && xrv.val < 1)

  # Different parms
  rv.fs <- rv(card.fit, parm = "fs")
  expect_true(is.numeric(rv.fs))

  rv.rf <- rv(card.fit, parm = "rf")
  expect_true(is.numeric(rv.rf))
})


test_that("ovb_contour_plot runs without error for all parm types", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  # IV lower limit
  expect_no_error(
    ovb_contour_plot(card.fit, benchmark_covariates = c("black", "smsa"),
                     lim = 0.08, sensitivity.of = "lwr")
  )

  # IV upper limit
  expect_no_error(
    ovb_contour_plot(card.fit, benchmark_covariates = c("black", "smsa"),
                     lim = 0.08, sensitivity.of = "upr")
  )

  # FS t-value
  expect_no_error(
    ovb_contour_plot(card.fit, parm = "fs", benchmark_covariates = c("black", "smsa"),
                     lim = 0.08, sensitivity.of = "t-value")
  )

  # RF t-value
  expect_no_error(
    ovb_contour_plot(card.fit, parm = "rf", benchmark_covariates = c("black", "smsa"),
                     lim = 0.08, sensitivity.of = "t-value")
  )

  # With manual bounds
  expect_no_error(
    ovb_contour_plot(card.fit, r2zw.x = 0.02, lim = 0.08)
  )
})


test_that("sensitivity_stats runs correctly", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  ss.iv <- sensitivity_stats(card.fit, parm = "iv")
  expect_true(is.data.frame(ss.iv))
  expect_true("estimate" %in% names(ss.iv))
  expect_true("rv_qa" %in% names(ss.iv))

  ss.fs <- sensitivity_stats(card.fit, parm = "fs")
  expect_true(is.data.frame(ss.fs))

  ss.rf <- sensitivity_stats(card.fit, parm = "rf")
  expect_true(is.data.frame(ss.rf))
})


test_that("coef.iv.sensemakr extracts correct estimates", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

  est.iv <- coef(card.sens, parm = "iv")
  expect_true(is.data.frame(est.iv) || is.numeric(est.iv))

  est.all <- coef(card.sens, parm = c("iv", "fs", "rf"))
  expect_true(nrow(est.all) == 3 || length(est.all) == 3)
})


test_that("card dataset loads and has expected structure", {
  data("card", package = "iv.sensemakr")
  expect_true(is.data.frame(card))
  expect_equal(nrow(card), 3010)
  expect_true(all(c("lwage", "educ", "nearc4") %in% names(card)))
})
