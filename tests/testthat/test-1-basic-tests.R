local_edition(2)
context("Card Example")

test_that("Card Example", {
  expect_true(TRUE)

  rm(list = ls())
  data("card")
  y <- card$lwage
  d <- card$educ
  z <- card$nearc4
  x <- model.matrix( ~ exper + expersq + black + south + smsa + reg661 + reg662 +
                       reg663 + reg664 + reg665+ reg666 + reg667 + reg668 + smsa66,
                     data = card)
  card.fit <- iv_fit(y,d, z, x)
  card.fit

  summary(card.fit)

  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black","smsa"))
  card.sens

  summary(card.sens)
  coef(card.sens, parm = c("iv", "fs", "rf"))
  coef(card.fit, parm = c("fs", "rf", "iv"))
  confint(card.fit, parm = "fs")

  coef(card.sens)

  summary(card.sens)
  plot(card.sens, lim = 0.09, alpha = 0.5)
  plot(card.sens, parm = "fs",
       sensitivity.of = "t-value",
       lim = 0.08)
  plot(card.sens, parm = "rf",
       sensitivity.of = "t-value",
       lim = 0.08)

  ovb_contour_plot(card.fit,  benchmark_covariates = c("black", "smsa"), lim = 0.08)

  do.call("rbind", card.sens$sensitivity_stats)

  # coef
  coef(card.fit, parm = "iv")
  coef(card.fit)
  confint(card.fit)
  confint(card.fit, parm = "iv", level = 0.9998)
  confint(card.fit, parm = "fs", level = 0.9998)

  rv(card.fit)
  robustness_value(card.fit)
  rv(card.fit, q = 1, parm = "rf")
  rv(card.fit, q = 1, parm = "fs")
  xrv(card.fit)
  xrv(card.fit, q = 1, parm = "rf")
  xrv(card.fit, q = 1, parm = "fs")
  # summary
  card.summ <- summary(card.fit)
  card.summ


  # par(mfrow = c(2, 2))
  ovb_contour_plot(card.fit, lim = 0.08, sensitivity.of = "lwr",
                   xlab = "mylab", y = "mylab",
                   alpha = 0.05,
                   benchmark_covariates = c("black", "smsa"))

  ovb_contour_plot(card.fit, lim = 0.08, sensitivity.of = "lwr",
                   alpha = 1,
                   benchmark_covariates = c("black", "smsa"))

  ovb_contour_plot(card.fit, lim = 0.08, sensitivity.of = "upr",
                   benchmark_covariates = c("black", "smsa"))

  plot(sensemakr(card.fit$models$fs, treatment= "z"), sensitivity.of = "t-value", lim = 0.08)

  # ovb_contour_plot(card.fit$models$fs, treatment= "z", sensitivity.of = "t-value", lim = 0.08)


  ovb_contour_plot(card.fit, parm = "fs", kz = 2, lim = 0.08, xlab = "mylab",
                   sensitivity.of = "upr", alpha = 0.05,
                   benchmark_covariates = c("black", "smsa"))

  ovb_contour_plot(card.fit$models$fs,
                   treatment = "z",
                   kd = 2, lim = 0.08,
                   sensitivity.of = "upr",
                   alpha = 0.05,
                   benchmark_covariates = c("black", "smsa"))


  ovb_contour_plot(card.fit, parm = "rf", kz = 1, lim = 0.08,
                   sensitivity.of = "t-value",
                   benchmark_covariates = c("black", "smsa"))

  ovb_contour_plot(card.fit, parm = "fs", kz = 2, lim = 0.08, lim.y = 0.08, benchmark = c("black", "smsa"), sensitivity.of = "lwr")
  ovb_contour_plot(card.fit, parm = "rf", lim = 0.08, lim.y = 0.08, benchmark = c("black", "smsa"))
  ovb_contour_plot(card.fit$models$rf,"z", lim = 0.08, benchmark = c("black", "smsa"))
  sensemakr::ovb_contour_plot(card.fit$models$fs, treatment= "z", benchmark = "black", sensitivity.of = "t-value", lim = 0.08)
})
