# Numerical validation tests for iv.sensemakr
# These tests verify computed quantities against hand-computed values.

# --- Helper: set up Card data consistently across tests ---
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


# ==========================================================================
# (a) IV estimation correctness
# ==========================================================================

test_that("IV point estimate equals RF/FS coefficient ratio from direct lm()", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Reconstruct the FS and RF regressions manually
  df <- card.fit$data
  fs <- lm(d ~ -1 + ., data = df[, c("d", "z", setdiff(names(df), c("y", "d")))])
  rf <- lm(y ~ -1 + ., data = df[, c("y", "z", setdiff(names(df), c("y", "d")))])

  # IV = RF_z / FS_z
  iv.manual <- unname(coef(rf)["z"] / coef(fs)["z"])
  expect_equal(coef(card.fit, parm = "iv")[["iv"]], iv.manual)
})

test_that("FS and RF coefficients match direct lm() output", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Extract FS estimates from the stored lm model
  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf

  fs.summ <- coef(summary(fs.lm))["z", ]
  rf.summ <- coef(summary(rf.lm))["z", ]

  # FS: coefficient, SE, t-value, p-value

  expect_equal(card.fit$estimates$fs$estimate, unname(fs.summ["Estimate"]))
  expect_equal(card.fit$estimates$fs$se,       unname(fs.summ["Std. Error"]))
  expect_equal(card.fit$estimates$fs$t.value,  unname(fs.summ["t value"]))
  expect_equal(card.fit$estimates$fs$p.value,  unname(fs.summ["Pr(>|t|)"]))

  # RF: coefficient, SE, t-value, p-value
  expect_equal(card.fit$estimates$rf$estimate, unname(rf.summ["Estimate"]))
  expect_equal(card.fit$estimates$rf$se,       unname(rf.summ["Std. Error"]))
  expect_equal(card.fit$estimates$rf$t.value,  unname(rf.summ["t value"]))
  expect_equal(card.fit$estimates$rf$p.value,  unname(rf.summ["Pr(>|t|)"]))
})

test_that("AR t-value matches lm(y - h0*d ~ z + x) under null", {
  dat <- card_setup()
  # Fit with h0 = 0 (default)
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  ar.lm <- card.fit$models$ar
  ar.summ <- coef(summary(ar.lm))["z", ]

  expect_equal(card.fit$estimates$iv$t.value, unname(ar.summ["t value"]))
  expect_equal(card.fit$estimates$iv$p.value, unname(ar.summ["Pr(>|t|)"]))
})

test_that("AR CI endpoints satisfy test-inversion property", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  ci <- confint(card.fit, parm = "iv")
  alpha <- card.fit$pars$alpha
  dof <- card.fit$estimates$iv$dof

  # Critical F value
  ts <- sqrt(qf(p = log(alpha), df1 = 1, df2 = dof, log.p = TRUE, lower.tail = FALSE))

  # At each CI endpoint tau0, the AR test stat should equal the critical value
  for (tau0 in ci) {
    if (is.finite(tau0)) {
      # Refit AR regression at this endpoint
      df <- card.fit$data
      df$ar_y <- df$y - tau0 * df$d
      covnames <- setdiff(names(df), c("y", "d", "ar_y"))
      fmla <- as.formula(paste0("ar_y ~ -1 + ", paste(paste0("`", covnames, "`"), collapse = " + ")))
      ar.lm <- lm(fmla, data = df)
      t.val <- abs(coef(summary(ar.lm))["z", "t value"])
      expect_equal(t.val, ts, tolerance = 1e-6)
    }
  }
})

test_that("confint for FS and RF match standard lm confint", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # FS
  ci.fs.pkg <- confint(card.fit, parm = "fs")
  ci.fs.lm  <- confint(card.fit$models$fs, parm = "z")
  expect_equal(unname(ci.fs.pkg), unname(ci.fs.lm))

  # RF
  ci.rf.pkg <- confint(card.fit, parm = "rf")
  ci.rf.lm  <- confint(card.fit$models$rf, parm = "z")
  expect_equal(unname(ci.rf.pkg), unname(ci.rf.lm))
})


# ==========================================================================
# (b) Robustness values
# ==========================================================================

test_that("RV for IV matches manual rv_from_summary_stats computation", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Extract summary stats manually
  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))
  dof <- fs.lm$df.residual

  # Manual RV computation (q=1)
  q <- 1
  tau.r <- rf.coef / fs.coef
  tau0 <- (1 - q) * tau.r  # = 0 when q=1
  phi.r <- rf.coef - tau0 * fs.coef
  se.phi.r <- sqrt(rf.se^2 + tau0^2 * fs.se^2 - 2 * tau0 * rho.val * fs.se * rf.se)
  t.phi.r <- abs(phi.r / se.phi.r)
  t.fs.r <- abs(fs.coef / fs.se)

  rv.phi <- sensemakr::rv(t_statistic = t.phi.r, dof = dof, q = 1, alpha = 0.05)
  rv.fs  <- sensemakr::rv(t_statistic = t.fs.r,  dof = dof, q = 1, alpha = 0.05)

  # RV_IV = min(rv.phi, rv.fs) when min=TRUE
  rv.iv.manual <- min(rv.phi, rv.fs)

  rv.iv.pkg <- as.numeric(rv(card.fit, parm = "iv", q = 1, alpha = 0.05))
  expect_equal(rv.iv.pkg, rv.iv.manual)
})

test_that("XRV for IV matches manual computation", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))
  dof <- fs.lm$df.residual

  tau.r <- rf.coef / fs.coef
  tau0 <- 0  # q=1
  phi.r <- rf.coef - tau0 * fs.coef
  se.phi.r <- sqrt(rf.se^2 + tau0^2 * fs.se^2 - 2 * tau0 * rho.val * fs.se * rf.se)
  t.phi.r <- abs(phi.r / se.phi.r)
  t.fs.r <- abs(fs.coef / fs.se)

  xrv.phi <- sensemakr::xrv(t_statistic = t.phi.r, dof = dof, q = 1, alpha = 0.05)
  xrv.fs  <- sensemakr::xrv(t_statistic = t.fs.r,  dof = dof, q = 1, alpha = 0.05)

  xrv.iv.manual <- min(xrv.phi, xrv.fs)
  xrv.iv.pkg <- as.numeric(xrv(card.fit, parm = "iv", q = 1, alpha = 0.05))
  expect_equal(xrv.iv.pkg, xrv.iv.manual)
})

test_that("RV with q=0.5 matches manual computation", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))
  dof <- fs.lm$df.residual

  q <- 0.5
  tau.r <- rf.coef / fs.coef
  tau0 <- (1 - q) * tau.r
  phi.r <- rf.coef - tau0 * fs.coef
  se.phi.r <- sqrt(rf.se^2 + tau0^2 * fs.se^2 - 2 * tau0 * rho.val * fs.se * rf.se)
  t.phi.r <- abs(phi.r / se.phi.r)
  t.fs.r <- abs(fs.coef / fs.se)

  rv.phi <- sensemakr::rv(t_statistic = t.phi.r, dof = dof, q = 1, alpha = 0.05)
  rv.fs  <- sensemakr::rv(t_statistic = t.fs.r,  dof = dof, q = 1, alpha = 0.05)
  rv.iv.manual <- min(rv.phi, rv.fs)

  rv.iv.pkg <- as.numeric(rv(card.fit, parm = "iv", q = q, alpha = 0.05))
  expect_equal(rv.iv.pkg, rv.iv.manual)
})

test_that("RV with min=FALSE returns rv.phi (not min)", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))
  dof <- fs.lm$df.residual

  q <- 1
  tau0 <- 0
  phi.r <- rf.coef - tau0 * fs.coef
  se.phi.r <- sqrt(rf.se^2 + tau0^2 * fs.se^2 - 2 * tau0 * rho.val * fs.se * rf.se)
  t.phi.r <- abs(phi.r / se.phi.r)

  rv.phi <- sensemakr::rv(t_statistic = t.phi.r, dof = dof, q = 1, alpha = 0.05)

  rv.iv.pkg <- as.numeric(rv(card.fit, parm = "iv", q = 1, alpha = 0.05, min = FALSE))
  expect_equal(rv.iv.pkg, as.numeric(rv.phi))
})

test_that("RV for FS and RF dispatch to sensemakr correctly", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # FS
  rv.fs.pkg <- as.numeric(rv(card.fit, parm = "fs"))
  rv.fs.direct <- as.numeric(sensemakr::rv(card.fit$models$fs, covariates = "z"))
  expect_equal(rv.fs.pkg, rv.fs.direct)

  # RF
  rv.rf.pkg <- as.numeric(rv(card.fit, parm = "rf"))
  rv.rf.direct <- as.numeric(sensemakr::rv(card.fit$models$rf, covariates = "z"))
  expect_equal(rv.rf.pkg, rv.rf.direct)
})


# ==========================================================================
# (c) Sensitivity analysis structure and bounds
# ==========================================================================

test_that("sensemakr() output has correct structure", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

  expect_s3_class(card.sens, "iv.sensemakr")
  expect_true(all(c("pars", "unadjusted", "sensitivity_stats", "bounds") %in% names(card.sens)))
  expect_s3_class(card.sens$unadjusted, "iv_fit")
  expect_true(all(c("iv", "fs", "rf") %in% names(card.sens$sensitivity_stats)))
  expect_true(all(c("iv", "fs", "rf") %in% names(card.sens$bounds)))
})

test_that("sensitivity_stats have expected columns", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

  expected_cols <- c("estimate", "lwr", "upr", "t.value", "xrv_qa", "rv_qa", "q", "min", "alpha", "dof")
  for (parm in c("iv", "fs", "rf")) {
    ss <- card.sens$sensitivity_stats[[parm]]
    expect_true(all(expected_cols %in% names(ss)),
                info = paste("Missing columns for", parm))
  }
})

test_that("bounds partial R2 values are in [0, 1]", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

  iv_bounds <- card.sens$bounds$iv
  expect_true(all(iv_bounds$r2zw.x >= 0 & iv_bounds$r2zw.x <= 1))
  expect_true(all(iv_bounds$r2y0w.zx >= 0 & iv_bounds$r2y0w.zx <= 1))
})

test_that("manual bounds (r2zw.x) produce valid adjusted CI", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, r2zw.x = 0.02)

  # Manual bounds should be in the bounds table
  expect_true(nrow(card.sens$bounds$iv) >= 1)
  manual_row <- card.sens$bounds$iv[card.sens$bounds$iv$bound_label == "Manual Bound", ]
  expect_equal(nrow(manual_row), 1)
  expect_equal(manual_row$r2zw.x, 0.02)
  # Adjusted CI should be wider than unadjusted (lower bound smaller, upper bound larger)
  unadj_ci <- confint(card.fit, parm = "iv")
  expect_true(manual_row$lwr <= unadj_ci[1] || manual_row$upr >= unadj_ci[2])
})


# ==========================================================================
# (d) AR CI edge cases
# ==========================================================================

test_that("irrelevant instrument produces unbounded CI", {
  set.seed(123)
  n <- 100
  z <- rnorm(n)
  d <- rnorm(n)  # instrument completely irrelevant
  y <- d + rnorm(n)
  fit <- iv_fit(y, d, z)
  ci <- confint(fit, parm = "iv")
  # Should be unbounded (-Inf, Inf)
  expect_true(any(is.infinite(ci)))
})

test_that("iv_fit works without covariates (x=NULL)", {
  set.seed(42)
  n <- 500
  z <- rbinom(n, 1, 0.5)
  d <- 0.5 * z + rnorm(n)
  y <- 2 * d + rnorm(n)
  fit <- iv_fit(y, d, z)
  expect_s3_class(fit, "iv_fit")
  expect_true(is.numeric(coef(fit, parm = "iv")))
  ci <- confint(fit, parm = "iv")
  expect_true(length(ci) >= 2)
})

test_that("confint with different levels works", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  ci_95 <- confint(card.fit, parm = "iv")
  ci_99 <- confint(card.fit, parm = "iv", level = 0.99)

  # 99% CI should be wider than 95%
  expect_true(ci_99[1] <= ci_95[1])
  expect_true(ci_99[2] >= ci_95[2])
})
