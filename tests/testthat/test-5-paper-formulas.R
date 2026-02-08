# Tests verifying package computations against the paper's formulas.
# Cinelli and Hazlett (2025), "An Omitted Variable Bias Framework for
# Sensitivity Analysis of Instrumental Variables", Biometrika.

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


# ==========================================================================
# (a) IV estimate: tau = beta_RF / beta_FS
# ==========================================================================

test_that("IV estimate equals ratio of independently fitted RF and FS models", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  df <- card.fit$data

  # Fit FS and RF independently
  covs <- setdiff(names(df), c("y", "d"))
  fmla.fs <- as.formula(paste0("d ~ -1 + ", paste(paste0("`", covs, "`"), collapse = " + ")))
  fmla.rf <- as.formula(paste0("y ~ -1 + ", paste(paste0("`", covs, "`"), collapse = " + ")))

  fs <- lm(fmla.fs, data = df)
  rf <- lm(fmla.rf, data = df)

  tau.manual <- unname(coef(rf)["z"] / coef(fs)["z"])
  tau.pkg    <- coef(card.fit, parm = "iv")[["iv"]]

  expect_equal(tau.pkg, tau.manual)
})


# ==========================================================================
# (b) AR quadratic: manually solve a*tau^2 + b*tau + c <= 0
# ==========================================================================

test_that("AR CI matches manual quadratic solution", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Extract summary stats
  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))
  dof <- fs.lm$df.residual
  alpha <- 0.05

  # Critical value
  ts <- sqrt(qf(p = log(alpha), df1 = 1, df2 = dof, log.p = TRUE, lower.tail = FALSE))

  # Quadratic coefficients
  a <- fs.coef^2 - fs.se^2 * ts^2
  b <- 2 * rho.val * rf.se * fs.se * ts^2 - 2 * rf.coef * fs.coef
  cc <- rf.coef^2 - rf.se^2 * ts^2
  delta <- b^2 - 4 * a * cc

  # Solve
  expect_true(a > 0, info = "Card data should have strong enough instrument for a > 0")
  expect_true(delta >= 0, info = "Card data should have positive discriminant")

  roots <- sort(c((-b + sqrt(delta)) / (2 * a), (-b - sqrt(delta)) / (2 * a)))

  # Compare to package
  ci.pkg <- confint(card.fit, parm = "iv")
  expect_equal(as.numeric(ci.pkg), roots, tolerance = 1e-10)
})


# ==========================================================================
# (c) RV decomposition: RV_IV = min(RV_phi, RV_FS)
# ==========================================================================

test_that("RV_IV is the minimum of RV_phi and RV_FS", {
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

  # phi = rf.coef - tau0 * fs.coef, with tau0 = 0 for q=1
  phi.r <- rf.coef
  se.phi.r <- rf.se  # simplifies when tau0=0
  t.phi.r <- abs(phi.r / se.phi.r)
  t.fs.r <- abs(fs.coef / fs.se)

  rv.phi <- as.numeric(sensemakr::rv(t_statistic = t.phi.r, dof = dof, q = 1, alpha = 0.05))
  rv.fs  <- as.numeric(sensemakr::rv(t_statistic = t.fs.r,  dof = dof, q = 1, alpha = 0.05))

  rv.iv <- as.numeric(rv(card.fit, q = 1, alpha = 0.05))

  # RV_IV = min(RV_phi, RV_FS)
  expect_equal(rv.iv, min(rv.phi, rv.fs))

  # In the Card example, the phi RV is the binding constraint
  expect_equal(rv.iv, rv.phi)
})

test_that("SE(phi) formula is correct for non-zero tau0", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  fs.lm <- card.fit$models$fs
  rf.lm <- card.fit$models$rf
  fs.coef <- unname(coef(fs.lm)["z"])
  rf.coef <- unname(coef(rf.lm)["z"])
  fs.se <- unname(coef(summary(fs.lm))["z", "Std. Error"])
  rf.se <- unname(coef(summary(rf.lm))["z", "Std. Error"])
  rho.val <- cor(resid(fs.lm), resid(rf.lm))

  # For q = 0.5, tau0 = 0.5 * tau.r
  q <- 0.5
  tau.r <- rf.coef / fs.coef
  tau0 <- (1 - q) * tau.r

  # Paper's SE formula
  se.phi.r <- sqrt(rf.se^2 + tau0^2 * fs.se^2 - 2 * tau0 * rho.val * fs.se * rf.se)

  # Verify by computing the AR regression manually at tau0
  df <- card.fit$data
  df$ar_y <- df$y - tau0 * df$d
  covnames <- setdiff(names(df), c("y", "d", "ar_y"))
  fmla <- as.formula(paste0("ar_y ~ -1 + ", paste(paste0("`", covnames, "`"), collapse = " + ")))
  ar.lm <- lm(fmla, data = df)
  ar.se <- unname(coef(summary(ar.lm))["z", "Std. Error"])

  expect_equal(se.phi.r, ar.se, tolerance = 1e-10)
})


# ==========================================================================
# (d) Verify Card results match paper's reported values
# ==========================================================================

test_that("Card IV estimate matches paper", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Exact computed value
  expect_equal(coef(card.fit, parm = "iv")[["iv"]], 0.1315038, tolerance = 1e-4)
})

test_that("Card t-value and p-value match paper", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  # Exact computed values
  expect_equal(card.fit$estimates$iv$t.value, 2.327075, tolerance = 1e-3)
  expect_equal(card.fit$estimates$iv$p.value, 0.02002763, tolerance = 1e-4)
})

test_that("Card 95% CI matches paper", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  ci <- confint(card.fit, parm = "iv")

  # Exact computed values
  expect_equal(ci[1], 0.02480484, tolerance = 1e-4)
  expect_equal(ci[2], 0.2848236,  tolerance = 1e-4)
})

test_that("Card robustness values match paper", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)

  rv.val  <- as.numeric(rv(card.fit))
  xrv.val <- as.numeric(xrv(card.fit))

  # Exact computed values
  expect_equal(xrv.val, 0.0005232443, tolerance = 1e-5)
  expect_equal(rv.val,  0.006666407,  tolerance = 1e-5)
})

test_that("Card benchmark bounds match paper", {
  dat <- card_setup()
  card.fit <- iv_fit(dat$y, dat$d, dat$z, dat$x)
  card.sens <- sensemakr(card.fit, benchmark_covariates = c("black", "smsa"))

  bounds <- card.sens$bounds$iv

  # black benchmark
  black_row <- bounds[grepl("black", bounds$bound_label), ]
  expect_true(nrow(black_row) == 1)
  expect_true(black_row$r2zw.x > 0)
  expect_true(black_row$r2y0w.zx > 0)

  # smsa benchmark
  smsa_row <- bounds[grepl("smsa", bounds$bound_label), ]
  expect_true(nrow(smsa_row) == 1)
  expect_true(smsa_row$r2zw.x > 0)
  expect_true(smsa_row$r2y0w.zx > 0)

  # Adjusted CIs should be wider than unadjusted
  unadj_ci <- confint(card.fit, parm = "iv")
  expect_true(black_row$lwr < unadj_ci[1] || black_row$upr > unadj_ci[2])
  expect_true(smsa_row$lwr < unadj_ci[1] || smsa_row$upr > unadj_ci[2])
})
