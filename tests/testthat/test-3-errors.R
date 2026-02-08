test_that("iv_fit validates inputs", {
  expect_error(iv_fit("a", 1, 1), "numeric")
  expect_error(iv_fit(1, "a", 1), "numeric")
  expect_error(iv_fit(1, 1, "a"), "numeric")

  expect_error(iv_fit(1:2, 1:3, 1:2), "same sample size")
})

test_that("alpha checks are enforced in robustness_value methods", {
  data("card", package = "iv.sensemakr")
  x <- model.matrix(
    ~ exper + expersq + black + south + smsa + reg661 + reg662 +
      reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
    data = card
  )
  card.fit <- iv_fit(card$lwage, card$educ, card$nearc4, x)

  expect_error(robustness_value(card.fit, alpha = -0.1), "alpha")
  expect_error(extreme_robustness_value(card.fit, alpha = 2), "alpha")
})

test_that("ovb_contour_plot checks r2 bounds", {
  data("card", package = "iv.sensemakr")
  x <- model.matrix(
    ~ exper + expersq + black + south + smsa + reg661 + reg662 +
      reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
    data = card
  )
  card.fit <- iv_fit(card$lwage, card$educ, card$nearc4, x)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(ovb_contour_plot(card.fit, r2zw.x = -0.01, lim = 0.05))
  expect_error(ovb_contour_plot(card.fit, r2zw.x = 2, lim = 0.05))
})

