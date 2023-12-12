local_edition(2)
context("Manual Bounds")

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
  ovb_contour_plot(card.fit, lim = 0.09)
  ovb_contour_plot(card.fit, r2zw.x = 0.02, lim = 0.09)
  capture.output(print(card.fit))
  capture.output(print(summary(card.fit)))
  capture.output(print(confint(card.fit)))

  card.sens <- sensemakr(card.fit, r2zw.x = 0.1)
  capture.output(print(card.sens))

  capture.output(print(summary(card.sens)))

})
