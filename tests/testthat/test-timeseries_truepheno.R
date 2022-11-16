## note: because of numerical imperfections, we expect minor differneces. Using
## a tolerance of 0.1, or 2.4 hours.
test_that("Gaussian gives right answer", {
  out = timeseries_truepheno(
    activity.type = "gauss",
    act.mean = 130,
    act.sd = 15)
  expect_equal(out$median, 130, tolerance = 0.1)
  expect_equal(out$onset, qnorm(0.1, mean = 130, sd = 15), tolerance = 0.1)
  expect_equal(out$end, qnorm(0.9, mean = 130, sd = 15), tolerance = 0.1)
  expect_equal(out$fp,
               diff(c(qnorm(0.1, mean = 130, sd = 15),
                    qnorm(0.9, mean = 130, sd = 15))),
                    tolerance = 0.1)
})
