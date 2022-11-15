
## check answers for uniform counts
## We expect the quantiles to be slightly "off" from the quantiles of a continuous distribution.
## We'll just tweak the tolerance very slightly to accomodate.
test_that("Metrics correct for uniform activity", {
  N = 100000
  out = gam_summarizer(count.pred = rep(1, N),
                       doy.pred = seq(0,365, length = N))
  expect_equal(out$abund, 365, tolerance = 0.001)
  expect_equal(out$med, 365/2, tolerance = 0.001)
  expect_equal(out$onset, 365/10, tolerance = 0.001)
  expect_equal(out$end, 365 * 9/10, tolerance = 0.001)
  expect_equal(out$fp, 365 * 8/10, tolerance = 0.001)
  expect_equal(out$boundary.reasonable.rel, FALSE)
  expect_equal(out$boundary.reasonable.abs, FALSE)
})

test_that("absolute goodfit works", {
  N=100
  out = gam_summarizer(count.pred = rep(.001, N),
                       doy.pred = seq(0,365, length = N))
  expect_equal(out$boundary.reasonable.rel, FALSE)
  expect_equal(out$boundary.reasonable.abs, TRUE)
}
)

test_that("Metrics correct for gaussian activity", {
  N = 100000
  activity.mean = 100; activity.sd =10
  out = gam_summarizer(count.pred = 100 * dnorm(seq(0,365, length = N),
                                          mean = activity.mean , sd = activity.sd),
                       doy.pred = seq(0,365, length = N))
  expect_equal(out$abund, 100, tolerance = 0.001)
  expect_equal(out$med, activity.mean, tolerance = 0.001)
  expect_equal(out$onset,
               qnorm(.1, mean = activity.mean, sd = activity.sd),
               tolerance = 0.001)
  expect_equal(out$end,
               qnorm(.9, mean = activity.mean, sd = activity.sd),
               tolerance = 0.001)
  expect_equal(out$fp,
               diff(qnorm(c(0.1, 0.9),mean = activity.mean, sd = activity.sd)),
               tolerance = 0.001)
  expect_equal(out$boundary.reasonable.rel, TRUE)
  expect_equal(out$boundary.reasonable.abs, TRUE)
})

test_that("errors correctly", {
  expect_error(gam_summarizer(1:10, "hello"))
  expect_error(gam_summarizer(TRUE, 1:10))
  expect_error(gam_summarizer(1:9, 1:10))
  expect_error(gam_summarizer(1:10, 1:10, bounds.reasonable = "hello"))
  expect_error(gam_summarizer(1:10, 1:10, bounds.thresh.rel = 1/(1:2)))
  expect_error(gam_summarizer(1:10, 1:10, bounds.thresh.abs = 1/(1:2)))
})

