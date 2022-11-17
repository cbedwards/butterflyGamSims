test_that("anchors work", {
  dat.test = data.frame(years = c(1991, 1992),
                        doy = c(130, 150),
                        count = c(5, 10))
  out = add_anchors(dat.test, anchor.dist = 10)
  expect_equal(sort(out$years), c(1991, 1991, 1991, 1992, 1992, 1992))
  expect_equal(sort(out$doy), c(120, 120, 130, 150, 160, 160))
  expect_equal(sort(out$count), c(0,0,0,0, 5, 10))
})
