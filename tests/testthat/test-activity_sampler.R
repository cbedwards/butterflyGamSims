test_that("pois behaves well", {
  expect_error(activity_samp_pos(1:10, theta = 5))
  expect_equal(mean(activity_samp_pois(rep(10, 100000))),
               10,
               tolerance = .01)

})

test_that("nb behaves well", {
  expect_error(activity_samp_nb(1:10))
  expect_error(activity_samp_nb(1:10, pzero = .1))
  expect_equal(mean(activity_samp_nb(rep(10, 100000), theta = 10E10)),
               10,
               tolerance = .01)

})


test_that("zinb behaves well", {
  expect_error(activity_samp_zinb(1:10))
  expect_error(activity_samp_zinb(1:10, theta = 5))
  expect_error(activity_samp_zinb(1:10, theta = 5, pzero = 5))
  expect_equal(sum(activity_samp_zinb(rep(10, 1000), theta = 5, pzero = 1)),
               0)
  expect_equal(mean(activity_samp_zinb(rep(10, 100000),
                                       theta = 10E10,
                                       pzero = 0)),
                   10,
               tolerance = .01)
})



