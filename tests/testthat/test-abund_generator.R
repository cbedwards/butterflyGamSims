test_that("Argument checking works for abund_generator_rlnorm", {
  expect_error(abund_generator_rlnorm(year = 1:10,
                                      mean = 10,
                                      sd = 5))
  #what if we mix up arguments?
  expect_error(abund_generator_rlnorm(years = 1:10,
                                      growth.rate = 10,
                                      init.size = 5))
})

test_that("correct arguments go through, abund_generator_rlnorm",{
  expect_length(abund_generator_rlnorm(years = 1:10,
                                       meanlog = 10,
                                       sdlog = 5),
                10)

})

test_that("Argument checking works for abund_generator_exp", {
  expect_error(abund_generator_exp(years = 1:10,
                                   mean = 10,
                                   sd = 5))
  expect_error(abund_generator_exp(years = 1:10,
                                   meanlog = 10,
                                   sdlog = 5))
})

test_that("results correct for basic cases, abund_generator_exp",{
  ## constant population
  expect_equal(abund_generator_exp(years = 1:10,
                                    growth.rate = 0,
                                    init.size=50),
                rep(50, 10))
  #10 year decline, starting year number is random
  expect_equal(tail(abund_generator_exp(years = 0:10+runif(1)*2000,
                                        growth.rate = -.1,
                                        init.size=50),1),
               50*exp(-.1*10))

})
