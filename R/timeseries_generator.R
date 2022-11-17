#' Generate time series
#'
#' @inheritParams timeseries_examples
#' @inherit timeseries_examples details
#' @return Data frame of simulated time series. Columns `years`, `doy`, and `count` contain the
#' the key data: the year and day of simulated observations, and the simulated count of
#' butterflies observed. `abund.true` gives the abundance index of the underlying activity curve
#' for that year, and `act` gives the value of the activity curve for that year and doy
#' (ie the expectation used to simulate `count` from.
#' `onset.true`, `median.true`, `end.true`, and `fp.true` give the phenology metrics calculated
#' on the underlying activity curve. The `*.true` columns can be used to compare
#' with estimates of these various metrics. Note that because of random sampling,
#' the count data will NOT support these estimates in any given instance,
#'  but repeated simulation and analysis using an unbiased method should produce
#'  metrics whose average converges on these `*.true` values. Note that when using
#'  `sample.type == "zinb"` option, `abund.true` is for the actual activity curve,
#'  but we expect our abundance estimates to be reduced by the zero inflation controlled by
#'  the "pzeros" argument.
#' @export
#' @examples
#' timeseries_generator(years = 1990:2000,
#'   doy.samples = seq(105,160, by = 7),
#'   abund.type = "exp",
#'   activity.type = "gauss",
#'   sample.type = "pois",
#'   growth.rate = -0.12,
#'   init.size = 500,
#'   act.mean = 130,
#'   act.sd = 15)
timeseries_generator = function(years,
                                doy.samples,
                                abund.type,
                                activity.type,
                                sample.type,
                                ...){
  dat = expand.grid(years = years,
                    doy = doy.samples)
  abund.merge = data.frame(years = years,
                           abund = abund_generator(unique(dat$years),
                                                   abund.type = abund.type,
                                                   ...))
  dat = dplyr::inner_join(dat, abund.merge, by = "years")
  dat$act = with(dat,
                 activity_gen(abund.vec = abund,
                              doy = doy,
                              activity.type = activity.type,
                              ...))
  dat$count = activity_sampler(dat$act,
                               sample.type,
                               ...)
  names(dat)[names(dat) == "abund"] = "abund.true"
  truepheno = timeseries_truepheno(activity.type,...)
  dat$onset.true = truepheno$onset
  dat$median.true = truepheno$median
  dat$end.true = truepheno$end
  dat$fp.true = truepheno$fp
  dat = dat[order(dat$years),]
  dat = dat[c("years", "doy", "count", "act", "abund.true","onset.true","median.true", "end.true","fp.true")]
  return(dat)
}
