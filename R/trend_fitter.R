#' Fit the abundance and phenology trends for summary statistics of time series
#'
#'
#' @param dat.filtered Data frame with phenology and abundance metrics for any number of years.
#' Generally just feed in `$summary` component of the results of `gam_fitter`, with or without
#' applying censorship filters.
#'
#' @return 1-row data frame with linear trends across years for each metric. For abundance,
#' `growth.rate` is the trend in log abundance across years. `$nyears` is the number of unique
#' years provided.
#' @export
#'
#' @examples
#' dat.sim = timeseries_generator(years = 1990:2000,
#' doy.samples = seq(105,160, by = 7),
#' abund.type = "exp",
#' activity.type = "gauss",
#' sample.type = "pois",
#' growth.rate = -0.12,
#' init.size = 500,
#' act.mean = 130,
#' act.sd = 15)
#' out = gam_fitter(years.vec = dat.sim$years,
#'                  doy.vec = dat.sim$doy,
#'                  count.vec = dat.sim$count,
#'                  doy.smooth = "cr",
#'                  doy.knots = 5,
#'                  years.smooth = "cr",
#'                  anchor.flag = FALSE
#' )
#' fit_plotter(dat.timeseries = dat.sim,
#'             dat.fitted = out$dat.fitted,
#'             activity.curve = out$activity.curve)
#' trend_fitter(out$summary)
trend_fitter = function(dat.filtered){
  res.cur = data.frame(growth.rate = -999,
                       median = -999,
                       onset = -999,
                       end = -999,
                       fp = -999)
  res.cur$nyears = length(unique(dat.filtered$years))

  ## abundance trend (log)
  dat.filtered$logabund=log(dat.filtered$abund)
  out = stats::lm(logabund ~ years, data=dat.filtered)
  res.cur$growth.rate = stats::coef(out)[2]

  # trend in median date
  out = stats::lm(median ~ years, data=dat.filtered)
  res.cur$median = stats::coef(out)[2]

  #onset
  out = stats::lm(onset ~ years, data=dat.filtered)
  res.cur$onset = stats::coef(summary(out))[2,2]

  #end
  out = stats::lm(end ~ years, data=dat.filtered)
  res.cur$end = stats::coef(summary(out))[2,2]

  #fp
  out = stats::lm(fp ~ years, data=dat.filtered)
  res.cur$fp = stats::coef(summary(out))[2,2]

  return(res.cur)
}

#' Calculate trends (and diagnostic metrics) for a given censorship method
#'
#' @param dat.summary Data frame with phenology and abundance metrics for any number of years.
#' Generally just feed in `$summary` component of the results of `gam_fitter`.
#' @param nobs.min Minimum number of total observations to include a year (REAL observations;
#' anchors are not counted). This is only relevant if years vary in the number of sampling days
#' @param nnzero.min Minimum number of non-zero observations to include a year
#' @param bound.reasonable.rel Should we only use years with "good" fits as identified
#' using bount.reasonable.rel? See `gam_summarizer` for details.
#' @param bound.reasonable.abs Should we only use years with "good" fits as identified
#' using bount.reasonable.abs? See `gam_summarizer` for details.
#'
#' @return 1-row data frame with estimated trends (`growth.rate`, `median`, `onset`,
#' `end`, `fp`), the number of years used to estimate the trends (`nyears`),
#' the filtering criterion used to censor years before fitting trends (arguments for
#' this function), and the original number of years before censoring (`nyear.original`))
#' @export
#'
#' @examples
#' set.seed(10)
#' dat.sim = timeseries_generator(years = 1990:2000,
#'                                doy.samples = seq(105,160, by = 10),
#'                                abund.type = "exp",
#'                                activity.type = "gauss",
#'                                sample.type = "pois",
#'                                growth.rate = -0.12,
#'                                init.size = 100,
#'                                act.mean = 130,
#'                                act.sd = 15)
#' out = gam_fitter(years.vec = dat.sim$years,
#'                  doy.vec = dat.sim$doy,
#'                  count.vec = dat.sim$count,
#'                  doy.smooth = "cr",
#'                  doy.knots = 5,
#'                  years.smooth = "cr",
#'                  anchor.flag = TRUE,
#'                  anchor.dist = 10
#' )
#' fit_plotter(dat.timeseries = dat.sim,
#'             dat.fitted = out$dat.fitted,
#'             activity.curve = out$activity.curve,
#'             xlim = 20)
#' trend_fitter(out$summary)
#' trend_method(out$summary,
#'              nobs.min = 0,
#'              nnzero.min = 3)

trend_method = function(dat.summary,
                        nobs.min,
                        nnzero.min,
                        bound.reasonable.rel = F,
                        bound.reasonable.abs = F){
  dat.use = dat.summary[dat.summary$n - dat.summary$nzero>=nnzero.min &
                          dat.summary$n >= nobs.min &
                          (!bound.reasonable.rel |
                             dat.summary$boundary.reasonable.rel) &
                          (!bound.reasonable.abs |
                             dat.summary$boundary.reasonable.abs), ]
  trend_fitter(dat.use)
  dat.res = cbind(trend_fitter(dat.use),
                  data.frame(nobs.min = nobs.min,
                             nnzero.min = nnzero.min,
                             bound.reasonable.rel = bound.reasonable.rel,
                             bound.reasonable.abs = bound.reasonable.abs,
                             nyear.original = nrow(dat.summary)))
  return(dat.res)
}
