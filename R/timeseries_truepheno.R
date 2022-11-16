#' Calculate "true" phenology metrics for a given parameter combination
#'
#' We are simulating distributions that are not changing in phenology across years,
#' so we only need to estimate the metrics from a single year. This function uses
#' numeric calculations rather than analytical ones (a) for generalizability to
#' distributions without closed form solutions for quantiles, and (b) to be exactly comparable
#' to GAM methods. These phenology metrics are showing what we would calculate if
#' our GAMs perfectly fit the activity curve.
#'
#' @inheritParams activity_gen
#' @inherit activity_gen details
#'
#' @return 1-row data frame of phenology metrics. `onset`, `median`, and `end` give
#' the days of the 0.1, 0.5, and 0.9 quantile respectively, `fp` gives the "flight period",
#' the number of days between onset and end.
#' @export
#'
#' @examples
#' out = timeseries_truepheno(
#'   activity.type = "gauss",
#'   act.mean = 130,
#'   act.sd = 15)
#'
timeseries_truepheno = function(activity.type,
                                ...){
  dat = data.frame(doy = seq(0,365, by = .1))
  dat$abund = 1000
  dat$act = with(dat,
                 activity_gen(abund.vec = abund,
                              doy = doy,
                              activity.type = activity.type,
                              ...))
  metrics.gam = gam_summarizer(count.pred = dat$act,
                               doy.pred = dat$doy)
  return(metrics.gam[,c("onset","median","end","fp")])
}


