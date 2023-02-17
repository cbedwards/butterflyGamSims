
#' Identify reasonable first and last sample days from activity curve
#'
#' Use parameterized activity curve to identify reasonable first and last days of sampling based on threshold
#' of activity (relative to maximum activity). This captures our expectation that sampling is
#' based on knowledgeable practitioners who time their sampling based on butterfly activity.
#'
#' @param height.threshold Identifies the minimum scaled activity for "reasonable sampling".
#'  This is based on the maximum of the activity curve (so independent of actual abundance). Reasonable values
#'  seem to be between 0.05 and 0.001, but use `make.plot` to explore.
#' @inheritParams activity_gen
#' @param make.plot If `TRUE`, provide a plot of the activity curve and the proposed first and last sampling days.
#'
#' @return named vector of the first and last day of year to be used
#' @export
#'
#' @examples
#' doy_selector(height.threshold = 0.004, activity.type = "gauss", act.mean = 100, act.sd = 14)
doy_selector = function(height.threshold, activity.type, make.plot = TRUE, ...){
  stopifnot(is.numeric(height.threshold) &
              length(height.threshold)==1 &
              height.threshold >= 0 &
              height.threshold < 1)
  doy.vec = seq(0, 365, by = .1)
  abund.vec = rep(1000, length(doy.vec))
  activity.curve = activity_gen(abund.vec = abund.vec, doy = doy.vec, activity.type = activity.type,...)
  activity.scale = activity.curve/max(activity.curve)
  above.thresh = activity.scale >= height.threshold
  doy.first = doy.vec[min(which(above.thresh))]
  doy.last = doy.vec[max(which(above.thresh))]
  if(make.plot){
    plot(x = doy.vec, y = activity.curve, type ='l')
    graphics::abline(v = c(doy.first, doy.last), col = 'blue', lty = 2)
  }
  return(c(doy.first = doy.first, doy.last = doy.last))
}
