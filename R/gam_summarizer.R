#' Summarizing activity curve of single year
#'
#' @param count.pred predicted count data
#' @param doy.pred day of year corresponding to count.pred
#' @param bounds.reasonable For identifying unreasonable fits, provide the days at which we really SHOULD see ~ 0 abundance.
#' Defaults to 1 and 364.
#' @param bounds.thresh.rel When the activity should be 0, what is "close enough", as a proportion of maximum count?
#' @param bounds.thresh.abs When the activity should be 0, what is "close enough" in absolute terms (for cases with flat activity)
#'
#' @details In addition to estimating our key metrics (abundance, quantile-based phenology landmarks, flight period),
#' this function also makes a simple pass at identifying problem fits.
#' It does this by taking days when we expect abundance to be approximately zero,
#' and testing whether the nearest doy.pred to those days have abundance that is "low enough". Currently
#' "low enough" is reasonable.thresh.rel x maximum abundance. We might also consider adding
#' in an absolute threshold, since "zero" years might be expected to be completely flat.
#'
#' Note that this was intended to map to the last and first day we expect zero abundance, but
#' a vector with additional days can be added if that is desired
#' (ie for bivoltine species with a known period of no abundance between peaks).
#'
#' @return Returns frame with summary information. `abund` = abundance index, `median` = median,
#' `onset` = 0.1 quantile, `end` = 0.9 quantile, `fp` = flight period (days between 0.9 and 0.1 quantile).
#' `boundary.reasonable.rel` and `boundary.reasonable.abs` will be true if all the provided `bounds.reasonable` values
#' are "sufficiently small" as defined by `bounds.thresh.rel` and `bounds.thresh.abs`, respectively. See "Details" for more... details.
#'
#' As a reminder, abundance index is NOT abundance, but is expected to be proportional to it,
#'  with that proportionality constant depending on detectability and lifespan. If those are
#'  generally constant across years, changes in log abundance index across years will be *equal* to changes
#'  in log abundance.
#'
#' @export
#'
#' @examples
#'
#' #imagining our fitted curve is Gaussian:
#' N = 100000
#' activity.mean = 100; activity.sd =10
#' out = gam_summarizer(count.pred = 100* dnorm(seq(0,365, length = N),
#'                                             mean = activity.mean , sd = activity.sd),
#'                     doy.pred = seq(0,365, length = N))
#'
gam_summarizer=function(count.pred,
                        doy.pred,
                        bounds.reasonable = c(1,364),
                        bounds.thresh.rel = 1/10,
                        bounds.thresh.abs = 1/100){
  stopifnot(is.numeric(count.pred),
            is.numeric(doy.pred),
            is.numeric(bounds.reasonable),
            is.numeric(bounds.thresh.rel),
            is.numeric(bounds.thresh.abs),
            length(count.pred)==length(doy.pred),
            length(bounds.thresh.rel) == 1,
            length(bounds.thresh.abs) == 1)
  #function to calculate whatever abundance and phenology metrics we desire
  #  from a single year of predicted gam
  #  count.pred =
  #  doy.pred =

  # calculate timestep h (we assume the doy are evenly spaced)
  h = diff(doy.pred)[1]

  res.cur = data.frame(abund = sum(count.pred)*h)
  # abundance metric as area under curve
  ## NOTE THAT THIS IS IN ACTIVITY DAYS!!
  ## for example, if an adult has a lifespan of 15 days on average, expected
  ## number of individuals would be abund/15.


  # first calulate the cumulative distribution function
  cdf = cumsum(count.pred)/sum(count.pred)
  res.cur$onset = doy.pred[min(which(cdf>=0.1))]
  res.cur$median = doy.pred[min(which(cdf>=0.5))]
  res.cur$end = doy.pred[min(which(cdf>=0.9))]
  res.cur$fp = res.cur$end - res.cur$onset

  ## testing for bad fits
  bounds.inds = numeric(length(bounds.reasonable))
  for(i.bound in 1:length(bounds.reasonable)){
    bounds.inds[i.bound] = which.min(abs(doy.pred - bounds.reasonable[i.bound]))
  }
  res.cur$boundary.reasonable.rel = all(count.pred[bounds.inds] < bounds.thresh.rel * max(count.pred))
  res.cur$boundary.reasonable.abs = all(count.pred[bounds.inds] < bounds.thresh.abs)

  return(res.cur)

}
