#' Fitting timeseries data with different gam fitting options
#'
#' @param years.vec Vector of the years for analysis, typically the `$years` column of the results of `timeseries_generator`.
#' @param doy.vec Vector of the doy for analysis, typically the `$doy` column of the results of `timeseries_generator`.
#' @param count.vec Vector of the counts for analysis, typically the `$count` column of the results of `timeseries_generator`.
#' @param doy.smooth What smoothing approach to use for the doy direction. Reasonable options to try include a cyclic spline (`"cc"`),
#' and cubic regression spline (`"cr"`). While thin plate regressions are powerful in unvariate cases,
#' for tensor product smooths as used here, they offer no advantage over cubic spliens (per `?mgcv::gam`, after `te`, `ti`, and `t2` are defined)
#' when choosing cubic regression splines, `gam_fitter` will add two bonus knots on day 0.5 and 364.5 so that the function wraps appropriately.
#' Shrinkage cubic regression splines (`"cs"`) are another possible approach, but likely to be suboptimal.
#' The shrinkage will tend to draw regions of low data density towards zero. This is good for
#' preventing unreasonable behavior outside the bounds of the data, but will tend to bias
#' estimated abundance downward, especially in years of low data.
#' @param doy.knots How many knots to use in the doy dimension? Not used (or not used normally?) for thin plate regression.
#' @param years.smooth What smoothing approach to use for the year direction. Reasonable options include cubic regression (`"cr"`), and random effect (`"re"`).
#' For cubic regression, `gam_fitter` uses a knot for
#' each year. Because of the shrinkage associated with random effects, using `year.smooth == "re"` is likely
#' to underestimate negative or positive trends in abundance across years (but this can be tested!).
#' @param anchor.flag Logical. Should add "anchor" zeroes (additional dummy observations with 0 counts) be added to each
#' year of data to help prevent unreasonable behavior outside of observed range. If TRUE, two additional observations are added
#' to each year, `anchor.dist` days before the first doy of observation in the data set, and
#' `anchor.dist` days after the last doy of observation in the data set.
#' some number of days before and after
#' @param anchor.dist How many days out from the real data should anchor zeroes be added? Defaults to NULL,
#' must be an integer if `anchor.dist==TRUE`.
#' @param limit.to.anchor For calculating metrics of the fitted gams, do we want to limit our calculations to only the time period between the anchors?
#' If `TRUE`, and `anchor.flag = FALSE`,  calculations will use only the period of time between the first and last observed data points.
#' @param ... additional arguments for identifying reasonable model fits. See `?gam_summarizer`
#'
#' @details `gam_fitter` fits the data using a tensor product smooth, with potentially different basis functions
#' for the doy and year dimensions. Fitting uses restricted maximum likelihood (`methods = "REML"` in `gam()`) and
#' a negative binomial distribution (`family = "nb"` in `gam()`).
#'
#' Note that in general we have found that using anchor zeroes improves model fit and reduces bias. Additionally, we expect that the "limit.to.anchor" should be TRUE.
#' Most butterflies are active for brief periods in the year, and spline predictions outside of those ranges are biologically meaningless, but given how
#' splines can work (especially 2-dimensional splines with changing population abundances), these extrapolations can qualitatively change predictions. We leave the limit.to.anchor as an argument (and have it default to FALSE) for the purposes of unit tests.
#'
#' @return List with four components. `$summary` is a data frame frame with estimated attributes for
#' each fitted year. See `gam_summarizer()` for details on most columns. Two additional columns are added.
#' `$nzeros` is the number of 0 counts for each year, EXCLUDING ANCHORS (if relevant). `$n` is the number of
#' observations in each year, again excluding anchors (if relevant). Currently `gam_fitter()` is estimating
#' phenology and abundance metrics using a resolution of 2.4 hours (ie increments of 0.1 doy).
#' `$dat.fitted` is the data frame that was used for fitting, including any additional anchor days added.
#' `$activity.curve` is the activity curve predicted by the fitted model, provided in 2.4 hour increments for each year.
#' `$model` is the fitted mgcv model, and can be used to create plots of estimated activity curves, etc.
#' `$parms` is the fitting parameters provided (the arguments except for the initial three vectors).
#' @export
#'
#' @import rlang
#'
#' @examples
#' dat.sim = timeseries_generator(years = 1990:2000,
#'                                doy.samples = seq(105,160, by = 7),
#'                                abund.type = "exp",
#'                                activity.type = "gauss",
#'                                sample.type = "pois",
#'                                growth.rate = -0.12,
#'                                init.size = 500,
#'                                act.mean = 130,
#'                                act.sd = 15)
#' out = gam_fitter(years.vec = dat.sim$years,
#'                  doy.vec = dat.sim$doy,
#'                  count.vec = dat.sim$count,
#'                  doy.smooth = "cr",
#'                  doy.knots = 5,
#'                  years.smooth = "cr",
#'                  anchor.flag = FALSE
#' )
gam_fitter = function(years.vec,
                      doy.vec,
                      count.vec,
                      doy.smooth,
                      doy.knots,
                      years.smooth,
                      anchor.flag,
                      anchor.dist = NULL,
                      limit.to.anchor = FALSE,
                      ...){
  stopifnot(is.numeric(years.vec),
            is.numeric(doy.vec),
            is.numeric(count.vec),
            is.character(doy.smooth),
            is.numeric(doy.knots),
            is.character(years.smooth),
            is.logical(anchor.flag)
  )
  stopifnot(length(years.vec)==length(doy.vec),
            length(years.vec)==length(count.vec),
            length(doy.smooth)==1,
            length(doy.knots)==1,
            length(years.smooth)==1,
            length(anchor.flag)==1)
  stopifnot(!anchor.flag |
              (is.numeric(anchor.dist) & length(anchor.dist) == 1))
  dat.fit = data.frame(years = years.vec,
                       doy = doy.vec,
                       count = count.vec
  )
  # dat.summary =  dat.fit |>
  #   dplyr::group_by(years) |>
  #   dplyr::summarise(n = dplyr::n(),
  #                    nzero = sum(count==0)) |>
  #   dplyr::ungroup()
  dat.summary =  dat.fit |>
    dplyr::group_by( .data[["years"]] ) |>
    dplyr::summarise(n = dplyr::n(),
              nzero = sum(.data[["count"]]==0)) |>
    dplyr::ungroup()

  if(anchor.flag){
    dat.fit = add_anchors(dat.fit, anchor.dist)
  }

  ## if years are fit as random effect, need to be factor
  if(years.smooth == "re"){
    dat.fit$years.use = as.factor(as.character(dat.fit$years))
  }else{
    dat.fit$years.use = dat.fit$years
  }
  ## Fitting the gam
  if(doy.smooth == "cc"){
    k.use = c(2 + doy.knots,
              length(unique(years.vec)))
    out = mgcv::gam(count ~ te(doy, years.use,
                               k = k.use,
                               bs = c(doy.smooth, years.smooth)),
                    knots = list(doy = c(.5, 364.5)),
                    method = "REML", family = "nb",
                    data = dat.fit)
  }else{
    k.use = c(doy.knots,
              length(unique(years.vec)))
    out = mgcv::gam(count ~ te(doy, years.use,
                               k = k.use,
                               bs = c(doy.smooth, years.smooth)),
                    method = "REML", family = "nb",
                    data = dat.fit)
  }
  ## estimating metrics
  activity.curve = data.frame(expand.grid(years.use = unique(dat.fit$years.use),
                                          doy = seq(0,365, by =.1)))
  activity.curve$act = mgcv::predict.gam(out, activity.curve, type = "response")
  names(activity.curve)[names(activity.curve)=="years.use"] = "years"

  ## okay, now cut out values outside of the anchor range, if using that.
  ##
  if(limit.to.anchor){
    lims.dat = range(dat.fit$doy)
    if(anchor.flag){
      lims.dat = lims.dat + c(-anchor.dist, anchor.dist)
    }
    activity.curve = activity.curve[activity.curve$doy >= lims.dat[1] & activity.curve$doy <= lims.dat[2],]
  }

  dat.sum = gam_summarize_all(activity.curve,
                              ...)


  return(list(summary = dplyr::inner_join(dat.summary, dat.sum, by = "years"),
              model = out,
              dat.fitted = dat.fit,
              activity.curve = activity.curve,
              parms = data.frame(doy.smooth = doy.smooth,
                                 doy.knots = doy.knots,
                                 years.smooth = years.smooth,
                                 anchor.flag = anchor.flag,
                                 anchor.dist = ifelse(is.null(anchor.dist), NA, anchor.dist),
                                 limit.to.anchor = limit.to.anchor)))
}


add_anchors = function(dat.fit, anchor.dist){
  anchors = data.frame(expand.grid(years = unique(dat.fit$years),
                                   doy = range(dat.fit$doy) + c(-anchor.dist, anchor.dist),
                                   count = 0
  ))
  return(rbind(dat.fit, anchors))

}
