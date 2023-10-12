#' Plot example timeseries
#'
#' Show actual activity curve (black line), expected number of butterflies on days of sampling (hollow blue circles),
#' and sampled numbers of butterflies incorporating sampling variation (black points). In full simulations,
#' black points are the simulated samples that are fit with smoothing splines.
#'
#' @param years Numeric vector of years for which which to simulate data
#' @param doy.samples Numeric vector of day of years for which to simulate censusing.
#' @inheritParams abund_generator
#' @inheritParams activity_gen
#' @inheritParams activity_sampler
#' @param xlim vector of limits for plotting
#' @param ... additional arguments for the various timeseries generation functions. See details.
#'
#' @details For all the additional arguments needed for the various component functions,
#' see `?abund_generator`, `activity_gen`, and `activity_sampler`. Unfortunately it is not
#' possible to import details from multiple functions.
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 aes_string
#'
#' @examples
#' timeseries_examples(years = 1990:2000,
#'   doy.samples = seq(105,160, by = 7),
#'   abund.type = "exp",
#'   activity.type = "gauss",
#'   sample.type = "pois",
#'   growth.rate = -0.12,
#'   init.size = 500,
#'   act.mean = 130,
#'   act.sd = 15)

timeseries_examples = function(years,
                               doy.samples,
                               abund.type,
                               activity.type,
                               sample.type,
                               xlim = c(0,365),
                               ...){
  stopifnot(is.numeric(xlim),
            length(xlim) == 2)
  if(abund.type == "both"){plot.simple = TRUE}
  dat = timeseries_generator(years = years,
                             doy.samples = doy.samples,
                             abund.type = abund.type,
                             activity.type = activity.type,
                             sample.type = sample.type,
                             ...)
  ## modifying to handle variability
  dat.detail = expand.grid(years = years,
                           doy = seq(xlim[1],xlim[2], by = .1))
  abund.merge = unique(dat[, c("years", "abund.true")])
  names(abund.merge) = c("years", "abund")
  dat.detail = dplyr::inner_join(dat.detail, abund.merge, by = "years")
  dat.detail$act = with(dat.detail,
                        activity_gen(abund.vec = abund,
                                     doy = doy,
                                     activity.type = activity.type,
                                     ...))




  dat.detail = timeseries_generator(years = years,
                                    doy.samples = seq(xlim[1],xlim[2], by = .1),
                                    abund.type = abund.type,
                                    activity.type = activity.type,
                                    sample.type = sample.type,
                                    ...)
  ggplot2::ggplot(data = dat, aes_string(x = "doy", y = "act"))+
    ggplot2::geom_point(shape = 1, col = 'blue')+
    ggplot2::geom_point(aes_string(y = "count"))+
    ggplot2::geom_line(data = dat.detail)+
    ggplot2::facet_wrap(~years)+
    ggplot2::xlab("day of year")+
    ggplot2::ylab("Activity or count")+
    ggplot2::xlim(xlim)
}

#' Plotted estimated activity curve from fitted model

#' @param dat.timeseries data frame of the original timeseries. Must contain `$count`, `$doy`, `$year` columns.
#' Typically this will be the output of `timeseries_generator()`.
#' @param dat.fitted data frame of the fitted data (which includes anchors). Typically the
#' `$dat.fitted` list entry of the output of `gam_fitter`. Must contain `$count`, `$doy`, `$year` columns.
#' @param activity.curve data frame of the fitted activity curve. Typically the `$activity.curve` list
#' entry of the output of `gam_fitter`.
#' Must contain `$act`, `$doy`, `$year` columns.
#' @param xlim.range Detaults to NULL. If given a numeric, will plot from `xlim.range` days before
#' the first doy of of observations from `dat.timeseries` to `xlim.range` days after the last doy.
#' Otherwise plot will span 0 to 365 days.
#'
#' @importFrom ggplot2 aes_string
#' @export
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
#' fit_plotter(dat.timeseries = dat.sim,
#'             dat.fitted = out$dat.fitted,
#'             activity.curve = out$activity.curve)
fit_plotter = function(dat.timeseries,
                       dat.fitted,
                       activity.curve,
                       xlim.range = NULL){
  stopifnot(all(c("count", "years", "doy") %in% names(dat.timeseries)))
  stopifnot(all(c("count", "years", "doy") %in% names(dat.fitted)))
  stopifnot(all(c("act", "years", "doy") %in% names(activity.curve)))
  stopifnot(is.null(xlim.range) | (is.numeric(xlim.range) & length(xlim.range) == 1))
  gfig = ggplot2::ggplot(data = activity.curve,
                         aes_string(x = "doy"))+
    ggplot2::geom_line(data = activity.curve,
                       aes_string(y = "act"))+
    ggplot2::geom_point(data = dat.fitted,
                        aes_string(y = "count"),
                        shape = 1)+
    ggplot2::geom_point(data = dat.timeseries,
                        aes_string(y = "count"))+
    ggplot2::facet_wrap(~years)+
    ggplot2::xlab("day of year")+
    ggplot2::ylab("Activity or count")

  if(!is.null(xlim.range)){
    gfig = gfig +
      ggplot2::xlim(range(dat.fitted$doy) + c(-xlim.range, xlim.range))
  }
  return(gfig)
}
