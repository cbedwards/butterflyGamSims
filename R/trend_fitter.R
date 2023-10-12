#' Fit the abundance and phenology trends for summary statistics of time series
#'
#'
#' @param dat.filtered Data frame with phenology and abundance metrics for any number of years.
#' Generally just feed in `$summary` component of the results of `gam_fitter`, with or without
#' applying censorship filters.
#'
#' @param nyear.min Minimum number of years to use when fitting, after censoring years
#' with insufficient data. If fewer years are provided trends are returned as NA. This value
#' should be no greater than the length of the `years` argument, or no simulations will meet
#' this criterion.
#'
#' @return 1-row data frame with linear trends across years for each metric. For abundance,
#' `growth.rate` is the trend in log abundance across years. `gr.ci.025` and `gr.ci.975` give the 95% confidence
#' limits for the estimated growth rate.  `$nyears` is the number of unique
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
#' trend_fitter(out$summary, nyear.min = 4)
trend_fitter = function(dat.filtered,
                        nyear.min){
  res.cur = data.frame(growth.rate = -999,
                       median = -999,
                       onset = -999,
                       end = -999,
                       fp = -999)
  res.cur$nyears = length(unique(dat.filtered$years))
  if(res.cur$nyears< nyear.min){
    res.cur$growth.rate = res.cur$gr.ci.025 = res.cur$gr.ci.975 =
      res.cur$onset = res.cur$median = res.cur$end = res.cur$fp = NA
  }else{

    ## abundance trend (log)
    dat.filtered$logabund=log(dat.filtered$abund+1)
    out = stats::lm(logabund ~ years, data=dat.filtered)
    res.cur$growth.rate = stats::coef(out)[2]
    ci = stats::confint(out)
    res.cur$gr.ci.025 = ci[2,1]
    res.cur$gr.ci.975 = ci[2,2]

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
  }
  return(res.cur)
}

#' Calculate trends (and diagnostic metrics) for a given censorship method
#'
#' Applies to all sim.id and gam.id separately
#'
#' @param dat.summary Data frame with phenology and abundance metrics for any number of years.
#' Generally just feed in `$summary` component of the results of `gam_fitter`. Must include `year`,
#' `n`, `nzero`, `abund`, `onset`, `median`, `end`, `fp`, `boundary.reasonable.rel`, `boundary.reasonable.abs`,
#' `sim.id`, and `gam.id`
#' @param nobs.min Minimum number of total observations to include a year (REAL observations;
#' anchors are not counted). This is only relevant if years vary in the number of sampling days.
#' Note that in the current iteration with constant number of days sampled (the vector `doy.samples`),
#' this should be no greater than the length of `doy.samples`.
#' @param nnzero.min Minimum number of non-zero observations to include a year. This should be
#' no greater than the length of doy.samples.
#' @param bound.reasonable.rel Should we only use years with "good" fits as identified
#' using bound.reasonable.rel? See `gam_summarizer` for details.
#' @param bound.reasonable.abs Should we only use years with "good" fits as identified
#' using bound.reasonable.abs? See `gam_summarizer` for details.
#' @inheritParams trend_fitter
#'
#' @return Data frame with estimated trends (`growth.rate`, `median`, `onset`,
#' `end`, `fp`), the number of years used to estimate the trends (`nyears`),
#' the filtering criterion used to censor years before fitting trends (arguments for
#' this function), and the original number of years before censoring (`nyear.original`))
#' @export
#'
#' @examples
#' set.seed(10)
trend_method = function(dat.summary,
                        nobs.min,
                        nnzero.min,
                        nyear.min,
                        bound.reasonable.rel = F,
                        bound.reasonable.abs = F){
  dat.use = dat.summary[(dat.summary$n - dat.summary$nzero) >= nnzero.min &
                          dat.summary$n >= nobs.min &
                          (!bound.reasonable.rel |
                             dat.summary$boundary.reasonable.rel) &
                          (!bound.reasonable.abs |
                             dat.summary$boundary.reasonable.abs), ]
  iteration.ids = unique(dat.summary[,c("sim.id", "gam.id")]) #based off of ALL data
  res.list = list(); list.ind = 1
  for(i.iter in 1:nrow(iteration.ids)){
    dat.cur = dat.use[dat.use$sim.id == iteration.ids$sim.id[i.iter] &
                        dat.use$gam.id == iteration.ids$gam.id[i.iter], ]
    cur.nyear.original = nrow(dat.summary[dat.summary$sim.id == iteration.ids$sim.id[i.iter] &
                                            dat.summary$gam.id == iteration.ids$gam.id[i.iter], ])
    out.trendfit = trend_fitter(dat.cur, nyear.min = nyear.min)
    dat.cur = cbind(out.trendfit,
                    data.frame(passed.filtering = out.trendfit$nyears >= nyear.min,
                               nyear.original = cur.nyear.original,
                               nyear.min = nyear.min,
                               nobs.min = nobs.min,
                               nnzero.min = nnzero.min,
                               bound.reasonable.rel = bound.reasonable.rel,
                               bound.reasonable.abs = bound.reasonable.abs,
                               sim.id = iteration.ids$sim.id[i.iter],
                               gam.id = iteration.ids$gam.id[i.iter]))
    res.list[[list.ind]] = dat.cur; list.ind = list.ind + 1
  }
  dat.res = do.call(rbind, res.list)
  return(dat.res)
}

#' Fit trend data and aggregate with other information
#'
#' This is the workhouse function for starting with yearly metrics and fitting parameters
#'  and producing a data frame with all the pieces needed to evaluate the methods.
#'
#' @inheritParams trend_method
#' @inheritParams gam_fitall
#'
#' @return Data frame including fitted trend summary (see `trend_method()`),
#' summary of the number of "unreasonable" gam fits, the real underlying trends,
#' the exclusion methods used, and the gam methods used. Trend terms starting with "true"
#' reflect the trend estimated from the (numerically estimated) abundance and phenology
#' metrics of the underlying activity curve. Because of the randomness associated with
#' sampling (e.g. Poisson or negative binomial counts), we don't expect the data to exactly
#' reflect this activity curve even if the gam fit + exclusion rules were perfect.
#'
#' @details Fill in column-by-column description here.
#'
#'
#'
#' @export

trend_aggregator = function(dat.summary,
                            gam.args,
                            timeseries,
                            nobs.min,
                            nnzero.min,
                            nyear.min,
                            bound.reasonable.rel = T,
                            bound.reasonable.abs = T){
  ## fitting trend to "true" pattern
  timeseries = timeseries[ , -which(names(timeseries) %in% c("doy", "count", "act"))]
  timeseries = timeseries[!duplicated(timeseries),]
  timeseries$gam.id = -1
  timeseries$boundary.reasonable.rel  = TRUE
  timeseries$boundary.reasonable.abs = TRUE
  timeseries$n = 999
  timeseries$nzero = 0
  names(timeseries) = gsub("[.]true","", names(timeseries))
  # print(head(timeseries))
  trends.true = trend_method(timeseries,
                             nobs.min = 1,
                             nnzero.min = 1,
                             nyear.min = nyear.min,
                             bound.reasonable.rel = T,
                             bound.reasonable.abs = T
  )
  trends.true = trends.true[, c("sim.id", "growth.rate", "median", "onset", "end", "fp")]
  names(trends.true)[-1] = c("true.growth.rate", "true.median", "true.onset", "true.end", "true.fp")
  exclusion.methods = data.frame(expand.grid(nobs.min = nobs.min,
                                             nnzero.min = nnzero.min,
                                             nyear.min = nyear.min))
  trends.list = list()
  for(i.exclusion in 1:nrow(exclusion.methods)){
    trends.list[[i.exclusion]] = trend_method(dat.summary = dat.summary,
                                              nobs.min = exclusion.methods$nobs.min[i.exclusion],
                                              nnzero.min = exclusion.methods$nnzero.min[i.exclusion],
                                              nyear.min = exclusion.methods$nyear.min[i.exclusion],
                                              bound.reasonable.rel = bound.reasonable.rel,
                                              bound.reasonable.abs = bound.reasonable.abs)
  }
  trends.est = do.call(rbind, trends.list)
  unreasonable = dat.summary |>
    dplyr::group_by( .data[["sim.id"]]) |>
    # dplyr::summarise(unreasonable.rel = sum(!boundary.reasonable.rel),
    #                  unreasonable.abs = sum(!boundary.reasonable.abs),
    #                  unreasonable.all = sum(!boundary.reasonable.abs |
    #                                           !boundary.reasonable.rel)) |>
    dplyr::summarise(unreasonable.rel = sum(!.data[["boundary.reasonable.rel"]]),
                     unreasonable.abs = sum(!.data[["boundary.reasonable.abs"]]),
                     unreasonable.all = sum(!.data[["boundary.reasonable.rel"]] |
                                              ! .data[["boundary.reasonable.abs"]])) |>
    dplyr::ungroup()
  res = dplyr::inner_join(trends.est, trends.true, by = "sim.id")
  res = dplyr::inner_join(res, unreasonable, by = "sim.id")
  res = dplyr::inner_join(res, gam.args, by = "gam.id")
  return(res)
}


trend_aggregator_gen = function(sim.name,
                                path,
                                dat.summary,
                                gam.args,
                                timeseries,
                                nobs.min,
                                nnzero.min,
                                nyear.min,
                                append = F,
                                bound.reasonable.rel = F,
                                bound.reasonable.abs = F){
  stopifnot(is.character(path),
            is.character(sim.name))
  path.use = paste0(path, "/", sim.name)
  res = trend_aggregator(dat.summary = dat.summary,
                         gam.args = gam.args,
                         timeseries = timeseries,
                         nobs.min = nobs.min,
                         nnzero.min = nnzero.min,
                         nyear.min = nyear.min,
                         bound.reasonable.rel = bound.reasonable.rel,
                         bound.reasonable.abs = bound.reasonable.abs)
  data.table::fwrite(res,
                     paste0(path.use,"/trends-aggregated.csv"), append = append)
  return(res)
}
