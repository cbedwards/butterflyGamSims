#' Plot example timeseries
#'
#' @param years Years  for which to simulate data
#' @param doy.samples Days of year to simulate observation events.
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
#'   activity.type = "Gauss",
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
  dat.detail = expand.grid(years = unique(dat$years),
                           doy = seq(xlim[1],xlim[2], by = .1))
  dat.detail = dplyr::inner_join(dat.detail, abund.merge, by = "years")
  dat.detail$act =  activity_gen(abund.vec = dat.detail$abund,
                                     doy = dat.detail$doy,
                                     activity.type = activity.type,
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
