#' Generate sample counts from expected activity
#'
#' @param activity.vec Vector of the expected activity (see `?activity_gen` to generate these).
#' @param sample.type Distribution to use when simulating censuses. The mean of this distribution
#' is determined by the underlying activity curve, and this distribution (and associated parameters)
#' determine the "sampling error" (doesn't need to strictly represent error in simulated observers).
#' Currently supported: Poisson (`"pois"`), negative binomial (`"nb"`), and zero-inflated negative binomial (`"zinb"`).
#' "nb" and "zinb" will require additional parameters (see details)
#'
#' @param ... additional arguments for the sample scheme. See details
#'
#' @details For Poisson sampling, no additional parameters are needed.
#'
#' For the negative binomial distribution, additional argument `theta` must be provided.
#' This uses the rnegbin() function of the `MASS` package, with variance of `mu + mu^2/theta`.
#'
#' For the zero-inflated negative binomial, `theta` and `pzero` must be provided. `pzero`
#' is the probability of a zero before the negative binomial is even calculated.
#' Note that for this distribution, the mean of generated counts *will be lower than
#' the expected abundance based on the activity curve*.
#'
#' @return Vector of simulated counts
#' @export
#'
#' @examples
#'
#'library(ggplot2)
#'library(tidyverse)
#'dat = expand.grid(years = seq(1:10),
#'                  doy = seq(100,150, by = 7))
#'abund.merge = data.frame(years = unique(dat$years),
#'                         abund = abund_generator(unique(dat$years),
#'                         abund.type = "exp",
#'                         growth.rate = -0.1, init.size=500)
#')
#'dat = merge(dat, abund.merge)
#'dat$act = activity_gen(abund.vec = dat$abund,
#'                       doy = dat$doy,
#'                       activity.type = "gauss",
#'                       act.mean = 130,
#'                       act.sd = 10)
#'dat$count = activity_sampler(dat$act,
#'                             sample.type = "pois")
#'
#'dat.detail = expand.grid(years = unique(dat$years),
#'                         doy = seq(50,250, by = .1))
#'dat.detail = merge(dat.detail, abund.merge)
#'dat.detail$act = activity_gen(abund.vec = dat.detail$abund,
#'                              doy = dat.detail$doy,
#'                              activity.type = "gauss",
#'                              act.mean = 130,
#'                              act.sd = 10)
#'
#' #ggplot(data = dat, aes(x = doy, y = act))+
#'  #geom_point(shape = 1, col = 'blue')+
#'  #geom_point(aes(y = count))+
#'  #geom_line(data = dat.detail)+
#'  #facet_wrap(~years)
#'
activity_sampler = function(activity.vec,
                            sample.type,
                            ...){
  stopifnot(is.numeric(activity.vec),
            sample.type %in% c("pois", "nb", "zinb"))
  switch(sample.type,
         pois = activity_samp_pois(activity.vec),
         nb = activity_samp_nb(activity.vec, ...),
         zinb = activity_samp_zinb(activity.vec, ...))
}

#' @rdname activity_sampler
activity_samp_pois = function(activity.vec){
  stats::rpois(length(activity.vec), activity.vec)
}

#' @rdname activity_sampler
activity_samp_nb = function(activity.vec, ...){
  parms.opt = list(...)
  stopifnot("theta" %in% names(parms.opt))
  MASS::rnegbin(n = length(activity.vec),
          mu = activity.vec,
          theta = parms.opt$theta)
}

#' @rdname activity_sampler
activity_samp_zinb = function(activity.vec, ...){
  parms.opt = list(...)
  stopifnot("theta" %in% names(parms.opt),
            "pzero" %in% names(parms.opt),
            parms.opt$pzero>=0,
            parms.opt$pzero <=1)
  nb = MASS::rnegbin(n = length(activity.vec),
                mu = activity.vec,
                theta = parms.opt$theta)
  zi = as.numeric(stats::runif(n = length(activity.vec))>parms.opt$pzero)
  return(nb * zi)
}
