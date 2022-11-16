#' Determine expected activity
#'
#' For a series of observation events defined by the abundance in that year and the doy of that
#' observation event, calculate the expected value of the activity curve. Combine with
#' `activity_sampler()` to generate
#' actual simulated counts. `activity_gen()` is the main function, which in turn can call
#' functions for specific types of activity curves.
#'
#' @param abund.vec vector of abundance index associated with the year of each sampling event.
#' Can generate with `abund_generator()`, note that it must be the same length as argument `doy`
#' @param doy day of year of each sampling event
#' @param activity.type What type activity curve to model? Currently supports the
#' Zonneveld model (`"zon"`),
#' a simple gaussian activity curve (`"gauss"`), and a bivoltine activity curve based
#' on two Gaussian curves (`"bivolt"`).
#' @param ... further arguments passed to the activity curves. See details.
#'
#' @return Numeric vector of expected activity for each observation event.
#' @export
#'
#' @import deSolve
#'
#' @details For the Gaussian distribution, the additional parameters are `act.mean` and `act.sd`,
#' defining the mean and standard deviation of the associated Guassian distribution.
#'
#' For the bivoltine distribution, the additional parameters are expanded to define the shape of
#' two activity curves, with `act.mean1` and `act.mean2` giving the means of the two gaussians,
#' and `act.sd1` and `act.sd2` giving the standard deviations of the two gaussians. Additionally,
#' `rel.size2` gives the size of thes second peak relative to the first. `rel.size = 1/2` means the
#' second peak has half the size of the first, `rel.size = 10` means the second peak has 10x the size, etc.
#' `rel.size` must be positive. Note that the code normalizes the sum of the gaussians, so rel.size will
#' not change the overall abundance.
#'
#' For the Zonneveld model,
#' peak emergence time (`zon.theta`), spread in emerge times (`beta`), and death rate (`alpha`).
#' These parameters
#' and the activity curve are taken from Zonneveld et al. 2003. For normal use of this
#' model,
#'  we would have to define initial number of larvae eclosing. However, here we set it such
#'  that the resulting abundance index equals the provided abundance index.
#'  Since area under curve = N / alpha, we set N = abundance index * alpha.
#'
#'
#' @examples
#' library(ggplot2)
#' library(tidyverse)
#' dat = expand.grid(years = seq(1:10),
#'                   doy = seq(100,150, by = 7))
#' abund.merge = data.frame(years = unique(dat$years),
#'                          abund = abund_generator(unique(dat$years),
#'                          abund.type = "exp",
#'                          growth.rate = -0.1, init.size=500)
#' )
#' dat = merge(dat, abund.merge)
#' dat$proj = activity_gen(abund.vec = dat$abund,
#'                              doy = dat$doy,
#'                              activity.type = "gauss",
#'                              act.mean = 130,
#'                              act.sd = 10)
#'
#' dat.detail = expand.grid(years = unique(dat$years),
#'                          doy = seq(50,250, by = .1))
#' dat.detail = merge(dat.detail, abund.merge)
#' dat.detail$proj = activity_gen(abund.vec = dat.detail$abund,
#'                                     doy = dat.detail$doy,
#'                                     activity.type = "gauss",
#'                                     act.mean = 130,
#'                                     act.sd = 10)
#'
#' #ggplot(data = dat, aes(x = doy, y = proj))+
#' #   geom_point()+
#' #   geom_line(data = dat.detail, aes(x = doy, y = proj))+
#' #   facet_wrap(~years)
#'
#'   #Replicating figure 1a from Zonneveld et al. 2003
#'   dat.test = data.frame(doy = seq(0,50, by =.1))
#' dat.test$abund = 55.1/0.096
#'
#' dat.test$act = activity_gen(abund.vec = dat.test$abund,
#'                                 doy = dat.test$doy,
#'                                 activity.type = "zon",
#'                                 zon.theta = 11.1,
#'                                 beta = 2.7,
#'                                 alpha = 0.096)
#' plot(dat.test$doy, dat.test$act, type='l')
activity_gen = function(abund.vec,
                        doy,
                        activity.type,
                        ...){
  stopifnot(is.numeric(abund.vec),
            is.numeric(doy),
            length(abund.vec)==length(doy),
            activity.type %in% c("zon", "gauss", "bivolt"))
  switch(activity.type,
         zon = activity_gen_zon(abund.vec, doy, ...),
         gauss = activity_gen_gaus(abund.vec, doy, ...),
         bivolt = activity_gen_bivolt(abund.vec, doy, ...))
}

## helper function for zonneveld model
zon_fun=function(t,y,parms) {
  # This function is based on the Lotka-Volterra competition model
  #state variables:
  x=y[1]
  # Parameters:
  N=parms$N; beta=parms$beta; zon.theta = parms$zon.theta; alpha = parms$alpha
  # Model:
  b = exp((t-zon.theta)/beta)
  dx=N * b / (beta * (1 + b)^2) - alpha * x
  dY=dx;
  return(list(dY));
}

#' @rdname activity_gen
activity_gen_zon = function(abund.vec, doy, ...){
  # library(deSolve)
  parms.opt = list(...)
  stopifnot("zon.theta" %in% names(parms.opt),
            "beta" %in% names(parms.opt),
            "alpha" %in% names(parms.opt))
  stopifnot(parms.opt$zon.theta>0,
            parms.opt$beta>0,
            parms.opt$alpha>=0
  )
  ## we need to streamline for odesolver, then merge back in and account for abund.vec
  dat.sub = data.frame(abund = abund.vec, doy = doy)
  doy.use = sort(unique(doy))
  parms.opt$N = 100*parms.opt$alpha #increasing by 100 for numerics reasons
  out.lsoda = data.frame(deSolve::lsoda(0, c(0, doy.use), zon_fun, parms.opt)[-1,])
  names(out.lsoda) = c("doy", "act.raw")
  dat.sub = dplyr::left_join(dat.sub, out.lsoda, by = "doy")
  dat.sub$act = dat.sub$act.raw / 100 * dat.sub$abund ## canceling out the 100 from above
  return(dat.sub$act.raw)
}

#' @rdname activity_gen
activity_gen_gaus = function(abund.vec, doy, ...){
  parms.opt = list(...)
  stopifnot("act.mean" %in% names(parms.opt),
            "act.sd" %in% names(parms.opt))
  abund.vec * stats::dnorm(doy,  mean = parms.opt$act.mean,  sd = parms.opt$act.sd)
}

#' @rdname activity_gen
activity_gen_bivolt = function(abund.vec, doy, ...){
  parms.opt = list(...)
  stopifnot("act.mean1" %in% names(parms.opt),
            "act.sd1" %in% names(parms.opt),
            "act.mean2" %in% names(parms.opt),
            "act.sd2" %in% names(parms.opt),
            "rel.size2" %in% names(parms.opt),
            parms.opt$rel.size2>0
  )
  abund.vec * 1/(1+parms.opt$rel.size2) *
    (stats::dnorm(doy,  mean = parms.opt$act.mean1,  sd = parms.opt$act.sd1) +
       parms.opt$rel.size2 * stats::dnorm(doy,  mean = parms.opt$act.mean2,  sd = parms.opt$act.sd2)
    )
}

