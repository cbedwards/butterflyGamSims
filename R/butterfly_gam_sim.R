#' Title
#'
#' @param sim.name
#' @param path
#' @param nsims
#' @param years
#' @param doy.samples
#' @param abund.type
#' @param activity.type
#' @param sample.type
#' @param sim.parms
#' @param gam.args
#' @param nobs.min
#' @param nnzero.min
#' @param nyear.min
#' @param append
#' @param bound.reasonable.rel
#' @param bound.reasonable.abs
#' @param append
#' @param timeseries.do
#' @param gam.do
#'
#' @return
#' @export
#'
#' @examples
#' #  path.res = "G:/Dropbox/academia/research projects/butterfly-gam-sims/4_res"
#' #  gam.args = data.frame(doy.smooth = c("cc", "cc", "cc"),
#' #                       doy.knots = c(5,3, 5),
#' #                       years.smooth = c("cr", "cr", "cr"),
#' #                       anchor.flag = c(TRUE, TRUE, FALSE),
#' #                       anchor.dist = c(10, 10, 10))
#' # sim.parms = list(growth.rate = -0.12,
#' #                  init.size = 500,
#' #                  act.mean = 130,
#' #                  act.sd = 15,
#' #                  theta = 5)
#' # res = butterfly_gam_sim(sim.name = "full-kaboodle",
#' #                   path = path.res,
#' #                   nsims = 10,
#' #                   years = 1990:2010,
#' #                   doy.samples = seq(0, 365, by = .1),
#' #                   abund.type = "exp",
#' #                   activity.type = "gauss",
#' #                   sample.type = "pois",
#' #                   sim.parms = sim.parms,
#' #                   gam.args = gam.args,
#' #                   nobs.min = 5,
#' #                   nnzero.min = 5,
#' #                   nyear.min = 3)
butterfly_gam_sim = function(sim.name,
                             path,
                             nsims,
                             years,
                             doy.samples,
                             abund.type,
                             activity.type,
                             sample.type,
                             sim.parms,
                             gam.args,
                             nobs.min,
                             nnzero.min,
                             nyear.min,
                             bound.reasonable.rel = F,
                             bound.reasonable.abs = F,
                             append = FALSE,
                             timeseries.do = TRUE,
                             gam.do = TRUE
){

  if(timeseries.do){
    print("generating time series")
    timeseries_gen(sim.name = sim.name,
                   path = path,
                   nsims = nsims,
                   years = years,
                   doy.samples = doy.samples,
                   abund.type = abund.type,
                   activity.type = activity.type,
                   sample.type = sample.type,
                   sim.parms = sim.parms)
  }
  dat.timeseries = timeseries_load(sim.name = sim.name,
                                   path = path)
  if(gam.do){
    print("fitting gams")
    gam_fitall_gen(sim.name = sim.name,
                   path = path,
                   dat.timeseries$timeseries,
                   gam.args = gam.args,
                   append = append)
  }
  dat.gamfits = gam_fitall_load(sim.name = sim.name,
                                path = path)
  print("fitting trends")
  res = trend_aggregator_gen(sim.name = sim.name,
                             path = path,
                             dat.summary = dat.gamfits$yearly.estimates,
                             gam.args = dat.gamfits$gam.args,
                             timeseries = dat.timeseries$timeseries,
                             nobs.min = nobs.min,
                             nnzero.min = nnzero.min,
                             nyear.min = nyear.min,
                             bound.reasonable.rel = bound.reasonable.rel,
                             bound.reasonable.abs = bound.reasonable.abs)
}

butterfly_gam_simulator = function(){
  #loop over years, doy.samples, abund.type, activity.type, sample.type, maybe other pieces,
  #reun each as a separate butterflygam_sim with a unique readable sim name.
}
