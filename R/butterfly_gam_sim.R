#' Generate time series, fit them with gams, summarize gam fits, and estimate associated trends.
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


## can take lists of doy.samples vectors, years vectors
## BIG NOTE: this function gives CRAZY errors if timeseries.do and gam.do are FALSE but
##   the appropriate timeseries and gam outputs have not been produced (ie changing
##   other arguments and leaving those at FALSE. gam.do and timeseries.do should
##   generally be TRUE; FALSE is only useful for rapid debugging of some steps in the code)

#example:
#
# path.res = "G:/Dropbox/academia/research projects/butterfly-gam-sims/4_res"
# gam.args = data.frame(doy.smooth = c("cc", "cc", "cc"),
#                       doy.knots = c(5,3, 5),
#                       years.smooth = c("cr", "cr", "cr"),
#                       anchor.flag = c(TRUE, TRUE, FALSE),
#                       anchor.dist = c(10, 10, 10))
# sim.parms = list(growth.rate = -0.12,
#                  init.size = 500,
#                  act.mean = 130,
#                  act.sd = 15,
#                  theta = 5)
#
# ## Notes to self: (a) getting crazy error messages when trying to run
# ##   simulator. Looks like it's interpretting filepath as a command somehow?
# ##   seems to be showing up specifically in the timeseries_load and gam_fitall_load step.
# ##   see locations of print items.
# ##   (b) using the sim + gam inds is producing what looks like misbehavior.
# ##   But check with different parameterization, which I will only find
# ##   once I've go the weird "taking input= as..." error cleared out.
#
# res = butterfly_gam_simulator(sim.name = "multi-parm-test2",
#                               path = path.res,
#                               nsims = 5,
#                               years = 1990:2010,
#                               doy.samples = seq(100,160, by = 10),
#                               abund.type = "exp",
#                               activity.type = "gauss",
#                               sample.type = "pois",
#                               sim.parms = sim.parms,
#                               gam.args = gam.args,
#                               nobs.min = c(8, 6),
#                               nnzero.min = 3,
#                               nyear.min = 6,
#                               timeseries.do = TRUE,
#                               gam.do = TRUE)


butterfly_gam_simulator = function(sim.name,
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
                                   timeseries.do = TRUE,
                                   gam.do = TRUE,
                                   safe = FALSE){
#  stopifnot(!file.exists(paste0(path, "/", sim.name) | ) #figure out stopifnot to ensure that safe causes error out instead of overwrite.
#  #figure out overwriting successfully.
  # if(timeseries.do & gam.do & !safe){
  #   unlink(paste0(path, "/", sim.name))
  #   dir.create(paste0(path, "/", sim.name))
  # }else{
    dir.create(paste0(path, "/", sim.name), showWarnings = FALSE)
  # }
  if(!is.list(doy.samples)){
    doy.samples = list(doy.samples)
  }
  if(!is.list(years)){
    years = list(years)
  }
  bigparms = as.data.frame(expand.grid(years.ind = 1:length(years),
                                       doy.samples.ind = 1:length(doy.samples),
                                       abund.type = unique(abund.type),
                                       activity.type = unique(activity.type),
                                       sample.type = unique(sample.type),
                                       nobs.min = unique(nobs.min),
                                       nnzero.min = unique(nnzero.min),
                                       nyear.min = unique(nyear.min),
                                       stringsAsFactors = FALSE))
  for(i in 1:nrow(bigparms)){
    sim.name.cur = paste0(sim.name,"/", sim.name, "-", paste0(bigparms[i,], collapse = '-'))
    print(bigparms[i,])
    butterfly_gam_sim(sim.name = sim.name.cur,
                      path = path,
                      nsims = nsims,
                      years = years[[bigparms$years.ind[i]]],
                      doy.samples = doy.samples[[bigparms$doy.samples.ind[i]]],
                      abund.type = bigparms$abund.type[i],
                      activity.type = bigparms$activity.type[i],
                      sample.type  = bigparms$sample.type[i],
                      sim.parms = sim.parms,
                      gam.args = gam.args,
                      nobs.min = bigparms$nobs.min[i],
                      nnzero.min = bigparms$nnzero.min[i],
                      nyear.min = bigparms$nyear.min[i],
                      bound.reasonable.rel = F,
                      bound.reasonable.abs = F,
                      append = FALSE,
                      timeseries.do = timeseries.do,
                      gam.do = gam.do)
  }
  write.csv(bigparms,
            file = paste0(path,"/", sim.name, "/parms.csv"),
            row.names = FALSE)
  file.doy = paste0(path,"/", sim.name, "/doy-list.txt")
  cat("Doy samples:\n", file = file.doy)
  for(i in 1:length(doy.samples)){
    cat(paste0("[[",i, "]] = ",paste0(doy.samples[[i]], collapse = ", "), "\n"), file = file.doy)
  }
  file.years = paste0(path,"/", sim.name, "/years-list.txt")
  cat("Years:\n", file = file.years)
  print(class(years))
  for(i in 1:length(years)){
    cat(paste0("[[",i, "]] = ", paste0(years[[i]], collapse = ", "), "\n"), file = file.years)
  }
}
