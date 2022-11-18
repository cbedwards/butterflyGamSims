#' Fit multiple time series with multiple gam methods
#'
#' Each separate timeseries (identified with `$sim.id`) will be fit separately, once for each
#' gam method. GAM methods are defined with separate rows of the `$gam.args` arguments.
#'
#' @param timeseries_sim Collection of timeseries, including `$years`, `$doy`, `$count`, and `$sim.id`.
#'   Typically the `$timeseries` output of `timeseries_sim()` or `timeseries_load()`. If speed is
#'   of particular concern, can pass the trimmed data frame (`timeseries_sim()` generates
#'   many more columns.)
#' @param gam.args Data frame, each row corresponding to a different set of gam method arguments.
#'  Must include columns `$doy.smooth`, `$doy.knots`,`$years.smooth`,`$anchor.flag`, `$anchor.dist`.
#'  See `gam_fitter()` for details. `$anchor.dist` can be any int when `$anchor.flag==FALSE`;
#'  I recommend using a large negative number to make it easy to catch mistakes.
#'
#' @return list of two data frames. `$yearly.estimates` contains the metrics estimated
#'  for each year of each simulation for each gam method (`$sim.id` and `$gam.id`). `$gam.args`
#'  gives the associated gam arguments that map to each gam.id in `$yearly.estimates`.
#' @export
#'
#' @examples
#' # path.res = "G:/Dropbox/academia/research projects/butterfly-gam-sims/4_res"
#' # parms = list(growth.rate = c(1,-1),
#' #              init.size = 500,
#' #              act.mean = 130,
#' #              act.sd = 15)
#' # timeseries_gen(sim.name = "testsim2",
#' #                path = path.res,
#' #                nsims = 5,
#' #                years = 1990:2000,
#' #                doy.samples = seq(100,160, by = 10),
#' #                abund.type = "exp",
#' #                activity.type = "gauss",
#' #                sample.type = "pois",
#' #                sim.parms = parms)
#' # out = timeseries_load(sim.name = "testsim2",
#' #                       path = path.res)
#' # gam.args = data.frame(doy.smooth = c("cc", "cc"),
#' #                       doy.knots = c(5,3),
#' #                       years.smooth = c("cr", "cr"),
#' #                       anchor.flag = c(TRUE, TRUE),
#' #                       anchor.dist = c(10, 10))
#' # out.fit = gam_fitall(out$timeseries, gam.args)
#' # out.fit = gam_fitall_gen(sim.name = "testsim2",
#' #                          path = path.res,
#' #                          out$timeseries,
#' #                          gam.args)
#' # temp = gam_fitall_load(sim.name = "testsim2",
#' #                        path = path.res)
gam_fitall = function(timeseries,
                      gam.args){
  stopifnot(is.data.frame(gam.args))
  stopifnot(all(names(gam.args) %in%
                  c("doy.smooth","doy.knots", "years.smooth","anchor.flag","anchor.dist")))

  gamfits.list = list()
  list.ind = 1
  for(i.gamfit in 1:nrow(gam.args)){
    yearly.list = list()
    sim.ind = 1
    for(cur.sim in unique(timeseries$sim.id)){
      ind.use = timeseries$sim.id == cur.sim
      out.fit = do.call(gam_fitter, c(list(years.vec = timeseries$years[ind.use],
                                           doy.vec = timeseries$doy[ind.use],
                                           count.vec = timeseries$count[ind.use]),
                                      gam.args[i.gamfit,]))
      fit.df = out.fit$summary
      fit.df$sim.id = cur.sim
      yearly.list[[sim.ind]] = fit.df
      sim.ind = sim.ind + 1
    }
    yearly.df = do.call(rbind, yearly.list)
    yearly.df$gam.id = i.gamfit
    gamfits.list[[i.gamfit]] = yearly.df
  }
  gamfits.df = do.call(rbind, gamfits.list)
  return(list(yearly.estimates = gamfits.df, gam.args = gam.args))
}

#' Helper funciton to run gam_fitall and save results in a systematic way.
#'
#' Path and sim name should match those used in companion `timeseries_gen()` call.
#'
#' @inheritParams timeseries_load
#' @inheritParams gam_fitall
#' @param append Defaults to FALSE. If true, add results to existing gam fits csvs instead
#' of replacing them.
#'
#' @return Saves csv of yearly estimates and gam-fitting parms
#' @export
#'
#' @inherit gam_fitall examples
gam_fitall_gen = function(sim.name,
                          path,
                          timeseries,
                          gam.args,
                          append = FALSE){
  stopifnot(is.character(path),
            is.character(sim.name))
  path.use = paste0(path, "/", sim.name)
  out.fit = gam_fitall(timeseries, gam.args)
  data.table::fwrite(out.fit$yearly.estimates, paste0(path.use,"/gam-yearly.estimates.csv"), append = append)
  data.table::fwrite(out.fit$gam.args, paste0(path.use,"/gam-parms.csv"), append = append)
  cat(paste0("Gam fit results saved in ", path.use, "\n"))
}

#' Load previously saved gam fits
#'
#' Helper function for using `gam_fitall_gen()`.
#'
#' @inherit gam_fitall examples
#' @inherit gam_fitall return
#' @export
gam_fitall_load = function(sim.name,
                           path){
  path.use = paste0(path, "/", sim.name)
  yearly.estimates = data.table::fread(paste0(path.use,"/gam-yearly.estimates.csv"))
  gam.args = data.table::fread(paste0(path.use,"/gam-parms.csv"))
  return(list(yearly.estimates = yearly.estimates,
              gam.args = gam.args))
}
