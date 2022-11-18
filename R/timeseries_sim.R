#' Simulate multiple time series for given scenario and parameter values.
#'
#' Dev note: should parallelize generation, I think.
#' Dev note: create helper function to run this and save results to 2 csv files,
#'   and second function to load them again.
#'
#' @param Nsims Numeber of simulations to run
#' @inheritParams timeseries_generator
#' @param sim.parms list with named components matching each of the required additional
#' parameters for the selected abundance type, activity type, and sample type.
#' See details for full list of parameter names. Each sim.parm component can either take
#' a single value (every simulation will use that value), or vector of the minimum and
#' maximum desired valeus (each simulation will use a different value drawn from a uniform
#' with those bounds).
#'
#' @return List of two dataframes, `$timeseries` and `$parms`. `$timeseries` contains all nsim
#' simulated timeseries, identified by sim.id. See `timeseries_generator()`
#' for description of all columns. `$parms` provides the individual parameter values
#' used in each simulation (relevant for cases when one or more element of `sim.parms` is
#' a 2-number vector, leading to uniform sampling of that parameter value).
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
#' #
#' # out = timeseries_load(sim.name = "testsim2",
#' #                       path = path.res)
#' # head(out$timeseries)
#' # head(out$parms)
timeseries_sim = function(nsims,
                          years,
                          doy.samples,
                          abund.type,
                          activity.type,
                          sample.type,
                          sim.parms){
  stopifnot(all(sapply(sim.parms, is.numeric)))
  stopifnot(all(sapply(sim.parms, function(x){length(x) == 1 | length(x)==2})))
  stopifnot(!is.null(names(sim.parms)))
  stopifnot(all(!(names(sim.parms) =="")))
  df.parms = data.frame(sim.id = 1:nsims)
  for(i.parm in 1:length(sim.parms)){
    if(length(sim.parms[[i.parm]]) == 2){
      df.parms[,names(sim.parms[i.parm])] = stats::runif(nsims,
                                                  min = min(sim.parms[[i.parm]]),
                                                  max = max(sim.parms[[i.parm]]))
    }else{
      df.parms[,names(sim.parms[i.parm])] = sim.parms[[i.parm]]
    }
  }
  df.list = list()
  for(i.sim in 1:nsims){
    parms.use = c(list(years = years,
                       doy.samples = doy.samples,
                       abund.type = abund.type,
                       activity.type = activity.type,
                       sample.type = sample.type),
                  unlist(df.parms[i.sim, names(df.parms) != "sim.id"]))
    cur = do.call(timeseries_generator,
                  parms.use)
    cur$sim.id = df.parms$sim.id[i.sim]
    df.list[[i.sim]]=cur
  }
  df = do.call(rbind, df.list)
  return(list(timeseries = df,
              parms = df.parms))
}

#' Simulate time series and save to csv files.
#'
#' Helper function for timeseries_sim, which also saves into a consistent structure.
#'
#' @param sim.name Name of this simulation instance.
#' @param path Filepath to save results into. Likely to remain constant across different calls of `timeseries_gen()`
#' @inheritParams timeseries_sim
#' @return Nothing is returned, but results are saved to csv files.
#' @export
#'
#' @inherit timeseries_sim examples
timeseries_gen = function(sim.name,
                          path,
                          nsims,
                          years,
                          doy.samples,
                          abund.type,
                          activity.type,
                          sample.type,
                          sim.parms){
  stopifnot(is.character(path),
            is.character(sim.name))
  path.use = paste0(path, "/", sim.name)
  dir.create(path.use)
  out = timeseries_sim(nsims = nsims,
                       years = years,
                       doy.samples = doy.samples,
                       abund.type = abund.type,
                       activity.type = activity.type,
                       sample.type = sample.type,
                       sim.parms = sim.parms)
  data.table::fwrite(out$timeseries, paste0(path.use,"/timeseries.csv"))
  data.table::fwrite(out$parms, paste0(path.use,"/parms.csv"))
  cat(paste0("Timeseries saved in ", path.use, "\n"))
}

#' Title
#'
#' @inheritParams timeseries_sim
#' @inherit timeseries_sim return
#' @export
#'
#' @inherit timeseries_sim examples
timeseries_load = function(sim.name,
                           path){
  stopifnot(is.character(path),
            is.character(sim.name))
  path.use = paste0(path, "/", sim.name)
  df.timeseries = data.table::fread(paste0(path.use,"/timeseries.csv"))
  df.parms = data.table::fread(paste0(path.use,"/parms.csv"))
  return(list(timeseries = as.data.frame(df.timeseries),
              parms = as.data.frame(df.parms)))
}

