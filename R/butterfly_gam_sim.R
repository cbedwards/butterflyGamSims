#' Run full simulation
#'
#' This is the main user function for this package. `butterfly_gam_sim()` takes timeseries parameters,
#' gam-fitting parameters, and censoring parameters, generates appropriate time series, fits
#' those time series with appropriate gams, estimates phenology metrics and abundance index
#' for each year, then looks for linear (log-linear in the case of abundance) trends across
#' years after applying the appropriate filters of `nobs.min`, `nnzero.min`, `nyear.min`.
#' Results are saved with a consistent file structure in the `sim.name` folder of the
#' `path` directory. For easy access to these simulated files, see `butterfly_load()`.
#'
#' @param sim.name Name for saving simulations of this combination of parameters
#' @param path File path to save simulated data and results into
#' @param nsims Number of simulations to carry out
#' @inheritParams timeseries_examples
#' @inheritParams timeseries_sim
#' @inheritParams gam_fitall
#' @inheritParams trend_fitter
#' @inheritParams trend_method
#' @param append For debugging purposes only. Leave at default value of `FALSE`.
#' @param timeseries.do For debugging purposes only. Leave at default value of `TRUE`
#' @param gam.do For debugging purposes only. Leave at default value of `TRUE`
#'
#' @details To allow simulations to capture a range of possible population patterns, activity curves and sampling
#' distributions, `sim.parms` has variable requirements based on `sample.type`, `abund.type`, and `activity.type`.
#' For `sample_type = "nb"`, theta must be included in `sim.parms`. For `sample_type = "zinb"`,
#' `theta` and `pzero` must be included in `sim.parms`.#' See `activity_sampler()` for details.
#'
#' For `activity.type = "gaus"`, `sim.parms` must include `act.mean` (peak day) and `act.sd` (peak width).
#' For `"zon"`, `sim.parms` must include `zon.theta` (peak day relative to day 0), `t0` (the "initial" day of the zonneveld model,
#' which is the final day that populations are 0), `beta` (spread), and `alpha` (death rate).
#' For "bivolt", `sim.parms` must include `act.mean1` and `act.mean2` (peak days), `act.sd1` and
#' `act.sd2` (peak widths), and `rel.size` (ratio of peak sizes). See `activity_gen()` for details.
#'
#' `sim.parms` must also contain information about the underlying population. For `abund.type = "exp"` (exponential growth),
#' this means the initial abundance index `init.size` and the growth rate `growth.rate`. For `abund.type = "rlnorm"`,
#' (log-nromal random draws), this means `meanlog` and `sdlog`. See `abund_generator` for details.
#'
#' @return No object is returned, but a series of files are generated in `path/sim.name`,
#'  which can be easily loaded with `butterfly_load`.
#' @export
#'
#' @inherit butterfly_load examples

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
  file.doy = paste0(path,"/", sim.name, "/doy-list.txt")
  cat("Doy samples:\n", file = file.doy)
  cat(paste0(paste0(doy.samples, collapse = ", "), "\n"), file = file.doy, append=T)
  file.years = paste0(path,"/", sim.name, "/years-list.txt")
  cat("Years:\n", file = file.years)
  print(class(years))
  cat(paste0(paste0(years, collapse = ", "), "\n"), file = file.years, append=T)
}


#' Load time series, yearly gam estimates, and trend information from simulated parameters
#'
#' @inheritParams butterfly_gam_sim
#'
#' @return List with all relevant data from butterfly simulations.
#' `$timeseries` gives the simulated data as if reported by butterfly survey teams (`years`, `doy`, `count`).
#' Additionally, because we simulated the data, we know what the TRUE underlying behavior was. The `act` column
#' shows the actual activity for each day based on the activity curve (`count` is a random sampling with that mean), and the `*.true` columns
#' show the underlying behavior for that year based on the actual activity curve:
#' the true total abundance for that year as determined by `abund.type` and associated parameters (`abund.true`),
#' the true onset (0.1 quantile, `onset.true`), median (`median.true`), end (0.9 quantile, `end.true`),
#' and flight period (days between 0.1 and 0.9 quantiles, `fp.true`) for that year, as determined by the
#' activity curve. The `sim.id` column differentiates different instances of the simulation (IE different
#' simulated data sets using the same parameters).
#'
#'  `$yearly.metrics` contains summary information for each year of simualted data for each simulation,
#'  as well as the estimated abundance index and phenology metrics based on the gam fit. `sim.id` distinguishes
#'  between different instances of the simulation (as in `$timeseries`), and `gam.id` distinguishes between
#'  different sets of gam-fitting parameters, so that we can compare how different approaches to fitting
#'  produce different yearly estimates. `years` gives each individual year, `n` goves the total number
#'  of simulated data points in that year, `nzero` gives the number of data points with 0 counts (high nzero
#'  can lead to unreliable gam fits), and `abund`, `onset`, `median`, `end` and `fp` are estimated metrics from
#'  the gam predictions of the activity curve for that year. The `boundary.reasonable.rel` and `boundary.reasonable.abs`
#'  are attempts at automatically detecting poor gam fits (in the extreme, predicted activity curves that are U-shaped).
#'  In general, `FALSE` indcates an unreasonable fit; see `gam_summarizer()` for details.
#'
#'  `$trend` contains the estimated trends in abundance and phenology, the true trends in abundance and phenology, a ton of additional summary information
#'  including gam-fitting parameters. `growth.rate` is the slope of the linear trend in log abundance index across years,
#'  which is generally a good approximation for population growth rate (assuming similar observation
#'  effort and detectability across years). `median`, `onset` `end`, and `fp` are the linear trends in
#'  estimated phenology metrics across years. `nyears` shows the number of years that passed the censoring filter
#'  and so were used in estimating these trends, and `passed.filtering` shows whether `nyears` was
#'  above the minimum number of years defined by the `nyear.min` argument (if `passed.filtering` is `FALSE`,
#'  estimated trends will all be `NA`). `nyear.original` gives the number of years originally simulated,
#'  and `nyear.min`, `nyobs.min`, and `nnzero.min` show the filtering criterion parameters used.
#'  `bound.reasonable.rel` and `bound.reasonable.abs` show whether these two metrics were used
#'  to further filter data-fitting. `sim.id` and `gam.id` identify the timeseries instance and the
#'  gam arguments instance for this row of data. `true.*` columns show the TRUE trends in log abundance index,
#'  median, onset, end, and fp. Note that because `true.fp` is calculated numberically, it is slightly off from zero.
#'  Currently simulations are designed to generate a constant phenology through time (with abundance shifting and
#'  variable sampling noise), so the true trends in phenology should all be ~0. `unreasonable.*` give the
#'  number of years that were identified as unreasonable based on relative thresholds, absolute thresholds, or either.
#'  These are calculated independently of whether these classifications were used as an additional filter
#'  before calculating the trend, and may be useful for identifying gam parameterizations that are more
#'  or less effective. The final columns are the gam parameters used: the type of smooth (cylic or cubic regression)
#'  used for smoothing in the day-of-year direction (`doy.smooth`) and the year direction (`years.smooth`),
#'  the number of knots used in the doy direction (the years direction always uses one knot per year), and
#'  whether or not "anchoring zeroes" should be added outside of the observed data (`anchor.flag`) and
#'  how far out to place them (`anchor.dist`)
#'
#'  `$parm.timeseries` gives the parameters used to generate each time series, and `$parms.gam` gives the
#'  parameters used for fitting the gams to each time series. Note that when `sim.parms` is given 2-value vectors for
#'  one or more list entry, each time series has a different value for that entry drawn from a uniform
#'  distribution with the specified values as bounds. When working with `$timeseries` and
#'  `$yearly.metrics`, `sim.id` and `gam.id` can be used to match individual simulations or gam-fitting parameters
#'  to the associated timeseries or yearly estimates.
#'
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
#' # butterfly_gam_sim(sim.name = "full-kaboodle",
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
#' # res = butterfly_load(sim.name = "full-kaboodle",
#' #                      path = path.res)
butterfly_load = function(sim.name,
                            path){
  path.use = paste0(path, "/", sim.name)
  res = list()
  ## load timeseries
  out.ts = timeseries_load(sim.name, path)
  out.gam = gam_fitall_load(sim.name, path)
  res[["timeseries"]] = out.ts$timeseries
  res[["yearly.metrics"]] = out.gam$yearly.estimates
  res[["trends"]] = data.table::fread(paste0(path.use,"/trends-aggregated.csv"))
  res[["parms.timeseries"]] = out.ts$parms
  res[["parms.gam"]] = out.gam$gam.args
  return(res)
}
