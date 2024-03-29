% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timeseries_sim.R
\name{timeseries_load}
\alias{timeseries_load}
\title{Read timeseries from files}
\usage{
timeseries_load(sim.name, path)
}
\arguments{
\item{sim.name}{Name of this simulation instance.}

\item{path}{Filepath to save results into. Likely to remain constant across different calls of \code{timeseries_gen()}}
}
\value{
List of two dataframes, \verb{$timeseries} and \verb{$parms}. \verb{$timeseries} contains all nsim
simulated timeseries, identified by sim.id. See \code{timeseries_generator()}
for description of all columns. \verb{$parms} provides the individual parameter values
used in each simulation (relevant for cases when one or more element of \code{sim.parms} is
a 2-number vector, leading to uniform sampling of that parameter value).
}
\description{
Read timeseries from files
}
\examples{
# path.res = "G:/Dropbox/academia/research projects/butterfly-gam-sims/4_res"
# parms = list(growth.rate = c(1,-1),
#              init.size = 500,
#              act.mean = 130,
#              act.sd = 15)
# timeseries_gen(sim.name = "testsim2",
#                path = path.res,
#                nsims = 5,
#                years = 1990:2000,
#                doy.samples = seq(100,160, by = 10),
#                abund.type = "exp",
#                activity.type = "gauss",
#                sample.type = "pois",
#                sim.parms = parms)
#
# out = timeseries_load(sim.name = "testsim2",
#                       path = path.res)
# head(out$timeseries)
# head(out$parms)
}
