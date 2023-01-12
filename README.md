
<!-- README.md is generated from README.Rmd. Please edit that file -->

# butterflyGamSims

## Installation

You can install the development version of butterflyGamSims from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cbedwards/butterflyGamSims")
```

## Description

**butterflyGamSims** can be used to simulate butterfly transect data
with known properties, fit possible generalized additive models to
estimate activity curves, calculate abundance index and phenology
metrics from these estimated activity curves, and using one or more data
filtering protocols, estimate across-year trends in abundance and
phenology. My primary motivation for creating this package was to
compare estimated and real trends in abundance and phenology for a range
of gam approaches and data censoring options to identify potential
biases when fitting real data. In particular, in exploring real data, I
found that for sufficiently sparse data in which abundance declined
across a time series, apparent trends in phenology became unreliable
without some form of removing years with insufficient data. This package
allows systematic exploration of *what* rules for removing data are
sufficient to avoid notable bias.

For this purpose, `butterfly_gam_sim()` and `butterfly_load()` are the
key functions. `butterfly_gam_sim()` takes all the arguments needed to
define the data to be simulated, the gam-fitting parameters, and the
data filtering protocol. See `?butterfly_gam_sim` for details on all the
arguments (and a working example). `butterfly_load()` reads in the files
saved by a run of `butterfly_gam_sim()` and organizes them into a single
list; `?butterfly_load` contains all the relevant information for
interpretting the simulated time series, estimated yearly metrics, and
estimated trends.

Because **butterflyGamSims** necessarily simulates butterfly time series
data, the timeseries simulation function `timeseries_sim()` may be
helpful for other projects (ie simulated data for teaching/training,
other methods testing projects). Note that `timeseries_sim()` returns a
set of simulated time series; depending on your needs,`timeseries_gen()`
(which saves them to a file) may be more helpful.

For visualizing simulated timeseries and associated fitted gam models,
`timeseries_examples()` and `fit_plotter()` may be helpful. As a simple
example, here we simulate a 10-year time series with a gaussian
underlying activity curve centered on day 130 with standard deviation
15, a starting abundance index of 500, a moderate exponential decline in
population across years, and 7 census per year with poisson sampling
error.

``` r
library(butterflyGamSims)
timeseries_examples(years = 1990:2000,
                    doy.samples = seq(105,160, by = 7),
                    abund.type = "exp",
                    activity.type = "gauss",
                    sample.type = "pois",
                    growth.rate = -0.12,
                    init.size = 500,
                    act.mean = 130,
                    act.sd = 15)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

And here we fit a similar timeseries with a simple gam (using a tensor
product smooth with a cubic regression smooth in both the day-of-year
and the year dimensions),

``` r
dat.sim = timeseries_generator(years = 1990:2000,
                               doy.samples = seq(105,160, by = 7),
                               abund.type = "exp",
                               activity.type = "gauss",
                               sample.type = "pois",
                               growth.rate = -0.12,
                               init.size = 500,
                               act.mean = 130,
                               act.sd = 15)
out = gam_fitter(years.vec = dat.sim$years,
                 doy.vec = dat.sim$doy,
                 count.vec = dat.sim$count,
                 doy.smooth = "cr",
                 doy.knots = 5,
                 years.smooth = "cr",
                 anchor.flag = FALSE
)
fit_plotter(dat.timeseries = dat.sim,
            dat.fitted = out$dat.fitted,
            activity.curve = out$activity.curve)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
