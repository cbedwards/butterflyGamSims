
<!-- README.md is generated from README.Rmd. Please edit that file -->

# butterflyGamSims

<!-- badges: start -->
<!-- badges: end -->

The goal of butterflyGamSims is to determine reasonable protocols for
fitting butterfly timeseries data with generalized additive models.

## Installation

You can install the development version of butterflyGamSims from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cbedwards/butterflyGamSims")
```

# Using this package:

# Project structure (internal use/trakcing)

## Overview

Simulate time series data from known yearly activity curves that change
across years in some predetermined fashion. Actual data will be counts
simulated from some probability distribution with mean determined by the
activity curve for that day, with the days of observation provided as
parameters.

## Parameters

## Key experimental designs / questions to explore/

## Code structure

### Time series simulation

In combination, the following functions simulate a time series:

#### Function to generate yearly population abundance index

`abund_generator.R` contains two functions for this.

#### Function to generate activity curve in each year

`activity_gen()`

#### Function to sample activity curve in each year

`activity_sampler()`

#### Function to calc actual pheno + abund metrics for each year from activity curve

`gam_summarizer` feeds into `timeseries_truepheno`, which feeds into
`timeseries_generator`. The end result is that `timeseries_generator`
provides the “correct” answers for abundance index (not accounting for
0s when using `sample.type="zinb"`) and phenology metrics.

#### Function to plot example time series

`timeseries_example`.

### Analysis

In combination, the following functions calculate metrics from a single
time series

#### Function to fit time series data with gam model

Note: want to be able to toggle on/off the use of anchors, the use of
different gam structures

#### Function to generate yearly metrics from FITTED activity curve

See `gam_summarizer`. Note that this also does basic “badfit”
diagnostics in the form of `boundary.reasonable.rel` and
`boundary.reasonable.abs`.

#### Function to fit across-year trends

Note: here we want to be able to try several criterion for avoiding
bias, including excluding some site-years for phenology trends based on
number of non-zero observations.

#### function to integrate the above for a given timeseries

### Simulating

Paralellized function to call the time series simulation and the
analyses. Propose saving the time series as CSV(s?) as well. Generate a
few example plots of the time series, and time series with fits.

Summarize effectiveness of gam fit: average error in abundance,
phenology measures. Bias in same.

#### Function to produce example plots of time series and fits

#### Function to summarize error and bias.

#### Function to create summary graphs

#### Save results in a meaningful way

#### Save parameters as metadata as well.

### 
