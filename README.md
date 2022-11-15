
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

# Project structure

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

#### Function to sample activity curve in each year

Notes: uses parameters for “sampling design”

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
