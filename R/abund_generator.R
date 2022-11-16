#' Determine yearly abundance index
#'
#' From a sequence of years, generates abundance indices to feed into activity_generator().
#' Currently two methods are implemented. `abund_generator_exp()` produces deterministic
#' exponential decline or growth from an initial abundnace. `abund_generator_rlnorm()` produces
#' abundances from a lognormal distribution.
#'
#'
#' @param years Years to generate abundance index for as a numeric vector.
#' For `abund_generator_exp()`, exponential growth/decline occurs between years,
#' so c(1,2,10) will give different results than c(1,2,3). For `abund_generator_rlnorm()`,
#' each year's value is IID.
#' @param abund.type What type of variation in abundance to model. Currently supported:
#' exponential growth (`"exp"`), log normal variation (`"rlnorm"`).
#' @param ... For `abund.type == "exp"`, needs `growth.rate` (yearly population growth rate) and `init.size` (abundance index in the first year).
#' For our purposes, recommend a negative growth rate. For `abund.type == "rlnorm"`,
#' needs `meanlog` and `sdlog`, the mean and standard deviation of the normal distribution on a log scale.
#'
#' @return Numeric vector of abundance indices for each year.
#' @export
#'
#' @examples
#' out = abund_generator(years = 1:10, abund.type = "exp", growth.rate = -0.1, init.size=500)
#' plot(x = 1:10, y = out)
#' out = abund_generator(years = 1:10, abund.type = "rlnorm", meanlog = 6, sdlog = 2)
#' plot(x = 1:10, y = out)
abund_generator = function(years, abund.type, ...){
  stopifnot(is.numeric(years),
            abund.type %in% c("exp", "rlnorm"))
  switch(abund.type,
         exp = abund_generator_exp(years, ...),
         rlnorm = abund_generator_rlnorm(years, ...))

}

#' @rdname abund_generator
abund_generator_exp = function(years, ...){
  parms.opt = list(...)
  stopifnot(is.numeric(years))
  stopifnot("growth.rate" %in% names(parms.opt),
            "init.size" %in% names(parms.opt))
   parms.opt$init.size * exp(parms.opt$growth.rate *(years-min(years)))
}



#' @rdname abund_generator
abund_generator_rlnorm = function(years, ...){
  parms.opt = list(...)
  stopifnot(is.numeric(years))
  stopifnot("meanlog" %in% names(parms.opt),
            "sdlog" %in% names(parms.opt))
  stats::rlnorm(length(years), parms.opt$meanlog, parms.opt$sdlog)
}
