
abund_generator_exp = function(years, ...){
  parms.opt = list(...)
  stopifnot(is.numeric(years))
  stopifnot("growth.rate" %in% names(parms.opt),
            "init.size" %in% names(parms.opt))
   parms.opt$init.size * exp(parms.opt$growth.rate *(years-min(years)))
}

abund_generator_exp(years = 1:10,
                    growth.rate = 0,
                    init.size=50)

abund_generator_rlnorm = function(years, ...){
  parms.opt = list(...)
  stopifnot(is.numeric(years))
  stopifnot("meanlog" %in% names(parms.opt),
            "sdlog" %in% names(parms.opt))
  rlnorm(length(years), parms.opt$meanlog, parms.opt$sdlog)
}
