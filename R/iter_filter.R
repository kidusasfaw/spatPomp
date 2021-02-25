mif2.cooling <- function (type, fraction, ntimes) {
  switch(
    type,
    geometric={
      factor <- fraction^(1/50)
      function (nt, m) {
        alpha <- factor^(nt/ntimes+m-1)
        list(alpha=alpha,gamma=alpha^2)
      }
    },
    hyperbolic={
      if (fraction < 1) {
        scal <- (50*ntimes*fraction-1)/(1-fraction)
        function (nt, m) {
          alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
          list(alpha=alpha,gamma=alpha^2)
        }
      } else {
        function (nt, m) {
          list(alpha=1,gamma=1)
        }
      }
    }
  )
}

perturbn.kernel.sd <- function (rw.sd, time, paramnames) {

  if (is.matrix(rw.sd)) return(rw.sd)
  if (is(rw.sd,"safecall")) {
    enclos <- rw.sd@envir
    rw.sd <- as.list(rw.sd@call)[-1L]
  } else {
    pStop_(sQuote("rw.sd")," should be specified using the ",sQuote("rw.sd"),
           " function. See ",sQuote("?mif2"),".")
  }
  if (is.null(names(rw.sd)) | any(names(rw.sd)==""))
    pStop("rw.sd","parameters must be referenced by name.")
  if (!all(names(rw.sd) %in% paramnames)) {
    unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
    pStop_("the following parameter(s), ",
           "given random walks in ",sQuote("rw.sd"),", are not present in ",
           sQuote("params"),": ",paste(sapply(unrec,sQuote),collapse=","),".")
  }
  ivp <- function (sd, lag = 1L) {
    sd*(seq_along(time)==lag)
  }
  sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp),enclos=enclos)
  for (n in names(sds)) {
    len <- length(sds[[n]])
    if (len==1) {
      sds[[n]] <- rep(sds[[n]],length(time))
    } else if (len!=length(time)) {
      pStop_(sQuote("rw.sd")," spec for parameter ",sQuote(n),
             " does not evaluate to a vector of the correct length (",
             sQuote("length(time(object))"),"=",length(time),").")
    }
  }
  do.call(rbind,sds)
}
