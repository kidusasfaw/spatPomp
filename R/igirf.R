##
## @references
## J. Park & E. L. Ionides (2018).
## Iterated GIRF.
## Arxiv 1708.08543v2
##' @include spatpomp_class.R girf.R
##

rw.sd <- pomp2:::safecall

setClass(
  "igirfd_spatpomp",
  contains="girfd_spatpomp",
  slots=c(
    Ngirf = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix'
  )
)

setGeneric(
  "igirf",
  function (data, ...)
    standardGeneric("igirf")
)

setMethod(
  "igirf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("igirf","data")
  }
)

setMethod(
  "igirf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("igirf",data)
  }
)

##' @name igirf-spatpomp
##' @aliases igirf,spatpomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="spatpomp"),
  definition=function (data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
                        Ninter,lookahead,Nguide,h,theta.to.v,v.to.theta,
                        tol = 1e-17, max.fail = Inf,save.states = FALSE,
                        ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      igirf.internal(data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
        Ninter,lookahead,Nguide,h,theta.to.v,v.to.theta,tol = tol,
        max.fail = Inf, save.states = FALSE,...,verbose=verbose),
      error = function (e) pomp2:::pStop("igirf",conditionMessage(e))
    )
  }
)

##' @name igirf-igirfd_spatpomp
##' @aliases igirf,igirfd_spatpomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="igirfd_spatpomp"),
  function (data,Ngirf,Np,rw.sd,cooling.type, cooling.fraction.50, Ninter,
            lookahead,Nguide,h,theta.to.v,v.to.theta,tol, ...,
            verbose = getOption("verbose", FALSE)) {
    if (missing(Ngirf)) Ngirf <- data@Ngirf
    if (missing(rw.sd)) rw.sd <- data@rw.sd
    if (missing(cooling.type)) cooling.type <- data@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- data@cooling.fraction.50
    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol
    if (missing(Ninter)) Ninter <- data@Ninter
    if (missing(Nguide)) Nguide <- data@Nguide
    if (missing(lookahead)) lookahead <- data@lookahead
    if (missing(h)) h <- data@h
    if (missing(theta.to.v)) theta.to.v <- data@theta.to.v
    if (missing(v.to.theta)) v.to.theta <- data@v.to.theta

    igirf(as(data,"spatpomp"), Ngirf=Ngirf, Np=Np,rw.sd = rw.sd, cooling.type = cooling.type,
         cooling.fraction.50 = cooling.fraction.50, tol=tol, Ninter=Ninter, Nguide=Nguide, lookahead=lookahead,
         h=h,theta.to.v=theta.to.v, v.to.theta=v.to.theta, ..., verbose=verbose)
  }
)

igirf.internal <- function (object,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
                            Ninter,lookahead,Nguide,h,theta.to.v, v.to.theta,
                            tol = 1e-17, max.fail = Inf,save.states = FALSE,
                            .ndone = 0L, .indices = integer(0),.paramMatrix = NULL,.gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)
  if (pomp2:::undefined(object@rprocess) || pomp2:::undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Ngirf) != 1 || !is.numeric(Ngirf) || !is.finite(Ngirf) || Ngirf < 1)
    pStop_(sQuote("Ngirf")," must be a positive integer.")
  Ngirf <- as.integer(Ngirf)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))

  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
               " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1) {
    Np <- rep(Np,times=ntimes)
  } else if (length(Np) > ntimes) {
    if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
      pWarn("igirf","Np[k] ignored for k > ",sQuote("length(time(object))"),".")
    }
    Np <- head(Np,ntimes)
  } else if (length(Np) < ntimes) {
    pStop_(sQuote("Np")," must have length 1 or ",
           sQuote("length(time(object))"),".")
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")

  Np <- as.integer(Np)
  Np <- c(Np,Np[1L])

  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  rw.sd <- pomp2:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- pomp2:::mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }

  traces <- array(dim=c(Ngirf+1,length(start)+1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Ngirf),
                                variable=c('loglik',names(start))))
  traces[1L,] <- c(NA,start)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)

  ## iterate the filtering
  for (n in seq_len(Ngirf)) {
    g <- igirf.girf(object=object,Ninter=Ninter,Nguide=Nguide,lookahead=lookahead,
      h = h,theta.to.v = theta.to.v,v.to.theta = v.to.theta,params=paramMatrix,
      Np=Np,girfiter=.ndone+n,cooling.fn=cooling.fn,rw.sd=rw.sd,tol=tol,max.fail=max.fail,
      verbose=verbose,.indices=.indices,.gnsi=gnsi)

    gnsi <- FALSE

    paramMatrix <- g@paramMatrix
    traces[n+1,-1] <- coef(g)
    traces[n+1,1] <- g@loglik
    .indices <- .indices

    if (verbose) cat("igirf iteration",n,"of",Ngirf,"completed\n")

  }

  g@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
                              .gnsi=gnsi)

  new(
    "igirfd_spatpomp",
    g,
    Ngirf=Ngirf,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

igirf.girf <- function (object, params, Ninter, lookahead, Nguide, h, theta.to.v, v.to.theta,
                        Np, girfiter, rw.sd, cooling.fn, tol = 1e-17, max.fail = Inf,
                        verbose, .indices = integer(0), .gnsi = TRUE) {

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  girfiter <- as.integer(girfiter)
  Np <- as.integer(Np)
  ep <- paste0("in ",sQuote("igirf"),": ")

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small positive number.")

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  loglik <- rep(NA,ntimes)
  eff.sample.size <- array(0, dim = c(ntimes, Ninter))
  cond.loglik <- array(0, dim = c(ntimes, Ninter))

  # initialize filter guide function
  filter_guide_fun <- array(1, dim = Np[1])
  for (nt in 0:(ntimes-1)) {
    ## perturb parameters
    pmag <- cooling.fn(nt,girfiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp2')
    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
    ## get initial states
    if (nt == 0) {
      x <- rinit(object,params=tparams)
      statenames <- rownames(x)
      paramnames <- rownames(tparams)
      nvars <- nrow(x)
    }
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    # four-dimensional array: nvars by nguide by ntimes by nreps
    Xg = array(0, dim=c(length(statenames), Nguide, lookahead_steps, Np[1]), dimnames = list(nvars = statenames, ng = NULL, lookahead = 1:lookahead_steps, nreps = NULL))
    ## for each particle get K guide particles, and fill in sample variance over K for each (lookahead value - unit - particle) combination
    fcst_samp_var <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
    for (p in 1:Np[1]){
      # find this particle's initialization and repeat in Nguide times
      xp = matrix(x[,p], nrow = nrow(x), ncol = Nguide, dimnames = list(nvars = statenames, ng = NULL))
      tparamsp = matrix(tparams[,p], nrow = nrow(tparams), ncol = Nguide, dimnames = list(params = paramnames, ng = NULL))
      # get all the guides for this particles
      Xg[,,,p] <- rprocess(object, x0=xp, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
                           params=tparamsp,.gnsi=gnsi)
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for(l in 1:lookahead_steps){
          hXg = apply(X=Xg[snames,,l,p, drop = FALSE], MARGIN = c(2,3,4), FUN = h, param.vec=tparams[,p])
          fcst_samp_var[u, l, p] = var(hXg)
        }
      }
    }
    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
      # get prediction simulations
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
                    params=tparams,.gnsi=gnsi)
      # X is now a nvars by nreps by 1 array
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp2::flow(object, x0=X.start, t0=tt[s+1], params=tparams, times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)])
      } else {
        skel <- X
      }
      # create measurement variance at skeleton matrix
      meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for (l in 1:lookahead_steps){
          hskel <- sapply(1:Np[1], function(i) apply(X=skel[snames,i,l, drop = FALSE], MARGIN = c(2,3), FUN = h, param.vec = tparams[,i]))
          dim(hskel) <- c(1,Np[1],1)
          meas_var_skel[u,l,] <- sapply(1:Np[1], function(i) theta.to.v(hskel[1,i,1],tparams[,i]))
        }
      }
      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        for(l in 1:lookahead_steps){
          fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
                                      FUN = function(x) x*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1]))
        }
      }
      mom_match_param <- array(0, dim = c(dim(params)[1], length(object@units), lookahead_steps, Np[1]), dimnames = list(variable = names(params[,1]), lookahead = NULL, J = NULL))
      inflated_var <- meas_var_skel + fcst_var_upd
      for(p in 1:Np[1]){
        mom_match_param[,,,p] = apply(X=inflated_var[,,p,drop=FALSE], MARGIN=c(1,2,3), FUN = v.to.theta, param.vec = tparams[,p])[,,,1]
      }
      guide_fun = vector(mode = "numeric", length = Np[1]) + 1

      for(l in 1:lookahead_steps){
        dmeas_weights <- tryCatch(
          vec_dmeasure(
            object,
            y=object@data[,nt+l,drop=FALSE],
            x=skel[,,l,drop = FALSE],
            times=times[nt+1+l],
            params=mom_match_param[,,l,],
            log=FALSE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        resamp_weights <- apply(dmeas_weights[,,1,drop=FALSE], 2, function(x) prod(x))
        guide_fun = guide_fun*resamp_weights
      }

      guide_fun[guide_fun < tol] <- tol
      s_not_1_weights <- guide_fun/filter_guide_fun
      if (!(s==1 & nt!=0)){
        weights <- s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        weights <- tryCatch(
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=tparams,
            log=FALSE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        weights <- as.numeric(weights)*s_not_1_weights
      }
      if (nt == ntimes-1 & s==Ninter) {
        if (any(weights>0)) {
          coef(object,transform=TRUE) <- apply(params,1L,weighted.mean,w=weights)
        } else {
          pomp2:::pWarn("igirf","filtering failure at last filter iteration; using ",
                "unweighted mean for point estimate.")
          coef(object,transform=TRUE) <- apply(params,1L,mean)
        }
      }

      xx <- tryCatch(
        .Call('girf_computations',
              x=X,
              params=params,
              Np=Np[nt+1],
              predmean=FALSE,
              predvar=FALSE,
              filtmean=FALSE,
              trackancestry=FALSE,
              doparRS=TRUE,
              weights=weights,
              gps=guide_fun,
              fsv=fcst_samp_var,
              tol=tol
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      all.fail <- xx$fail
      eff.sample.size[nt+1, s] <- xx$ess
      cond.loglik[nt+1, s] <- xx$loglik
      x <- xx$states
      filter_guide_fun <- xx$filterguides
      params <- xx$params
      fcst_samp_var <- xx$newfsv
    }
  }
  print(sum(cond.loglik))
  new(
    "girfd_spatpomp",
    object,
    Ninter=Ninter,
    Nguide=Nguide,
    lookahead=lookahead,
    h = h,
    theta.to.v = theta.to.v,
    v.to.theta = v.to.theta,
    paramMatrix=params,
    eff.sample.size=eff.sample.size,
    saved.states=x,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}

