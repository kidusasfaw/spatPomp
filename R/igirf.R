##' Iterated guided intermediate resampling filter (iGIRF)
##'
##' An implementation of a parameter estimation algorithm combining
##' GIRF with IF2, proposed by Park and Ionides (2019), following the pseudocode in Asfaw, Ionides and King (2019).
##'
##' @name igirf
##' @rdname igirf
##' @include spatPomp_class.R generics.R spatPomp.R girf.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @inheritParams spatPomp
##' @inheritParams pomp::mif2
##'
##' @param Ngirf the number of iterations of perturbed GIRF.
##'
##' @return
##' Upon successful completion, \code{igirf} returns an object of class
##' \sQuote{igirfd.spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ gives the Monte Carlo estimate of the maximum likelihood. }
##' }
##'
##' @references
##' \park2019
##'
##' \asfaw2019
NULL

rw.sd <- pomp:::safecall

setClass(
  "igirfd_spatPomp",
  contains="girfd_spatPomp",
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

##' @name igirf-spatPomp
##' @aliases igirf,spatPomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="spatPomp"),
  definition=function (data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
                        Ninter,lookahead,Nguide,
                        tol = 1e-300, max.fail = Inf,save.states = FALSE,
                        ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      igirf.internal(data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
        Ninter,lookahead,Nguide,tol = tol,
        max.fail = Inf, save.states = FALSE,...,verbose=verbose),
      error = function (e) pomp:::pStop("igirf",conditionMessage(e))
    )
  }
)

##' @name igirf-igirfd_spatPomp
##' @aliases igirf,igirfd_spatPomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="igirfd_spatPomp"),
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

    igirf(as(data,"spatPomp"), Ngirf=Ngirf, Np=Np,rw.sd = rw.sd, cooling.type = cooling.type,
         cooling.fraction.50 = cooling.fraction.50, tol=tol, Ninter=Ninter, Nguide=Nguide, lookahead=lookahead,
         h=h,theta.to.v=theta.to.v, v.to.theta=v.to.theta, ..., verbose=verbose)
  }
)

igirf.internal <- function (object,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
                            Ninter,lookahead,Nguide,
                            tol, max.fail = Inf,save.states = FALSE,...,
                            .ndone = 0L, .indices = integer(0),.paramMatrix = NULL,.gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)
  if (pomp:::undefined(object@rprocess) || pomp:::undefined(object@dmeasure))
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
  rw.sd <- pomp:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- pomp:::mif2.cooling(
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
      params=paramMatrix,
      Np=Np,girfiter=.ndone+n,cooling.fn=cooling.fn,rw.sd=rw.sd,tol=tol,max.fail=max.fail,
      verbose=verbose,.indices=.indices,...,.gnsi=gnsi)

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
    "igirfd_spatPomp",
    g,
    Ngirf=Ngirf,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

igirf.girf <- function (object, params, Ninter, lookahead, Nguide,
                        Np, girfiter, rw.sd, cooling.fn, tol, max.fail = Inf,
                        verbose, .indices = integer(0), ...,.gnsi = TRUE) {

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
  cond.loglik <- array(0, dim = c(ntimes, Ninter))

  znames <- object@accumvars

  # initialize filter guide function
  log_filter_guide_fun <- array(0, dim = Np[1])
  for (nt in 0:(ntimes-1)) {
    ## perturb parameters
    pmag <- cooling.fn(nt,girfiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')
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
    # test code
    x_with_guides <- x[,rep(1:Np[1], rep(Nguide, Np[1]))]
    tp_with_guides <- tparams[,rep(1:Np[1], rep(Nguide, Np[1]))]
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
                   params=tp_with_guides,.gnsi=gnsi)
    xx <- tryCatch(
      .Call('do_fcst_samp_var',
            object=object,
            X=Xg,
            Np = as.integer(Np[1]),
            times=times[(nt+2):(nt+1+lookahead_steps)],
            params=tp_with_guides,
            gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    fcst_samp_var <- xx
    dim(fcst_samp_var) <- c(length(spat_units(object)), lookahead_steps, Np[1])
    # end test code
    # for (p in 1:Np[1]){
    #   # find this particle's initialization and repeat in Nguide times
    #   xp = matrix(x[,p], nrow = nrow(x), ncol = Nguide, dimnames = list(nvars = statenames, ng = NULL))
    #   tparamsp = matrix(tparams[,p], nrow = nrow(tparams), ncol = Nguide, dimnames = list(params = paramnames, ng = NULL))
    #   # get all the guides for this particles
    #   Xg[,,,p] <- rprocess(object, x0=xp, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
    #                        params=tparamsp,.gnsi=gnsi)
    #   for(u in 1:length(object@units)){
    #     snames = paste0(object@unit_statenames,u)
    #     for(l in 1:lookahead_steps){
    #       hXg = apply(X=Xg[snames,,l,p, drop = FALSE], MARGIN = c(2,3,4), FUN = h, param.vec=tparams[,p])
    #       fcst_samp_var[u, l, p] = var(hXg)
    #     }
    #   }
    # }
    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      #cat(paste0("nt = ", nt, ", s = ", s, "\n"))
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
      # get prediction simulations
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
                    params=tparams,.gnsi=gnsi)
      # X is now a nvars by nreps by 1 array

      if(s>1 && length(znames)>0){
        x.znames <- x[znames,]; dim(x.znames) <- c(dim(x.znames),1)
        X[znames,,] <- X[znames,,,drop=FALSE] + x.znames
      }
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        #print(X.start)
        skel <- pomp::flow(object, x0=X.start, t0=tt[s+1], params=tparams, times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],...)
        #if(s>1 && length(znames) > 0){
          skel.start <- skel[,,1]
          X.start.znames <- X.start[znames,]
          skel.start.znames <- skel.start[znames,]
          skel.end.znames <- X.start.znames + skel.start.znames
          skel[znames,,1] <- skel.end.znames
        #}
      } else {
        skel <- X
      }
      # begin test code
      meas_var_skel <- tryCatch(
        .Call('do_theta_to_v',
              object=object,
              X=skel,
              Np = as.integer(Np[1]),
              times=times[(nt+2):(nt+1+lookahead_steps)],
              params=tparams,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      dim(meas_var_skel) <- c(length(spat_units(object)), lookahead_steps, Np[1])
      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(l in 1:lookahead_steps) fcst_var_upd[,l,] <- fcst_samp_var[,l,]*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1])
      inflated_var <- meas_var_skel + fcst_var_upd
      array.tparams <- array(NA, dim = c(dim(tparams)[1], length(spat_units(object)), Np[1], lookahead_steps), dimnames = list(tparams = rownames(tparams)))
      for(i in 1:length(spat_units(object))) array.tparams[,i,,1:lookahead_steps] <- tparams
      mmp <- tryCatch(
        .Call('do_v_to_theta',
              object=object,
              X=skel,
              vc=inflated_var,
              Np = as.integer(Np[1]),
              times=times[(nt+2):(nt+1+lookahead_steps)],
              params=array.tparams,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      mom_match_param <- mmp
      dim(mom_match_param) <- c(dim(tparams)[1], length(spat_units(object)), lookahead_steps, Np[1])
      dimnames(mom_match_param) <- list(tparam = rownames(tparams))
      # end test code
      # create measurement variance at skeleton matrix
      # meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      # for(u in 1:length(object@units)){
      #   snames = paste0(object@unit_statenames,u)
      #   for (l in 1:lookahead_steps){
      #     hskel <- sapply(1:Np[1], function(i) apply(X=skel[snames,i,l, drop = FALSE], MARGIN = c(2,3), FUN = h, param.vec = tparams[,i]))
      #     dim(hskel) <- c(1,Np[1],1)
      #     meas_var_skel[u,l,] <- sapply(1:Np[1], function(i) theta.to.v(hskel[1,i,1],tparams[,i]))
      #   }
      # }
      # fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      # for(u in 1:length(object@units)){
      #   for(l in 1:lookahead_steps){
      #     fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
      #                                 FUN = function(x) x*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1]))
      #   }
      # }
      # mom_match_param <- array(0, dim = c(dim(params)[1], length(object@units), lookahead_steps, Np[1]), dimnames = list(variable = names(params[,1]), lookahead = NULL, J = NULL))
      # inflated_var <- meas_var_skel + fcst_var_upd
      # for(p in 1:Np[1]){
      #   mom_match_param[,,,p] = apply(X=inflated_var[,,p,drop=FALSE], MARGIN=c(1,2,3), FUN = v.to.theta, param.vec = tparams[,p])[,,,1]
      # }
      log_guide_fun = vector(mode = "numeric", length = Np[1]) + 1

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead_steps <= 0) discount_denom_init = object@t0
        else discount_denom_init = times[nt+1+l - lookahead_steps]
        discount_factor = 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)
        log_dmeas_weights <- tryCatch(
          (vec_dmeasure(
            object,
            y=object@data[,nt+l,drop=FALSE],
            x=skel[,,l,drop = FALSE],
            times=times[nt+1+l],
            params=mom_match_param[,,l,],
            log=TRUE,
            .gnsi=gnsi
          )),
          error = function (e) {
            stop(ep,"error in calculation of log_dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        log_resamp_weights <- apply(log_dmeas_weights[,,1,drop=FALSE], 2, function(x) sum(x))*discount_factor
        log_guide_fun = log_guide_fun + log_resamp_weights
      }
      ## log_guide_fun[log_guide_fun < log(tol)] <- log(tol)
      log_s_not_1_weights <- log_guide_fun - log_filter_guide_fun
      if (!(s==1 & nt!=0)){
        log_weights <- log_s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        log_meas_weights <- tryCatch(
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=tparams,
            log=TRUE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of log_meas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        log_weights <- as.numeric(log_meas_weights) + log_s_not_1_weights
      }
      if (nt == ntimes-1 & s==Ninter) {
        if (any(log_weights>-Inf)) {
          coef(object,transform=TRUE) <- apply(params,1L,weighted.mean,w=exp(log_weights))
        } else {
          pomp:::pWarn("igirf","filtering failure at last filter iteration; using ",
                "unweighted mean for point estimate.")
          coef(object,transform=TRUE) <- apply(params,1L,mean)
        }
      }
      max_log_weights <- max(log_weights)
      if(max_log_weights > -Inf){
        log_weights <- log_weights - max_log_weights
        weights <- exp(log_weights)
        xx <- tryCatch(
          .Call('girf_computations',
                x=X,
                params=params,
                Np=Np[nt+1],
                trackancestry=FALSE,
                doparRS=TRUE,
                weights=weights,
                lgps=log_guide_fun,
                fsv=fcst_samp_var,
                tol=tol
          ),
          error = function (e) {
            stop(ep,conditionMessage(e),call.=FALSE) # nocov
          }
        )
        cond.loglik[nt+1, s] <- xx$loglik
        x <- xx$states
        log_filter_guide_fun <- xx$logfilterguides
        params <- xx$params
        fcst_samp_var <- xx$newfsv
      }
      else{
        cond.loglik[nt+1, s] <- log(tol)
      }
    }
  }
  #print(sum(cond.loglik))
  new(
    "girfd_spatPomp",
    object,
    Ninter=Ninter,
    Nguide=Nguide,
    lookahead=lookahead,
    paramMatrix=params,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}

