##
## @references
## J. Park & E. L. Ionides (2018).
## A guided intermediate resampling particle filter for inference on high dimensional systems.
## Arxiv 1708.08543v2
##

setClass(
  "girfd_pomp",
  contains="pomp",
  slots=c(
    S="integer",
    K="integer",
  # remaining slots are from pfilter
    pred.mean="array",
    pred.var="array",
    filter.mean="array",
    filter.traj="array",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    saved.states="list",
    Np="integer",
    tol="numeric",
    nfail="integer",
    loglik="numeric"
  ),
  prototype=prototype(
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    pred.var=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    saved.states=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=as.integer(NA),
    loglik=as.double(NA)
  )
)

setGeneric(
  "girf",
  function (data, ...)
    standardGeneric("girf")
)

setMethod(
  "girf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("girf","data")
  }
)

setMethod(
  "girf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("girf",data)
  }
)

##' @name girf-spatpomp
##' @aliases girf,spatpomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(data="spatpomp"),
  definition=function (
    data,
    Np,
    Ninter,
    lookahead,
    Nguide,
    h,
    theta_to_v,
    v_to_theta,
    tol = 1e-17, max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      girf.internal(
        data,
        Np,
        Ninter,
        lookahead,
        Nguide,
        h,
        theta_to_v,
        v_to_theta,
        tol = tol, max.fail = Inf,
        pred.mean = FALSE,
        pred.var = FALSE,
        filter.mean = FALSE,
        filter.traj = FALSE,
        save.states = FALSE,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp2:::pStop("girf",conditionMessage(e))
    )
  }
)

##' @name girf-girfd_pomp
##' @aliases girf,girfd_pomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(data="girfd_pomp"),
  function (data, Np, tol,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    girf(as(data,"pomp"),Np=Np,tol=tol,
      ...,verbose=verbose)

  }
)

girf.internal <- function (object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        h,
        theta_to_v,
        v_to_theta,
        tol = 1e-17, max.fail = Inf,
        pred.mean = FALSE,
        pred.var = FALSE,
        filter.mean = FALSE,
        filter.traj = FALSE,
        save.states = FALSE,
        cooling, cooling.m,
        .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)
  ep <- paste0("in ",sQuote("girf"),": ")

  if (pomp2:::undefined(object@rprocess) || pomp2:::undefined(object@dmeasure))
    pomp2:::pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  save.states <- as.logical(save.states)

  params <- coef(object)
  params_matrix <- matrix(params,nrow=length(params), ncol = Np[1])
  rownames(params_matrix) <- names(params)
  times <- time(object,t0=TRUE)
  t0 <- times[1]
  ntimes <- length(times)-1

  if (missing(Np) || is.null(Np)) {
    pomp2:::pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pomp2:::pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pomp2:::pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pomp2:::pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pomp2:::pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pomp2:::pStop_(sQuote("tol")," should be a small positive number.")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states || filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- array(0, dim = c(ntimes, Ninter))
  eff.sample.size <- array(0, dim = c(ntimes, Ninter))
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(data=numeric(1),dim=c(nvars,1,ntimes+1),
      dimnames=list(variable=statenames,rep=1,time=times))
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  # initialize filter guide function
  filter_guide_fun <- array(1, dim = Np[1])
  for (nt in 0:(ntimes-1)) { ## main loop
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
      # get all the guides for this particles
      Xg[,,,p] <- rprocess(object, xstart=xp, times=times[(nt+1):(nt+1+lookahead_steps)],
               params=params,offset=1L,.gnsi=gnsi)
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for(l in 1:lookahead_steps){
          hXg = apply(X=Xg[snames,,l,p, drop = FALSE], MARGIN = c(2,3,4), FUN = h, obj=object)
          fcst_samp_var[u, l, p] = var(hXg)
        }
      }
    }
    # print("=============================\n")
    # print("Guide simulations")
    # print(Xg)
    # print("=============================\n")
    # print("Forecast sample variance")
    # print(fcst_samp_var)
    # print("=============================\n")
    # print("xstart")
    # print(x)
    # print("=============================\n")
    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      # get prediction simulations
      X <- rprocess(object,xstart=x,times=c(tt[s], tt[s+1]),
                    params=params,offset=1L,.gnsi=gnsi)
      # print("Prediction particles\n")
      # print(paste0("s = ",s, ", nt = ",nt, "\n"))
      # print(paste0("X sub", nt, ",", s, ","))
      # print(X)
      # print("=============================\n")

      # X is now a nvars by nreps by 1 array
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp2::flow(object, xstart=X.start, params=params_matrix, times = c(tt[s+1], times[(nt + 1 + 1):(nt + 1 + lookahead_steps)]), offset = 1)
      } else {
        skel <- X
      }
      # print("Skeleton\n")
      # print(paste0("s = ",s, ", nt = ",nt, "\n"))
      # print(paste0("Skel = ", "\n"))
      # print(skel)
      # print("=============================\n")
      # create measurement variance at skeleton matrix
      meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for (l in 1:lookahead_steps){
          hskel <- apply(X=skel[snames,,l, drop = FALSE], MARGIN = c(2,3), FUN = h, obj = object)
          meas_var_skel[u,l,] <- apply(X=hskel, MARGIN = c(1,2), FUN = theta_to_v, obj = object)
        }
      }
      # print("Measurement Variance at Skeleton\n")
      # print(paste0("s=",s, ", nt = ",nt, "\n"))
      # print(paste0("Measurement Variance at Skeleton = ", "\n"))
      # print(meas_var_skel)
      # print("=============================\n")
      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        for(l in 1:lookahead_steps){
          fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
                                      FUN = function(x) x*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1]))
        }
      }
      # print("Forecast Variance Update\n")
      # print(paste0("s=",s, ", nt = ",nt, "\n"))
      # print(paste0("Forecast Variance Update = ", "\n"))
      # print(fcst_var_upd)
      # print("=============================\n")
      mom_match_param <- array(0, dim = c(length(params), length(object@units), lookahead_steps, Np[1]), dimnames = list(params = names(params), lookahead = NULL, J = NULL))
      inflated_var <- meas_var_skel + fcst_var_upd
      mom_match_param = apply(X=inflated_var, MARGIN=c(1,2,3), FUN = v_to_theta, obj = object)
      # print("Moment Matched Parameter\n")
      # print(paste0("s=",s, ", nt = ",nt, "\n"))
      # print(paste0("Moment Matched Parameter = ", "\n"))
      # print(mom_match_param)
      # print("=============================\n")
      # guide functions as product (so base case is 1)
      guide_fun = vector(mode = "numeric", length = Np[1]) + 1
      # print("Guide function evolution\n")
      # print(paste0("s=",s, ", nt = ",nt, "\n"))
      # print(paste0("Guide function starts off at = ", "\n"))
      # print(guide_fun)
      # print("=============================\n")
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
        # print(paste0("s=",s, ", nt = ",nt, "l = ", l, "\n"))
        # print(paste0("Guide update for l = ", l, "="))
        # print(dmeas_weights)
        # print("=============================\n")
        resamp_weights <- apply(dmeas_weights[,,1,drop=FALSE], 2, function(x) prod(x))
        # print(paste0("s=",s, ", nt = ",nt, "l = ", l, "\n"))
        # print(paste0("Resampling weights for l = ", l, "="))
        # print(resamp_weights)
        # print("=============================\n")
        guide_fun = guide_fun*resamp_weights
      }
      # print(paste0("s=",s, ", nt = ",nt, "Pre-tolerance imputation guide function", "\n"))
      # print(paste0("Guide function for", "s=",s, ", nt = ",nt, "="))
      # print(guide_fun)
      # print("=============================\n")
      guide_fun[guide_fun < tol^(lookahead*length(object@units))] <- tol^(lookahead*length(object@units))
      # print(paste0("s=",s, ", nt = ",nt, "Post-tolerance imputation guide function", "\n"))
      # print(paste0("Guide function for", "s=",s, ", nt = ",nt, "="))
      # print(guide_fun)
      # print("=============================\n")
      s_not_1_weights <- guide_fun/filter_guide_fun
      # print(paste0("s_not_1 weight for s=",s, ", nt = ",nt, "=", "\n"))
      # print(s_not_1_weights)
      # print("=============================\n")
      if (!(s==1 & nt!=0)){
        weights <- s_not_1_weights
        # print(paste0("weight for s=",s, "nt = ",nt, "=", "\n"))
        # print(weights)
        # print("=============================\n")
      }
      else {
        print("hi")
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        # print(paste0("x for s=",s, "nt = ",nt, "=", "for dmeasure = ", "\n"))
        # print(x_3d)
        # print("=============================\n")
        # print(paste0("y for s=",s, "nt = ",nt, "=", "for dmeasure = ", "\n"))
        # print(object@data[,nt,drop=FALSE])
        # print("=============================\n")
        weights <- tryCatch(
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=params,
            log=FALSE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        # print(paste0("dmeasure weights for s=",s, "nt = ",nt, "=", "\n"))
        # print(weights)
        # print("=============================\n")
        gnsi <- FALSE

        weights <- as.numeric(weights)*s_not_1_weights
        # print(paste0("weight for s=",s, "nt = ",nt, "after dmeasure and s_not_1 product", "=", "\n"))
        # print(weights)
        # print("=============================\n")
      }
      xx <- tryCatch(
        .Call('girf_computations',
              x=X,
              params=params,
              Np=Np[nt+1],
              predmean=pred.mean,
              predvar=pred.var,
              filtmean=filter.mean,
              trackancestry=filter.traj,
              doparRS=FALSE,
              weights=weights,
              gps=guide_fun,
              fsv=fcst_samp_var,
              tol=tol^(lookahead*length(object@units))
              ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      all.fail <- xx$fail
      eff.sample.size[nt+1, s] <- xx$ess
      loglik[nt+1, s] <- xx$loglik
      x <- xx$states
      filter_guide_fun <- xx$filterguides
      params <- xx$params[,1]
      fcst_samp_var <- xx$newfsv
      # print(paste0("loglik for s=",s, "nt = ",nt, "after resampling step", "=", "\n"))
      # print(xx$loglik)
      # print(paste0("states for s=",s, "nt = ",nt, "after resampling step", "=", "\n"))
      # print(xx$states)
      # print(paste0("filter guide functions for s=",s, "nt = ",nt, "after resampling step", "=", "\n"))
      # print(xx$filterguides)
      # print(paste0("params for s=",s, "nt = ",nt, "after resampling step", "=", "\n"))
      # print(xx$params[,1])
      # print(paste0("forecast sample variance for s=",s, "nt = ",nt, "after resampling step", "=", "\n"))
      # print(xx$newfsv)
      # print("=============================\n")
    }
  }
  return(loglik)
}

nonfinite_dmeasure_error <- function (time, lik, datvals, states, params) {
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",preserve.width="common")
  paste0(
    sQuote("dmeasure")," returns non-finite value.\n",
    "likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}
