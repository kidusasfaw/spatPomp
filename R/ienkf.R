setClass(
  "ienkfd_spatPomp",
  contains="enkfd_spatPomp",
  slots=c(Nenkf = 'integer',
          rw.sd = 'matrix',
          cooling.type = 'character',
          cooling.fraction.50 = 'numeric',
          traces = 'matrix'
      )
)


## ITERATIVE ENSEMBLE KALMAN FILTER (ienkf)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$
## Forecast: $Y_t=h(X_t)$
## Forecast mean: $N_t=\langle Y \rangle$.
## Forecast variance: $S_t=\langle\langle Y \rangle\rangle$
## State/forecast covariance: $W_t=\langle\langle X,Y\rangle\rangle$
## Kalman gain: $K_t = W_t\,S_t^{-1}$
## New observation: $y_t\in \mathbb{R}^{n\times 1}$
## Updated ensemble: $X^u_{t}=X_t + K_t\,(y_t - Y_t)$
## Filter mean: $m_t=\langle X^u_t \rangle = \frac{1}{q} \sum\limits_{i=1}^q x^{u_i}_t$

##' @name ienkf-spatPomp
##' @aliases ienkf,spatPomp-method
##' @rdname ienkf
##' @export
setMethod(
  "ienkf",
  signature=signature(data="spatPomp"),
  definition=function (data,
                       Nenkf = 1, rw.sd,
                       cooling.type = c("geometric", "hyperbolic"), cooling.fraction.50,
                       Np, tol = 0, max.fail = Inf,
                       ..., verbose = getOption("verbose", FALSE)) {
    tryCatch(
      ienkf.internal(
        data,
        Nenkf=Nenkf,
        rw.sd=rw.sd,
        cooling.type=match.arg(cooling.type),
        cooling.fraction.50=cooling.fraction.50,
        Np=Np,
        max.fail=max.fail,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp:::pStop("ienkf",conditionMessage(e))
    )
  }
)

ienkf.internal <- function (object, Nenkf, rw.sd,
                           cooling.type, cooling.fraction.50,
                           Np, max.fail = Inf,
                           ..., verbose,
                           .ndone = 0L, .indices = integer(0), .paramMatrix = NULL,
                           .gnsi = TRUE) {

  verbose <- as.logical(verbose)

  if (pomp:::undefined(object@rprocess) || pomp:::undefined(object@unit_emeasure) || pomp:::undefined(object@unit_vmeasure))
    pomp:::pStop_(paste(sQuote(c("rprocess","unit_emeasure","unit_vmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Nenkf) != 1 || !is.numeric(Nenkf) || !is.finite(Nenkf) || Nenkf < 1)
    pomp:::pStop_(sQuote("Nenkf")," must be a positive integer.")
  Nenkf <- as.integer(Nenkf)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {  ## if '.paramMatrix' is supplied, 'start' is ignored
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))

  if (is.null(Np)) {
    pomp:::pStop_(sQuote("Np")," must be specified.")
  }  else if (!is.numeric(Np)) {
    pomp:::pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }

  Np <- as.integer(Np)

  if (missing(rw.sd))
    pomp:::pStop_(sQuote("rw.sd")," must be specified!")
  rw.sd <- pomp:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    pomp:::pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pomp:::pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- pomp:::mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }

  traces <- array(dim=c(Nenkf+1,length(start)+1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Nenkf),
                                variable=c("loglik",names(start))))
  traces[1L,] <- c(NA,start)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)

  ## iterate the filtering
  for (n in seq_len(Nenkf)) {

    es <- ienkf.filter(
      object=object,
      params=paramMatrix,
      Np=Np,
      enkfiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      max.fail=max.fail,
      verbose=verbose,
      .indices=.indices,
      .gnsi=gnsi
    )

    gnsi <- FALSE
    paramMatrix <- es@paramMatrix
    traces[n+1,-1L] <- coef(es)
    traces[n,1L] <- es@loglik
    .indices <- es@indices

    if (verbose) cat("ienkf iteration",n,"of",Nenkf,"completed\n")

  }

  es@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
                              .gnsi=gnsi)

  new(
    "ienkfd_spatPomp",
    es,
    Nenkf=Nenkf,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

###################################################################
###################ienkf.filter()##################################
###################################################################
ienkf.filter <- function (object, params, Np, enkfiter, rw.sd, cooling.fn,
                 max.fail = Inf, verbose, .indices = integer(0),
                 .gnsi = TRUE) {

  verbose <- as.logical(verbose)
  gnsi <- as.logical(.gnsi)
  enkfiter <- as.integer(enkfiter)
  Np <- as.integer(Np)

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np)
    pomp:::pStop_(sQuote(".indices")," has improper length.")

  times <- tt <- time(object,t0=TRUE)
  t <- time(object)
  ntimes <- length(times)-1

  y <- obs(object)
  nobs <- nrow(y)

  loglik <- rep(NA,ntimes)

  for (nt in seq_len(ntimes)) {
    ## perturb parameters
    pmag <- cooling.fn(nt,enkfiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')
    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)

    ## get initial states
    if (nt == 1L) {
      X <- rinit(object,params=tparams)
      xnames <- rownames(X)
      tpnames <- rownames(tparams)
      pnames <- rownames(params)
    }

    ######################ENKF FROM HERE ON DOWN #################

    ## advance ensemble according to state process
    X <- rprocess(object,x0=X,t0=times[nt],times=times[nt+1],params=tparams,.gnsi=gnsi)

    # data
    yk <- y[,nt]

    # ensemble of forecasts
    Y <- tryCatch(
      .Call('do_theta_to_e',
            object=object,
            X=X,
            Np = as.integer(Np),
            times=times[nt+1],
            params=tparams,
            gnsi=gnsi),
      error = function (e) {
        stop("ep",conditionMessage(e),call.=FALSE) # nocov
      }
    )
    Y <- Y[,,1]

    # variance of artificial noise (i.e. R) computed using vmeasure
    meas_var <- tryCatch(
      .Call('do_theta_to_v',
            object=object,
            X=X,
            Np = Np,
            times=times[nt+1],
            params=tparams,
            gnsi=gnsi),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    dim(meas_var) <- c(length(spat_units(object)),  Np)
    R <- diag(rowMeans(meas_var))
    sqrtR <- tryCatch(
      t(chol(R)),                     # t(sqrtR)%*%sqrtR == R
      error = function (e) {
        pomp:::pStop_("degenerate ",sQuote("R"), "at time ", sQuote(nt), ": ",conditionMessage(e))
      }
    )

    # expand the state space
    XT <- rbind(X[,,1],params)
    pm <- rowMeans(XT) # prediction mean

    # forecast mean
    ym <- rowMeans(Y)

    # center prediction and forecast ensembles
    XT <- XT-pm
    Y <- Y-ym

    fv <- tcrossprod(Y)/(Np-1)+R  # forecast variance
    vyx <- tcrossprod(Y,XT)/(Np-1)   # forecast/state covariance

    svdS <- svd(fv,nv=0)            # singular value decomposition
    Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d) # transpose of Kalman gain
    Ek <- sqrtR%*%matrix(rnorm(n=nobs*Np),nobs,Np) # artificial noise
    resid <- y[,nt]-ym

    XT <- XT+pm+crossprod(Kt,resid-Y+Ek)
    params <- XT[pnames,,drop = FALSE]
    X <- XT[xnames,,drop = FALSE]
    loglik[nt] <- sum(dnorm(x=crossprod(svdS$u,resid),mean=0,sd=sqrt(svdS$d),log=TRUE))
    # print(rowMeans(partrans(object,params,dir="fromEst",.gnsi=gnsi)))

    ## compute mean at last timestep
    if (nt == ntimes) {
      coef(object,transform=TRUE) <- apply(params,1,mean)
    }
  }
  new("enkfd_spatPomp",
      object,
      Np=Np,
      cond.loglik=loglik,
      loglik=sum(loglik),
      indices=.indices,
      paramMatrix=params,
      unit_rmeasure = object@unit_rmeasure,
      unit_dmeasure = object@unit_dmeasure,
      unit_emeasure = object@unit_emeasure,
      unit_vmeasure = object@unit_vmeasure,
      unit_mmeasure = object@unit_mmeasure,
      units=object@units,
      unit_statenames=object@unit_statenames,
      obstypes = object@obstypes)
}



