##' Adapted Simulation Island Filter with Intermediate Resampling (ASIF-IR)
##'
##' An algorithm for estimating the likelihood of a spatiotemporal partially-observed
##' Markov process (SpatPOMP for short).
##' Running \code{asifir} causes the algorithm to run independent island jobs which
##' each carry out an adapted simulation using intermediate resampling.
##' Adapted simulation is an easier task than filtering, since particles in each island
##' remain close to each other. Intermediate resampling further assists against
##' the curse of dimensionality (COD) problem for importance sampling.
##' The adapted simulations are then weighted in a way that tries to avert COD by
##' making a weak coupling assumption to get an approximate filter distribution.
##' As a by-product, we also get an approximation to the likelihood of the data.
##'
##' @name asifir
##' @rdname asifir
##' @include spatPomp_class.R generics.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @inheritParams asif
##' @inheritParams girf
##' @inheritParams pomp::pfilter
##' @param object A \code{spatPomp} object.
##' @param Np The number of particles for the adapted simulations within each island.
##' @param islands The number of islands for the adapted simulations.
##' @examples
##' # Create a simulation of a BM using default parameter set
##' b <- bm(U=3, N=10)
##'
##' # Create a neighborhood function mapping a point in space-time to a list of ``neighboring points" in space-time
##' bm_nbhd <- function(object, time, unit) {
##'   nbhd_list = list()
##'   if(time > 1 && unit > 1) nbhd_list = c(nbhd_list, list(c(unit - 1, time - 1)))
##'   return(nbhd_list)
##' }
##' # Specify the expected measurement, theta.to.v and v.to.theta
##' girf.h <- function(state.vec, param.vec){
##'   # find index matching unit_statename
##'   ix<-grep('X',names(state.vec))
##'   state.vec[ix]
##' }
##' girf.theta.to.v <- function(meas.mean, param.vec){
##'   param.vec['tau']^2
##' }
##' girf.v.to.theta <- function(var, state.vec, param.vec){
##'   param.vec['tau'] <- sqrt(var)
##'   param.vec
##' }
##' # Run ASIFIR specified number of Monte Carlo islands and particles per island
##' asifird.b <- asifir(b,
##'                    islands = 50,
##'                    Np=20,
##'                    nbhd = bm_nbhd,
##'                    Ninter = length(spat_units(bm3)),
##'                    h = girf.h,
##'                    theta.to.v = girf.theta.to.v,
##'                    v.to.theta = girf.v.to.theta)
##' # Get the likelihood estimate from ASIFIR
##' logLik(asifird.b)
##'
##' # Compare with the likelihood estimate from Particle Filter
##' pfd.b <- pfilter(b, Np = 500)
##' logLik(pfd.b)
##' @return
##' Upon successful completion, \code{asifir} returns an object of class
##' \sQuote{asifird_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields a biased estimate of the log-likelihood of
##' the data under the model. }
##' }
##'
NULL

#setGeneric("asifir",function(object,...)standardGeneric("asifir"))

setClass(
  "asifird_spatPomp",
  contains="spatPomp",
  slots=c(
    Ninter="integer",
    Np="integer",
    islands="integer",
    tol="numeric",
    loglik="numeric",
    nbhd="function"
  ),
  prototype=prototype(
    Ninter=as.integer(NA),
    Np=as.integer(NA),
    islands=as.integer(NA),
    tol=as.double(NA),
    loglik=as.double(NA),
    nbhd=function(){}
  )
)
asifir.internal <- function (object, params, Np, nbhd,
      Ninter, tol, .gnsi = TRUE,...) {
  ep <- paste0("in ",sQuote("asifir"),": ")
  verbose <- FALSE
  if(missing(nbhd))
    stop(ep,sQuote("nbhd")," must be specified for the spatPomp object",call.=FALSE)
  object <- as(object,"spatPomp")
  pompLoad(object,verbose=verbose)
  gnsi <- as.logical(.gnsi)

  if (length(params)==0) stop(ep,sQuote("params")," must be specified",call.=FALSE)
  if (missing(tol)) stop(ep,sQuote("tol")," must be specified",call.=FALSE)

  times <- time(object,t0=TRUE)
  N <- length(times)-1
  U <- length(object@units)

  if (missing(Np)) stop(ep,sQuote("Np")," must be specified",call.=FALSE)
  if (is.function(Np)) stop(ep,"Functions for Np not supported by asifir",call.=FALSE)
  if (length(Np)!=1) stop(ep,"Np should be a length 1 vector",call.=FALSE)
  Np <- as.integer(Np)

  if (NCOL(params)>1) stop(ep,"does not accept matrix parameter input",call.=FALSE)

  coef(object) <- params
  paramnames <- names(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have names",call.=FALSE)

  ## ideally, we shouldn't need param_matrix in asifir
  param_matrix <- matrix(params,nrow=length(params),ncol=Np,
    dimnames=list(names(params),NULL))

  x_init <- rinit(object,params=params,nsim=1,.gnsi=gnsi) # Nx x 1 matrix
  statenames <- rownames(x_init)
  Nx <- nrow(x_init)
  xas <- as.numeric(x_init) # adapted simulation state vector

  loglik <- rep(NA,N)
  cond.densities <- array(data = numeric(0), dim=c(U,Np,N))
  dimnames(cond.densities) <- list(unit = 1:U, rep = 1:Np, time = 1:N)

  for (n in seq_len(N)) {
    ## assimilate observation n given filter at n-1
    ## note that times[n+1] is the time for observation n

    ## xg: Nx x Np x 1 matrix of guide simulations
    ## also used to calculate local prediction weights
    xf <- matrix(xas,nrow=Nx,ncol=Np,dimnames=list(statenames,NULL))
    xg <- tryCatch(
      rprocess(
        object,
        x0=xf,
        t0=times[n],
        times=times[n+1],
        params=params,
        .gnsi=gnsi
      ),
      error = function (e) stop(ep,"process simulation error: ",
        conditionMessage(e),call.=FALSE)
    )
    xg_2dim <- xg[,,1]
    # print(xg_2dim)
    # print("xg_2dim done")
    xg_with_rep <- do.call(cbind,replicate(Np, xg_2dim, simplify=FALSE))
    dim(xg_with_rep) <- c(dim(xg_with_rep)[1],dim(xg_with_rep)[2],1)
    dimnames(xg_with_rep) <- list(states = rownames(xg_2dim))
    xx <- tryCatch(
      .Call('do_fcst_samp_var',
            object=object,
            X=xg_with_rep,
            Np = as.integer(Np),
            times=times[n+1],
            params=params,
            gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    fcst_samp_var <- xx
    dim(fcst_samp_var) <- c(length(spat_units(object)), Np)

    ## determine the weights
    weights <- tryCatch(
      vec_dmeasure(
        object,
        y=object@data[,n,drop=FALSE],
        x=xg,
        times=times[n+1],
        params=param_matrix,
        log=FALSE,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    weights[weights < tol] <- tol
    cond.densities[,,n] <- weights[,,1]
    ## adapted simulation via intermediate resampling
    # tt has S+1 (or Ninter+1) entries
    tt <- seq(from=times[n],to=times[n+1],length.out=Ninter+1)
    log_gf <- rep(0,Np) # filtered guide function
    for (s in 1:Ninter){
      xp <- rprocess(object,x0=xf, t0 = tt[s], times= tt[s+1],
        params=params,.gnsi=gnsi) # an Nx by nreps by 1 array
      if(s < Ninter){
        skel <- pomp::flow(object, x0=xp[,,1], t0=tt[s+1],
          params=param_matrix, times = times[n + 1],...)
      } else {
        skel <- xp
      }

      # create measurement variance at skeleton matrix
      # meas_var_skel <- array(0, dim = c(U, Np))
      # for(u in 1:U){
      #   snames = paste0(object@unit_statenames,u)
      #   for(np in 1:Np){
      #     hskel <- h(state.vec=skel[snames,np,1],param.vec=params)
      #     meas_var_skel[u,np] <- theta.to.v(meas.mean=hskel,param.vec=params)
      #   }
      # }
      meas_var_skel <- tryCatch(
        .Call('do_theta_to_v',
              object=object,
              skel=skel,
              Np = as.integer(Np[1]),
              times=times[n+1],
              params=params,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      meas_var_skel <- meas_var_skel[,,1]
      fcst_var_upd <- fcst_samp_var*(times[n+1] - tt[s+1])/(times[n+1] - times[n])
      # print(xg)
      # print("done with xg")
      # print(meas_var_skel)
      # print("done with mvs")
      # print(fcst_samp_var)
      # print("done with fsv")
      # print(fcst_var_upd)
      # print("done with fvu")
      # return(1)
      # fcst_var_upd <- array(0, dim = c(U, Np))
      # for(u in 1:U){
      #   fcst_var_upd[u,] <- fcst_samp_var[u] *
      #     (times[n+1] - tt[s+1])/(times[n+1] - times[n])
      # }
      inflated_var <- meas_var_skel + fcst_var_upd
      dim(inflated_var) <- c(U, Np, 1)
      array.params <- array(params, dim = c(length(params), length(spat_units(object)), Np, 1), dimnames = list(params = names(params)))

      mmp <- tryCatch(
        .Call('do_v_to_theta',
              object=object,
              X=skel,
              vc=inflated_var,
              Np = as.integer(Np[1]),
              times=times[n+1],
              params=array.params,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      mom_match_param <- mmp[,,,1]

#       mom_match_param <- array(0, dim = c(length(params), U, Np),
#         dimnames = list(params = names(params), unit = NULL, J = NULL))
#       inflated_var <- meas_var_skel + fcst_var_upd
#       for(u in 1:U){
#         for(np in 1:Np){
#           snames = paste0(object@unit_statenames,u)
#           mom_match_param[,u,np] <- v.to.theta(var=inflated_var[u,np],
# 	    param.vec=params, state.vec=skel[snames,np,1])
#         }
#       }

      # U x Np x 1 matrix of skeleton prediction weights
      log_wp <- tryCatch(
        vec_dmeasure(
          object,
          y=object@data[,n,drop=FALSE],
          x=skel,
          times=times[n+1],
          params=mom_match_param,
          log=TRUE,
          .gnsi=gnsi
        ),
        error = function (e) stop(ep,"error in calculation of wp: ",
          conditionMessage(e),call.=FALSE)
      )
      log_gp <- apply(log_wp[,,1,drop=FALSE],2,sum)
      log_gp <- log_gp - max(log_gp)
      log_gp[log_gp < log(tol)] <- log(tol)
      weights <- exp(log_gp - log_gf)

      gnsi <- FALSE

      xx <- tryCatch(
          .Call('asifir_resample', xp, Np, weights, log_gp, tol),
          error = function (e) stop(ep,conditionMessage(e),call.=FALSE)
      )
      xf <- xx$states
      log_gf <- xx$filterguides
    } ## end of intermediate resampling loop s = 1:Ninter

    # resample down to one particle, making Np copies of, say, particle #1.
    xas <- xf[,1]

    if (verbose && (n%%5==0)) cat("asif timestep",n,"of",N,"finished\n")

  } ## end of main loop n = 1:N

  # compute locally combined pred. weights for each time, unit and particle
  #
  # matches asif.R except
  #   pp -> np
  #   Np is assumed scalar
  #   ntimes -> N , nt -> n
  #   nunits -> U , unit -> u
  loc.comb.pred.weights = array(data = numeric(0), dim=c(U,Np, N))
  wm.times.wp.avg = array(data = numeric(0), dim = c(U, N))
  wp.avg = array(data = numeric(0), dim = c(U, N))
  for (n in seq_len(N)){
      for (u in seq_len(U)){
          full_nbhd <- nbhd(object, n, u)
          prod_cond_dens_nt  <- rep(1, Np)
          prod_cond_dens_not_nt <- matrix(1, Np, n-1)
          for (neighbor in full_nbhd){
              neighbor_u <- neighbor[1]
              neighbor_n <- neighbor[2]
              if (neighbor_n == n)
                  prod_cond_dens_nt  <- prod_cond_dens_nt * cond.densities[neighbor_u, ,neighbor_n]
              else
                  prod_cond_dens_not_nt[, neighbor_n] <- prod_cond_dens_not_nt[, neighbor_n] * cond.densities[neighbor_u, ,neighbor_n]
          }
          loc.comb.pred.weights[u,,n]  <- prod(apply(prod_cond_dens_not_nt, 2, mean))*prod_cond_dens_nt
      }
  }
  wm.times.wp.avg = apply(loc.comb.pred.weights * cond.densities, c(1,3), FUN = mean)
  wp.avg = apply(loc.comb.pred.weights, c(1,3), FUN = mean)

  pompUnload(object,verbose=verbose)
  new(
    "island.spatPomp",
    wm.times.wp.avg = wm.times.wp.avg,
    wp.avg = wp.avg,
    Np=as.integer(Np),
    tol=tol
  )


}
##' @name asifir-spatPomp
##' @aliases asifir,spatPomp-method
##' @rdname asifir
##' @export
setMethod(
  "asifir",
  signature=signature(object="spatPomp"),
  function (object, params, Np, islands, nbhd,
            Ninter, tol = (1e-18), ...) {
  if (missing(params)) params <- coef(object)
    # set.seed(396658101,kind="L'Ecuyer")
  # begin single-core
  # single_island_output <- spatPomp:::asifir.internal(
  #   object=object,
  #   params=params,
  #   Np=Np,
  #   nbhd=nbhd,
  #   Ninter=Ninter,
  #   tol=tol,
  #   ...
  # )
  # return(single_island_output)
  # end single-core
  mult_island_output <- foreach::foreach(i=1:islands,
     .packages=c("pomp","spatPomp"),
     .options.multicore=list(set.seed=TRUE)) %dopar% spatPomp:::asifir.internal(
       object=object,
       params=params,
       Np=Np,
       nbhd=nbhd,
       Ninter=Ninter,
       tol=tol,
       ...
     )
   # compute sum (over all islands) of w_{d,n,i}^{P} for each (d,n)
   N <- length(object@times)
   U <- length(object@units)
   island_mp_sums = array(data = numeric(0), dim = c(U,N))
   island_p_sums = array(data = numeric(0), dim = c(U, N))
   cond.loglik = array(data = numeric(0), dim=c(U, N))
   for (u in seq_len(U)){
     for (n in seq_len(N)){
       mp_sum = 0
       p_sum = 0
       for (k in seq_len(islands)){
         mp_sum = mp_sum + mult_island_output[[k]]@wm.times.wp.avg[u,n]
         p_sum = p_sum + mult_island_output[[k]]@wp.avg[u,n]
       }
       cond.loglik[u,n] = log(mp_sum) - log(p_sum)
     }
   }
   new(
      "asifird_spatPomp",
      object,
      Np=as.integer(Np),
      tol=tol,
      loglik=sum(cond.loglik),
      Ninter = as.integer(Ninter),
      islands = as.integer(islands),
      nbhd = nbhd
      )
  }
)

setMethod(
  "asifir",
  signature=signature(object="asifird_spatPomp"),
  function (object, params, Np, islands, nbhd,
            Ninter, tol = (1e-18), ...) {
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    if (missing(Ninter)) Ninter <- object@Ninter
    if (missing(params)) params <- coef(object)
    if (missing(islands)) islands <- object@islands
    if (missing(nbhd)) nbhd <- object@nbhd

    asifir(as(object,"spatPomp"),
         Np=Np,
         Ninter=Ninter,
         islands=islands,
         nbhd=nbhd,
         params = params,
         tol=tol,
         ...)
  }
)

