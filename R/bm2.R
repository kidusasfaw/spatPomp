#' Brownian motion spatPomp generator with shared or unit-specific parameters
#'
#' An extension of \code{bm} allowing for shared or unit-specific parameters.
#' Generate a class \sQuote{spatPomp} object representing a \code{U}-dimensional
#' Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @name bm2
#' @rdname bm2
#' @author Edward L. Ionides
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @param delta_t Process simulations are performed every \code{delta_t} time units
#' whereas observations occur every one time unit
#' @param shared_names identifies parameters that have common shared value for all units, which by default is all parameters.
#' @param unit_specific_names determines which parameters take a different value
#' for each unit. Cannot be specified if shared_names is specified.
#' each unit. Other parameters are considered shared between all units.
#' @param unit_params parameter values used to build the object, copied across 
#' each unit for unit-specific parameters
#' @importFrom utils data
#' @return An object of class \sQuote{spatPomp} representing a simulation from a
#' \code{U}-dimensional Brownian motion
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' b <- bm2(U=4, N=20,shared_names="rho",unit_specific_names=c("sigma","tau"))
#' # See all the model specifications of the object
#' spy(b)
#' # Examples of methodologies applied to this model
#' # are provided in the tests directory
#' }
#' @export

bm2 <- function(U=5,N=100,delta_t=0.1,
  unit_specific_names="rho",
  shared_names=NULL,
  unit_params =c(rho=0.4,sigma=1,tau=1,X_0=0)){
  dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
  dmat <- matrix(0,U,U)
  for(u in 1:U) {
    for(v in 1:U) {
      dmat[u,v] <- dist(u,v)
    }
  }
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  dist_C_rows <- apply(dmat,1,to_C_array)
  dist_C_array <- to_C_array(dist_C_rows)
  dist_C <- paste0("const int dist[",U,"][",U,"] = ",dist_C_array,"; ")

  obs_names <- paste0("U",1:U)
  bm2_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)
  bm2_unitnames <- unique(bm2_data[["unit"]])
  bm2_unitnames_level <- paste("U",sort(as.numeric(stringr::str_remove(bm2_unitnames, "U"))),sep='')

  bm2_unit_statenames <- "X"
  bm2_unit_IVPnames <- "X_0"
  bm2_unit_RPnames <- c("rho","sigma","tau")
  bm2_unit_paramnames <- c(bm2_unit_RPnames,bm2_unit_IVPnames)


  if(!missing(shared_names)) {
    if(!missing(unit_specific_names)) {
      stop ("both shared_names and unit_specific names cannot be given to bm2")
    } else
      unit_specific_names <- bm2_unit_paramnames[!bm2_unit_paramnames%in%shared_names]
  }
  
  if(missing(shared_names)) {
    if(missing(unit_specific_names)) {
      shared_names <- bm2_unit_paramnames
      unit_specific_names <- NULL
    } else shared_names <- bm2_unit_paramnames[!bm2_unit_paramnames%in%unit_specific_names]
  }

  set_unit_specific <- Csnippet(paste0("const int ", unit_specific_names,
    "_unit = 1;\n", collapse=" "))
  set_shared <- Csnippet(paste0("const int ", shared_names,
    "_unit = 0;\n", collapse=" "))

  bm2_globals <- Csnippet(
    paste(dist_C,set_unit_specific, set_shared, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  bm2_paramnames <- c(
    if(length(shared_names)>0){
      paste0(shared_names, "1")
    },
    if(length(unit_specific_names)>0){
      paste0(rep(unit_specific_names, each=U), 1:U)
    }
  )

  bm2_rprocess <- spatPomp_Csnippet("
    const double *rho=&rho1;
    const double *sigma=&sigma1;
    double dW[U];
    int u,v;

    for (u = 0 ; u < U ; u++) {
      dW[u] = rnorm(0,sigma[u*sigma_unit]*sqrt(dt));
    }
    for (u = 0 ; u < U ; u++) {
      for (v=0; v < U ; v++) {
        X[u] += dW[v]*pow(rho[u*rho_unit],dist[u][v]);
      }
    }
  ", unit_statenames = c("X"))

  bm2_skel <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_vfnames = c("X"),
    code = "
      for (int u = 0 ; u < U ; u++) {
        DX[u] = 0;
      }
    "
  )

  bm2_rinit <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    code = "
      const double *X_0 = &X_01;
      for (int u = 0; u < U; u++) {
        X[u]=X_0[u*X_0_unit];
      }
    "
  )

  bm2_dmeasure <- Csnippet("
    const double *tau = &tau1;
    const double *X = &X1;
    const double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    lik=0;
    for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau[u*tau_unit],1);
    if(!give_log) lik = exp(lik) + tol;
  ")

  bm2_eunit_measure <- Csnippet("
    ey = X;
  ")

  bm2_vunit_measure <- Csnippet("
    const double *tau = &tau1;   
    vc = tau[u*tau_unit]*tau[u*tau_unit];
  ")

  bm2_rmeasure <- Csnippet("
    const double *tau = &tau1;
    const double *X = &X1;
    double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau[u*tau_unit]+tol);
  ")

  bm2_dunit_measure <- Csnippet("
    const double *tau = &tau1;
    lik = dnorm(Y,X,tau[u*tau_unit],1);
    if(!give_log) lik = exp(lik);
  ")

  bm2_runit_measure <- Csnippet("
    const double *tau = &tau1;
    Y = rnorm(X,tau[u*tau_unit]);
  ")

log_unit_names <- c("sigma", "tau")
logit_unit_names <- c("rho")
log_names <- unlist(lapply(log_unit_names,
  function(x,y,U){if(x%in%y)paste0(x,"1") else paste0(x,1:U)},
  y=shared_names,U=U))
logit_names <- unlist(lapply(logit_unit_names,
  function(x,y,U){if(x%in%y)paste0(x,"1") else paste0(x,1:U)},
  y=shared_names,U=U))

bm2_partrans <- parameter_trans(log=log_names,logit=logit_names)


  bm2_spatPomp <- spatPomp(bm2_data %>% dplyr::arrange(time, factor(.data$unit, levels = bm2_unitnames_level)),
                 times="time",
                 t0=0,
                 units="unit",
                 unit_statenames = bm2_unit_statenames,
                 rprocess=euler(bm2_rprocess,delta.t = delta_t),
                 skeleton=vectorfield(bm2_skel),
                 paramnames=bm2_paramnames,
                 globals=bm2_globals,
                 rmeasure=bm2_rmeasure,
                 dmeasure=bm2_dmeasure,
                 eunit_measure=bm2_eunit_measure,
                 vunit_measure=bm2_vunit_measure,
                 dunit_measure=bm2_dunit_measure,
                 runit_measure=bm2_runit_measure,
                 partrans = bm2_partrans,
                 rinit=bm2_rinit
    )

  ## We need a parameter vector.
  bm2_params <- rep(0,length=length(bm2_paramnames))
  names(bm2_params) <- bm2_paramnames
  
  for(p in unit_specific_names) {
    bm2_params[paste0(p,1:U)] <- unit_params[p]
  }
  for(p in shared_names) {
    bm2_params[paste0(p,1)] <- unit_params[p]
  }
  simulate(bm2_spatPomp,params=bm2_params)
}

