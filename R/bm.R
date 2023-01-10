#' Brownian motion spatPomp simulator
#'
#' Generate a class \sQuote{spatPomp} object representing a \code{U}-dimensional
#' Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @name bm
#' @rdname bm
#' @author Edward L. Ionides
#' @family spatPomp model generators
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @param delta_t Process simulations are performed every \code{delta_t} time units
#' whereas observations occur every one time unit
#' @importFrom utils data
#' @return An object of class \sQuote{spatPomp} representing a simulation from a \code{U}-dimensional
#' Brownian motion
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' b <- bm(U=4, N=20)
#' # See all the model specifications of the object
#' spy(b)
#' # Examples of methodologies applied to this model
#' # are provided in the tests directory
#' }
#' @export

bm <- function(U=5,N=100,delta_t=0.1){
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
  bm_globals <- Csnippet(paste0(dist_C))


  obs_names <- paste0("U",1:U)
  bm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)
  bm_unitnames <- unique(bm_data[["unit"]])
  bm_unitnames_level <- paste("U",sort(as.numeric(stringr::str_remove(bm_unitnames, "U"))),sep='')

  bm_unit_statenames <- c("X")
  bm_statenames <- paste0(bm_unit_statenames,1:U)

  bm_IVPnames <- paste0(bm_statenames,"_0")
  bm_RPnames <- c("rho","sigma","tau")
  bm_paramnames <- c(bm_RPnames,bm_IVPnames)

  bm_rprocess <- spatPomp_Csnippet("
    double dW[U];
    double pow_rho[U];
    int u,v;

    pow_rho[0] = 1;
    for (u=1 ; u < U ; u++) {
      pow_rho[u] = pow_rho[u-1]*rho;
    }

    for (u = 0 ; u < U ; u++) {
      dW[u] = rnorm(0,sigma*sqrt(dt));
    }
    for (u = 0 ; u < U ; u++) {
      for (v=0; v < U ; v++) {
        X[u] += dW[v]*pow_rho[dist[u][v]];
      }
    }
  ", unit_statenames = c("X"))

  bm_skel <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_vfnames = c("X"),
    code = "
      for (int u = 0 ; u < U ; u++) {
        DX[u] = 0;
      }
    "
  )

  bm_rinit <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_ivpnames = c("X"),
    code = "
      for (int u = 0; u < U; u++) {
        X[u]=X_0[u];
      }
    "
  )

  bm_dmeasure <- Csnippet("
    const double *X = &X1;
    const double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    lik=0;
    for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau,1);
    if(!give_log) lik = exp(lik) + tol;
  ")

  bm_eunit_measure <- Csnippet("
    ey = X;
  ")

  bm_munit_measure <- Csnippet("
    M_tau = sqrt(vc);
  ")

  bm_vunit_measure <- Csnippet("
    vc = tau*tau;
  ")

  bm_rmeasure <- Csnippet("
    const double *X = &X1;
    double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau+tol);
  ")

  bm_dunit_measure <- Csnippet("
    lik = dnorm(Y,X,tau,1);
    if(!give_log) lik = exp(lik);
  ")

  bm_runit_measure <- Csnippet("
    Y = rnorm(X,tau);
  ")

  bm_spatPomp <- spatPomp(bm_data %>% dplyr::arrange(time, factor(.data$unit, levels = bm_unitnames_level)),
                 times="time",
                 t0=0,
                 units="unit",
                 unit_statenames = bm_unit_statenames,
                 rprocess=euler(bm_rprocess,delta.t = delta_t),
                 skeleton=vectorfield(bm_skel),
                 paramnames=bm_paramnames,
                 globals=bm_globals,
                 rmeasure=bm_rmeasure,
                 dmeasure=bm_dmeasure,
                 eunit_measure=bm_eunit_measure,
                 munit_measure=bm_munit_measure,
                 vunit_measure=bm_vunit_measure,
                 dunit_measure=bm_dunit_measure,
                 runit_measure=bm_runit_measure,
                 partrans = parameter_trans(log = c("sigma", "tau"), logit = c("rho")),
                 rinit=bm_rinit
    )

  ## We need a parameter vector. For now, we initialize the process at zero.
  test_ivps <- rep(0,U)
  names(test_ivps) <- bm_IVPnames
  test_params <- c(rho=0.4, sigma=1, tau=1, test_ivps)
  simulate(bm_spatPomp,params=test_params)
}
