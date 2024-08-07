## ----packages,include=F,echo=F,cache=F----------------------------------------
library("spatPomp")
library("gridExtra")
library("grid")
library("ggplot2")
library("xtable")
library("tidyverse")
library("cowplot")
library("knitr")
library("kableExtra")






## ----knitr-opts-for-purl,include=F,cache=F,purl=T-----------------------------
opts_chunk$set(
    warning=FALSE,
    message=FALSE
)


## ----setup,cache=F,include=F--------------------------------------------------
myround<- function (x, digits = 1) {
  # adapted from the broman package
  # solves the bug that round() kills significant trailing zeros
  if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  if (digits < 1) {
    as.character(round(x,digits))
  } else {
    tmp <- sprintf(paste("%.", digits, "f", sep = ""), x)
    zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
    tmp[tmp == paste0("-", zero)] <- zero
    tmp
  }
}
mysignif <- function (x, digits = 1) {
  myround(x, digits - ceiling(log10(abs(x))))
}


## ----set-seed,eval=T,echo=F---------------------------------------------------
bm_sim_U <- 10
bm_sim_N <- 20
set.seed(76)

## ----example_bm4,eval=T,echo=F------------------------------------------------
bm10 <- bm(U = 10, N = 20)


## ----bm10-simplot,eval=T,echo=F,fig.height=2.5,fig.width=6.5,out.width="100%",fig.cap=paste0(r"{Result of executing \code{plot(bm10)}, where \code{bm10} is the \class{spatPomp} object representing a simulation from a }",bm_sim_U,r"{-dimensional correlated Brownian motions model with }",bm_sim_N,r"{ observations that are one unit time apart (see text).}")----
plot(bm10,nrow=2)+
  theme(plot.margin=unit(c(0,0.1,0,0),"cm"))




## ----abf_parallel,eval=T,echo=T-----------------------------------------------
library("doParallel")
registerDoParallel(detectCores()) 

## ----cores,eval=T,echo=F------------------------------------------------------
cores <- detectCores()


## ----example_nbhd-------------------------------------------------------------
example_nbhd <- function(object, unit, time){
  nbhd_list = list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}


## ----bm10-params,echo=T-------------------------------------------------------
coef(bm10)




## ----bm_set_i,cache=F,echo=F,eval=T-------------------------------------------
i <- 2


## ----bm-settings,cache=FALSE,echo=F-------------------------------------------
bm_nbhd <- function(object, time, unit) {
  nbhd_list = list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
  return(nbhd_list)
}

bm_files_dir <- paste0("bm_",i,"/")
if(!dir.exists(bm_files_dir)) dir.create(bm_files_dir)
bm_U <- switch(i,c(2,3),c(40,30,20,10))
bm_N <- switch(i,2,15)
bm_lik_reps=switch(i,2,5)
bm_jobs <- rep(seq_along(bm_U),each=bm_lik_reps)


## ----bm_sim,echo=F------------------------------------------------------------
bm_list <- foreach(U=iter(bm_U))%do% bm(U=U,N=bm_N)


## ----bm_kf,cache=F,echo=F-----------------------------------------------------
bm_kalman_loglik <- function(spo,theta){
  U <- length(unit_names(spo))
  rows <- matrix(seq(U),nrow=U,ncol=U)
  cols <- matrix(seq(U),nrow=U,ncol=U,byrow=TRUE)
  bm_dmat <- pmin(abs(rows-cols), abs(rows-cols+U), abs(rows-cols-U))
  rootQ <- theta["rho"]^bm_dmat * theta["sigma"]
  kalmanFilter(spo,
    X0=rinit(spo),
    A= diag(U),
    Q= rootQ %*% rootQ,
    C=diag(U),
    R=diag(theta["tau"]^2, nrow=U)
  )$logLik
}
bm_kf <- bake(file=paste0(bm_files_dir,"bm_kf.rds"),{
  foreach(bm_obj=iter(bm_list),.combine=c) %do% {
    bm_kalman_loglik(bm_obj,coef(bm_obj))
  }
})


## ----bm_girf,cache=F,echo=F---------------------------------------------------
bm_girf <- bake(file=paste0(bm_files_dir,"bm_girf.rds"),{
  girf_list <- foreach(j=iter(bm_jobs)) %dopar% girf(
    bm_list[[j]],
    Np=switch(i,10,500),
    Ninter=switch(i,2,5),
    lookahead=1,
    Nguide=switch(i,5,50),
    tol=1e-300)
  girf_loglik <-  sapply(girf_list,logLik)
  attr(girf_loglik,"Np") <- girf_list[[1]]@Np
  attr(girf_loglik,"Nguide") <- girf_list[[1]]@Nguide
  attr(girf_loglik,"Ninter") <- girf_list[[1]]@Ninter
  attr(girf_loglik,"cores") <- cores
  girf_loglik
})
compute_cores <- attr(bm_girf,"cores")


## ----bm_abf,cache=F,echo=F----------------------------------------------------
bm_abf <- bake(file=paste0(bm_files_dir,"bm_abf.rds"),{
  abf_list <- foreach(j=iter(bm_jobs)) %do% abf(
    bm_list[[j]],
    Nrep = switch(i,3,500),
    Np = switch(i,10,100),
    nbhd = bm_nbhd,
    tol=1e-300)
  abf_loglik <- sapply(abf_list,logLik)
  attr(abf_loglik,"Np") <- abf_list[[1]]@Np
  attr(abf_loglik,"Nrep") <- abf_list[[1]]@Nrep
  abf_loglik
})


## ----bm_pfilter,cache=F,echo=F------------------------------------------------
bm_pfilter <- bake(file=paste0(bm_files_dir,"bm_pfilter.rds"),{
  pfilter_list <- foreach(j=iter(bm_jobs)) %dopar% pfilter(
    bm_list[[j]],
    Np=switch(i,50,2000))
  pfilter_loglik <-  sapply(pfilter_list,logLik)
  attr(pfilter_loglik,"Np") <- pfilter_list[[1]]@Np[1]
  pfilter_loglik
})


## ----bm_enkf,cache=F,echo=F---------------------------------------------------
bm_enkf <- bake(file=paste0(bm_files_dir,"bm_enkf.rds"),{
  enkf_list <- foreach(j=iter(bm_jobs)) %dopar% enkf(
    bm_list[[j]], Np=switch(i,50,2000))
  enkf_loglik <-  sapply(enkf_list,logLik)
  attr(enkf_loglik,"Np") <- enkf_list[[1]]@Np
  enkf_loglik
})


## ----bm_bpf,cache=F,echo=F----------------------------------------------------
bm_bpf <- bake(file=paste0(bm_files_dir,"bm_bpf.rds"),{
  bpf_list <- foreach(j=iter(bm_jobs)) %dopar% bpfilter(
    bm_list[[j]],
    Np=switch(i,50,2000),
    block_size = switch(i,1,2))
  bpf_loglik <-  sapply(bpf_list,logLik)
  attr(bpf_loglik,"Np") <- bpf_list[[1]]@Np[1]
  attr(bpf_loglik,"block_size") <- max(sapply(bpf_list[[1]]@block_list,length))
  bpf_loglik  
})


## ----bm_results,echo=F--------------------------------------------------------
bm_methods <- c("ABF", "Particle Filter", "GIRF", "EnKF", "BPF", "KF")
bm_results <- data.frame(
  method=rep(bm_methods,each=length(bm_jobs)),
  U=rep(rep(bm_U,each=bm_lik_reps),6),
  logLik=c(bm_abf,bm_pfilter,bm_girf,bm_enkf,bm_bpf,rep(bm_kf,each=bm_lik_reps))
)
bm_results$dodge_units <- bm_results$U + rep(seq(from=-1,to=1,length=length(bm_methods)),each=bm_lik_reps*length(bm_U))
bm_results$logLik_per_unit <- bm_results$logLik/bm_results$U
bm_results$logLik_per_obs <- bm_results$logLik_per_unit/bm_N
bm_max <- max(bm_results$logLik_per_obs)


## ----bm_lik_plot, echo=F, eval = T, fig.height=3,fig.width=5.5,fig.align='center',out.width="75%",fig.cap = "Log-likelihood estimates for 5 replications of ABF, BPF, EnKF, GIRF and particle filter on correlated Brownian motions of various dimensions. The Kalman filter (KF) provides the exact likelihood in this case."----
ggplot(bm_results,mapping = aes(x=dodge_units, y=logLik_per_obs)) +
  geom_point(mapping = aes(color=method))+
  scale_color_manual(values=c("Particle Filter" = "#0072B2", "ABF" = "#E69F00", "GIRF" = "#009E73", "EnKF" = "#CC79A7", "BPF" = "#D55E00", "KF" = "#000000")) +
  scale_x_continuous(breaks=bm_U,
        labels=bm_U) +
  labs(x = "Number of spatial units", y = "Log-likelihood per unit per time")+
  coord_cartesian(ylim=c(bm_max-0.8,bm_max))+
  theme(legend.key.width = unit(1,"cm"))


## ----bm_table_df,echo=F,eval=T,cache=F----------------------------------------
bm_table_df <- data.frame(Method=c('Particle Filter', 'ABF', 'GIRF', 'EnKF', 'BPF'),
  Runtime=c(mysignif(attr(bm_pfilter, "system.time")["elapsed"]/60*cores,2),
            mysignif(attr(bm_abf, "system.time")["elapsed"]/60*cores,2),
            mysignif(attr(bm_girf, "system.time")["elapsed"]/60*cores,2),
            mysignif(attr(bm_enkf, "system.time")["elapsed"]/60*cores,2),
            mysignif(attr(bm_bpf, "system.time")["elapsed"]/60*cores,2)),
  Particles = c(attr(bm_pfilter,"Np"),
                attr(bm_abf,"Np"),
		attr(bm_girf,"Np"),
                attr(bm_enkf,"Np"),
                attr(bm_bpf,"Np")),				
  Replicates = c('-',
                 attr(bm_abf,"Nrep"),
                 '-',
                 '-',
                 '-'),
  Guides = c('-',
             '-',
             attr(bm_girf,"Nguide"),
             '-',
             '-'),
  Lookaheads = c('-',
                 '-',
                 1,
                 '-',
                 '-')
  )
bm_table_df %>% kable(booktabs=T,escape=F,
  align="lccccc",
  caption = 'Comparison of computational resources of the filtering algorithms',
  col.names = linebreak(c("Method \n { }",
                          "Resources \n (core-minutes)",
                          "Particles \n (per replicate)",
                          "Replicates \n { }",
                          "Guide \n particles",
                          "Lookahead \n { }")))


## ----bm_start_params,echo=T,eval=T--------------------------------------------
start_params <- c(rho = 0.8, sigma = 0.4, tau = 0.2,
  X1_0 = 0, X2_0 = 0, X3_0 = 0, X4_0 = 0, X5_0 = 0,
  X6_0 = 0, X7_0 = 0, X8_0 = 0, X9_0 = 0, X10_0 = 0)



## ----bm_igirf_eval, eval=T, echo=F--------------------------------------------
set.seed(5282122)
ig1 <- bake(file=paste0(bm_files_dir,"bm_igirf.rds"),{
i <- 2
ig1 <- igirf(
  bm10,
  params=start_params,
  Ngirf=switch(i,2,50),
  Np=switch(i,10,1000),
  Ninter=switch(i,2,5),
  lookahead=1,
  Nguide=switch(i,5,50),
  rw.sd=rw_sd(rho=0.02,sigma=0.02,tau=0.02),
  cooling.type = "geometric",
  cooling.fraction.50=0.5
)
attr(ig1,"cores") <- cores
ig1
})


## ----bm_igirf_enkf, eval = T, echo = F----------------------------------------
ig1_loglik <- bm_kalman_loglik(bm10,coef(ig1))

ig1_enkf <- logLik(enkf(bm10, params=coef(ig1), Np=1000)) 


## ----bm_mle,echo=F,eval=T-----------------------------------------------------
bm10_negloglik <- function(cf) -bm_kalman_loglik(bm10,cf)
mle <- optim(coef(bm10), bm10_negloglik)
kf_ml <- -mle$value



## ----bm_igirf_convergence2, echo = F, eval = T, fig.align = "center", out.width="65%",fig.width=6, fig.height=3.5, fig.cap=paste0(r"{The output of the \code{plot()} method on the object of \class{igirfd\_spatPomp} that encodes our model for correlated Brownian motions produces convergence traces for $\rho$, $\sigma$ and $\tau$, and the corresponding log-likelihoods. Over }", ig1@Ngirf, r"{ iterations \code{igirf()} has allowed us to get within a neighborhood of the maximum likelihood.}")----
plot(ig1, params = c("rho", "sigma", "tau"), ncol = 2) +
 theme(axis.text.x = element_text(size = 12),
       axis.text.y = element_text(size = 12),
       axis.title.x = element_text(size = 14),
       axis.title.y = element_text(size = 14),
       strip.text = element_text(size = 16))




## ----bm_rho_prof_bounds,echo=T------------------------------------------------
theta_lo_trans <- partrans(bm10,coef(bm10),dir="toEst") - log(2)
theta_hi_trans <- partrans(bm10,coef(bm10),dir="toEst") + log(2)
profile_design(
  rho=seq(from=0.2,to=0.6,length=10),
  lower=partrans(bm10,theta_lo_trans,dir="fromEst"),
  upper=partrans(bm10,theta_hi_trans,dir="fromEst"),
  nprof=switch(i,2,10)
) -> pd




## ----bm_rho_prof_eval,eval=T,cache=F,echo=F-----------------------------------
rho_prof <- bake(file = paste0(bm_files_dir,"rho-profile.rds"),{
foreach (p=iter(pd,"row"),.combine=dplyr::bind_rows) %dopar% {
  library(spatPomp)
  ig2 <- igirf(ig1,params=p,rw.sd=rw_sd(sigma=0.02,tau=0.02))
  ef <- replicate(switch(i,2,10),enkf(ig2,Np=switch(i,50,2000)))
  ll <- sapply(ef,logLik)
  ll <- logmeanexp(ll,se=TRUE)
  data.frame(as.list(coef(ig2)),loglik=ll[1],loglik.se=ll[2])
} -> rho_prof
attr(rho_prof,"cores") <- cores
rho_prof
})


## ----bm_rho_prof_mcap---------------------------------------------------------
rho_mcap <- mcap(rho_prof[,"loglik"],parameter=rho_prof[,"rho"])
rho_mcap$ci


## ----measles_cases_covar,eval=TRUE,echo=FALSE,cache=F-------------------------
measles_U <- 6
birth_lag <- 4*26  # delay until births hit susceptibles, in biweeks
data(measlesUK)
measlesUK$city<-as.character(measlesUK$city)
cities <- unique(measlesUK$city)[1:measles_U]
measles_cases <- measlesUK[measlesUK$city %in% cities,c("year","city","cases")]
measles_cases <- measles_cases[measles_cases$year>1949.99,]
measles_covar <- measlesUK[measlesUK$city %in% cities,c("year","city","pop","births")]
u <- split(measles_covar$births,measles_covar$city)
v <- sapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
measles_covar$lag_birthrate <- as.vector(v[,cities])*26
measles_covar$births<- NULL
measles_covar$P <- measles_covar$pop
measles_covar$pop <- NULL
measles_covarnames <- paste0(rep(c("P","lag_birthrate"),each=measles_U),1:measles_U)


## ----measles_print_cases,echo = F,cache=F-------------------------------------
print(head(measles_cases %>% arrange(year),8), row.names = FALSE)


## ----measles_first_spatpomp,eval=TRUE,cache=F---------------------------------
measles6 <- spatPomp(
  data=measles_cases,
  units='city',
  times='year',
  t0=min(measles_cases$year)-1/26
)


## ----measles_print_covar,echo=F,cache=F---------------------------------------
print(head(measles_covar %>% arrange(year) %>% filter(year >= 1950.000),3), row.names = FALSE)


## ----measles_globals,eval=T,cache=F-------------------------------------------
measles_globals <- spatPomp_Csnippet("
  const double V[6][6] = {
  {0,2.42,0.950,0.919,0.659,0.786},
  {2.42,0,0.731,0.722,0.412,0.590},
  {0.950,0.731,0,1.229,0.415,0.432},
  {0.919,0.722,1.229,0,0.638,0.708},
  {0.659,0.412,0.415,0.638,0,0.593},
  {0.786,0.590,0.432,0.708,0.593,0}
  }; 
")


## ----measles_rinit,echo=T, eval=T,cache=F-------------------------------------
measles_rinit <- spatPomp_Csnippet(
  unit_statenames = c('S','E','I','C'),
  unit_ivpnames = c('S','E','I'),
  unit_covarnames = c('P'),
  code = "
    for (int u=0; u<U; u++) {
      S[u] = round(P[u]*S_0[u]);
      E[u] = round(P[u]*E_0[u]);
      I[u] = round(P[u]*I_0[u]);
      C[u] = 0;
    }
  "
)


## ----measles_rprocess,echo=T, eval=T,cache=F----------------------------------
measles_rprocess <- spatPomp_Csnippet(
  unit_statenames = c('S','E','I','C'),
  unit_covarnames = c('P','lag_birthrate'),
  code = "
    double beta, seas, Ifrac, mu[7], dN[7];
    int u, v;
    int BS=0, SE=1, SD=2, EI=3, ED=4, IR=5, ID=6;

    beta = R0*(muIR+muD);
    t = (t-floor(t))*365.25;
    seas = (t>=7&&t<=100)||(t>=115&&t<=199)||(t>=252&&t<=300)||(t>=308&&t<=356)
      ? 1.0 + A * 0.2411/0.7589 : 1.0 - A;

    for (u = 0 ; u < U ; u++) {
      Ifrac = I[u]/P[u];
      for (v=0; v < U ; v++) if(v != u)
        Ifrac += g * V[u][v]/P[u] * (I[v]/P[v] - I[u]/P[u]);

      mu[BS] = lag_birthrate[u];   
      mu[SE] = beta*seas*Ifrac*rgammawn(sigmaSE,dt)/dt; 
      mu[SD] = muD;                
      mu[EI] = muEI;             
      mu[ED] = muD;  
      mu[IR] = muIR; 
      mu[ID] = muD;  

      dN[BS] = rpois(mu[BS]*dt);
      reulermultinom(2,S[u],&mu[SE],dt,&dN[SE]);
      reulermultinom(2,E[u],&mu[EI],dt,&dN[EI]);
      reulermultinom(2,I[u],&mu[IR],dt,&dN[IR]);

      S[u] += dN[BS] - dN[SE] - dN[SD];
      E[u] += dN[SE] - dN[EI] - dN[ED];
      I[u] += dN[EI] - dN[IR] - dN[ID];
      C[u] += dN[EI];           
    }
  "
)


## ----measles_dunit_measure,echo = TRUE, eval = TRUE,cache=F-------------------
measles_dunit_measure <- spatPomp_Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  lik = dnorm(cases,m,sqrt(v),give_log);
")


## ----measles_dmeasure, echo = FALSE, eval = TRUE,cache=F----------------------
measles_dmeasure <- spatPomp_Csnippet(
  unit_statenames = c("C"),
  unit_obsnames = c("cases"),
  code = "
  double m,v;
  double tol = 1e-300;
  double mytol = 1e-5;
  int u;

  lik= 0;
  for (u = 0; u < U; u++) {
    m = rho*(C[u]+mytol);
    v = m*(1.0-rho+psi*psi*m);
    if (cases[u] > 0.0) {
      lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
    } else {
        lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
    }
  }
  if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
  "
)

## ----measles_runit_measure, echo = TRUE, eval = T,cache=F---------------------
measles_runit_measure <- spatPomp_Csnippet("
  double cases;
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  cases = rnorm(m,sqrt(v));
  if (cases > 0.0) cases = nearbyint(cases);
  else cases = 0.0;
")

## ----measles_rmeasure, echo = FALSE, eval = TRUE,cache=F----------------------
measles_rmeasure <- spatPomp_Csnippet("
  const double *C = &C1;
  double *cases = &cases1;
  double m,v;
  int u;
  for (u= 0; u < U; u++) {
    m = rho*C[u];
    v = m*(1.0-rho+psi*psi*m);
    cases[u] = rnorm(m,sqrt(v));
    if (cases[u] > 0.0) {
      cases[u] = nearbyint(cases[u]);
    } else {
      cases[u] = 0.0;
    }
  }
")


## ----measles_emeasure, echo = TRUE, eval = T,cache=F--------------------------
measles_eunit_measure <- spatPomp_Csnippet("ey = rho*C;")
measles_vunit_measure <- spatPomp_Csnippet("
  double m = rho*C;
  vc = m*(1.0-rho+psi*psi*m);
")


## ----meales_sim_par,echo=T,eval=T---------------------------------------------
IVPs <- rep(c(0.032,0.00005,0.00004,0.96791),each=6) 
names(IVPs) <- paste0(rep(c('S','E','I','R'),each=6),1:6,"_0")
measles_params <- c(R0=30,A=0.5,muEI=52,muIR=52,muD=0.02,
  alpha=1,sigmaSE=0.01,rho=0.5,psi=0.1,g=1500,IVPs)


## ----measles_full_spatpomp, echo=T, eval=T,cache=F----------------------------
measles6 <- spatPomp(
  data = measles6,
  covar = measles_covar,
  unit_statenames = c('S','E','I','R','C'),
  unit_accumvars = c('C'),
  paramnames = names(measles_params),
  rinit = measles_rinit,
  rprocess = euler(measles_rprocess, delta.t=1/365),
  dunit_measure = measles_dunit_measure,
  eunit_measure = measles_eunit_measure,
  vunit_measure = measles_vunit_measure,
  runit_measure = measles_runit_measure,
  dmeasure = measles_dmeasure,
  rmeasure = measles_rmeasure,
  globals = measles_globals
)


## ----measles_seed,echo=F,eval=T-----------------------------------------------
set.seed(8)
measles_sim <- simulate(measles6,params=measles_params)


## ----measles_sim_plot,fig.width=7, fig.height=3.5,out.width = '90%', fig.align='center',fig.show='hold',echo=F,eval=T,cache=F,fig.cap=r"{A: reported measles cases in two week intervals for the six largest cities in England, \code{plot(measles6,log=TRUE)}. B: simulated data, \code{plot(simulate(measles6),log=TRUE)}. The vertical scale is \code{log10(cases+1)}.}"----

measles_plot_data <- plot(measles6,log=TRUE) +
  labs(y="") +
  theme(axis.text.y=element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ylim(0,4) +
  scale_x_continuous(breaks = c(1950,1960))

measles_plot_sim <- plot(measles_sim,log=TRUE) +
  theme(axis.title.y=element_blank(),
    plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  ylim(0,4) +
  scale_x_continuous(breaks = c(1950,1960))

prow <- plot_grid(
  measles_plot_data + theme(legend.position="none",
    axis.text.y = element_text(size = 8)),
  NULL,
  measles_plot_sim + theme(legend.position="none",
    axis.text.y = element_text(size = 8)),
  labels = c("A", "", "B"),
  nrow = 1,
  rel_widths = c(2,0.15,2),
  vjust = c(1.5,NA,1.5),
  hjust = c(-0.2,NA,0.8)
) + theme(plot.margin=margin(0,0,0,0))
prow

