# spatPomp to-do list

15. 23-06-08. Extend pomp-style list classes to spatPomp.

14. 23-06-08. Some development of diagnostic plots at https://kidusasfaw.github.io/spatPomp/vignettes/ibpf.pdf should be incorporated into the package. Related to [15].

13. 23-05-26. munit_measure does not yet have a Csnippet argument for unit-specific variables of the form Mtheta.

12. 23-05-24. bpfilter has log-likelihood NaN when all particles have zero likelihood. Better to return -Inf.

11. 23-05-18. Currently, spatPomp inference methodologies are not set up with a "continue" method. This should be added. As part of this, test that paramMatrix is treated consistently: iterated, perturbed filter methods should have a paramMatrix slot in the corresponding class, which is used by continue. Aim to follow the choices of pomp::mif2, where mif2_internal has a .paramMatrix argument used for continue applied to a pfilterd_pomp, but otherwise supplying a matrix of starting values is not supported---in principle, this variation can be controlled by rw_sd.

10. 23-04-17 applying mif2 to a spatPomp gives a pomp (more specifically, a mif2d_pomp). This is annoying if one uses mif2 on spatPomps. Any useful pomp function could have a wrapper so that it can be applied to a spatPomp while preserving the class, but preserving pomp class functionality within spatPomp is nontrivial. 

9. 23-03-29. Allow spatPomp to use the "order" pomp feature for covariate interpolation. Perhaps have a spatPomp_covariate_table function analogous to pomp::covariate_table?

8. 22-12-06. rmeasure and dmeasure could be built using runit_measure and dunit_measure if those are supplied.

7. 22-12-06. note that eunit_measure, vunit_measure and munit_measure calculate for all units and subset out the requested ones. this could lead to inefficient code when not used in a vectorized way.

6. 22-12-06. dunit_measure and runit_measure could be vectorized to prevent the need for vec_dmeasure and vec_rmeasure.

5. 22-09-02. add more unit tests (current coverage 83.1%)

4. 22-09-02. add localization to enkf and ienkf

3. 22-09-02. girf updates: looks for ways to get it working better on the measles example, for example Joonha Park's ideas:

(a). Since the population dynamics show exponential growth or decay, the pseudo-simulations (obtained by adding suitable noises to a deterministic skeleton) can be made on the log scale.

(b). The current bootstrap guide method scales the prediction variance linearly with the length of the prediction time interval for making pseudo-simulations.  However, this would not give accurate predictions when the size of the process noise varies over time.  A simple way to address this issue would be to make pseudo-prediction from time t_c to t in the following way:

pseudo-simulation at time t = (deterministic skeleton prediction at t) 
  + { (reference random prediction at time t)
     - (reference deterministic prediction at time t) }
  - { (reference random prediction at time t_c)
     - (reference deterministic prediction at time t_c) }

2. 20-11-05. handle multivariate observations

1. 20-11-05. output filter particles


