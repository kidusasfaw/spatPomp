# spatPomp to-do list

23-05-24 bpfilter has log-likelihood NaN when all particles have zero likelihood. Better to return -Inf.

23-05-18 Currently, spatPomp inference methodologies are not set up with a "continue" method. This should be added. As part of this, test that paramMatrix is treated consistently: iterated, perturbed filter methods should have a paramMatrix slot in the corresponding class, which is used by continue. Aim to follow the choices of pomp::mif2, where mif2_internal has a .paramMatrix argument used for continue applied to a pfilterd_pomp, but otherwise supplying a matrix of starting values is not supported---in principle, this variation can be controlled by rw_sd.

23-05-17 added method argument to spatPomp_Csnippet to allow it to figure out whether the variable needs a "const" specificiation. This is done in a backwards-compatible way - at some point, it could give warnings if this is not specified when it is needed.

23-04-17 applying mif2 to a spatPomp gives a pomp (more specifically, a mif2d_pomp). This is annoying if one uses mif2 on spatPomps. Any useful pomp function could have a wrapper so that it can be applied to a spatPomp while preserving the class, but preserving pomp class functionality within spatPomp is nontrivial. 

23-03-29. Allow spatPomp to use the "order" pomp feature for covariate interpolation. Perhaps have a spatPomp_covariate_table function analogous to pomp::covariate_table?

22-12-06. rmeasure and dmeasure could be built using runit_measure and dunit_measure if those are supplied.

22-12-06. note that eunit_measure, vunit_measure and munit_measure calculate for all units and subset out the requested ones. this could lead to inefficient code when not used in a vectorized way.

22-12-06. dunit_measure and runit_measure could be vectorized to prevent the need for vec_dmeasure and vec_rmeasure.

22-09-02. add more unit tests (current coverage 83.1%)

22-09-02. add localization to enkf and ienkf

22-09-02. girf updates: looks for ways to get it working better on the measles example, for example Joonha Park's ideas:

1. Since the population dynamics show exponential growth or decay, the pseudo-simulations (obtained by adding suitable noises to a deterministic skeleton) can be made on the log scale.

2. The current bootstrap guide method scales the prediction variance linearly with the length of the prediction time interval for making pseudo-simulations.  However, this would not give accurate predictions when the size of the process noise varies over time.  A simple way to address this issue would be to make pseudo-prediction from time t_c to t in the following way:

pseudo-simulation at time t = (deterministic skeleton prediction at t) 
  + { (reference random prediction at time t)
     - (reference deterministic prediction at time t) }
  - { (reference random prediction at time t_c)
     - (reference deterministic prediction at time t_c) }

20-11-05. handle multivariate observations

20-11-05. output filter particles


