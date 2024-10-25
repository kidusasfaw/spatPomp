---
title: 'spatPomp: An R package for spatiotemporal partially observed Markov process models'
tags:
  - R
  - statistics
  - POMP models
  - sequential Monte Carlo
  - particle filter
authors:
  - name: Kidus Asfaw
    affiliation: "1"
    orcid: 0000-0001-8625-3743
  - name: Joonha Park
    affiliation: "2"
    orcid: 0000-0002-4493-7730
  - name: Aaron A. King
    affiliation: "3, 4"
    orcid: 0000-0001-6159-3207
  - name: Edward L. Ionides
    corresponding: true 
    orcid: 0000-0002-4190-0174
    affiliation: "1"
affiliations:
 - name: University of Michigan, Department of Statistics.
   index: 1
 - name: University of Kansas, Department of Mathematics.
   index: 2
 - name: University of Michigan, Departments of Ecology & Evolutionary Biology and Complex Systems.
   index: 3
 - name: Santa Fe Institute, Santa Fe, New Mexico.
   index: 4
date: 3 June 2024
bibliography: paper.bib

---

# Summary

The development of `spatPomp` was motivated by the goal of investigating dynamics arising from a collection of spatially distributed, interacting biological populations.
The entire population, consisting of the union of these sub-populations over all the spatial locations, is called a metapopulation.
Each sub-population may have its own structure, which could correspond to disease status in an epidemiological model or abundance of several species in an ecosystem model. 
The `spatPomp` package embeds this goal in a more general problem: inference for spatiotemporal partially observed Markov process (SpatPOMP) models.
A POMP model consists of a latent Markov process model, together with a measurement model describing how the data arise from noisy and/or incomplete observation of this latent state.
The latent Markov process may be constructed in discrete or continuous time, taking scalar or vector values in a discrete or continuous space.
POMP models are also known as state space models, or hidden Markov models.
A SpatPOMP model extends the POMP model formulation by adding an index set corresponding to spatial location, so that the state of the SpatPOMP is comprised of a value for each location.
We say "unit" rather than "spatial location" to build our framework in the general context of an arbitrary index set. 
Measurements are made on each unit, and are assumed to depend only on the latent state value for that unit.
The `spatPomp` R package provides a computational framework for modeling and statistical inference on SpatPOMP models.

# Statement of Need

The `spatPomp` package provides statistical methodology for a broad class of nonlinear and non-Gaussian SpatPOMP models.
This gives scientists the freedom to construct and analyze scientifically motivated mechanistic models.
`spatPomp` emphasizes likelihood-based inference, using scalable Monte Carlo methods to evaluate and maximize the likelihood function.
Previous approaches for evaluating the likelihood function for SpatPOMP models required specific model assumptions: linear Gaussian SpatPOMP models can be investigated using the Kalman filter [@kalman60]; SpatPOMP models with sufficiently minor deviations from linearity and Gaussianity  can be effectively analyzed using the extended Kalman filter or the ensemble Kalman filter [@evensen22].
Likelihood evaluation for highly nonlinear low-dimensional POMP models can be carried using the particle filter, also known as sequential Monte Carlo [@chopin20].
However, the particle filter suffers from a curse of dimensionality that makes it inapplicable on SpatPOMP models.
Recent algorithmic developments have addressed this limitation, permitting consideration of the general class of nonlinear non-Gaussian SpatPOMP models.
`spatPomp` provides implementations of such algorithms, including bagged particle filters [@ionides23], block particle filters [@ning23], guided particle filters [@park20], and ensemble Kalman filters [@evensen22].

SpatPOMP models with high nonlinearity and stochasticity can arise when investigating the ecological dynamics of a spatially distributed collection of interacting biological populations, known as a metapopulation.
Existing demonstrations of `spatPomp` have arisen from studying the ecology of infectious diseases [@wheeler24;@li24;@zhang22].
In such epidemiological settings, the state of each unit may be comprised of the abundance of a pathogen species, a host species, and perhaps also a vector species.
We anticipate that epidemiology will continue to be a major application area for SpatPOMP models.
However, this is a general model class with potential applications across the biological and social sciences, healthcare, engineering, industry and government. 

The `spatPomp` package is designed for researchers who aim to develop scientifically plausible dynamic models to describe spatiotemporal systems.
The package assists with the application of existing models, modification of such models, or the development of entirely new models.
It provides methodologies to carry out statistical inference on these models, involving parameter estimation, model selection, and model criticism.
It focuses on algorithms with the plug-and-play property, meaning that the dynamic model can be specified by code to simulate the latent process for this model.
A consequence of the plug-and-play property is that the data analyst is not required to provide explicit specification of transition probabilities.
This makes `spatPomp` a flexible tool to assist model development.

The `spatPomp` package builds on `pomp` [@king16] which is a successful software package for low-dimensional POMP models.
Other packages with similar capabilities to `pomp` include `nimble` [@michaud21], `LiBBi` [@murray15] and `mcstate` with `odin` and `dust` [@fitzjohn20].
All these packages enable plug-and-play inference based on sequential Monte~Carlo.
Markov chain Monte Carlo packages, such as `stan`, have been found to be effective for inference on some POMP models [@li18] but they lack the plug-and-play property.
Perhaps for that reason, sequential Monte Carlo methods have found broader applicability for this model class.
We are not aware of alternative packages to `spatPomp` that provide statistically efficient, plug-and-play inference for the general class of SpatPOMP models.


# Package Design

The `spatPomp` package provides a standardized interface between SpatPOMP models and statistical inference methods.
This approach is designed to provide an environment for data analysis using existing algorithms as well as the development of new algorithms.
New methods can readily be tested on existing models, since the models have defined operations (such as simulation, or evaluation of the measurement density) that the methods can access.
For the same reason, new models can readily be investigated using a range of methods.
Currently, all the methods in `spatPomp` have the plug-and-play property, i.e., they require a simulator for the SpatPOMP model under investigation, but not an evaluator of its transitions densities.
The functionality of `spatPomp` permits specification of transition probabilities, so it is possible to implement algorithms without the plug-and-play property.
However, based on the development trajectory of `pomp` [@king16], we anticipate that most use of `spatPomp` will focus on plug-and-play methods.

The development of `spatPomp` has focused on likelihood-based inference.
However, the framework also permits Bayesian inference and consideration of non-likelihood-based model fitting criteria.

A distance may be defined between units, and algorithms may assume that distant units have only weak interactions.
Such assumptions may involve a bias/variance tradeoff specific to the choice of model and the choice of inference algorithm.
Therefore, it may be beneficial to evaluate various different algorithms when investigating a specific model of scientific interest.
The inter-operability of methods across models, provided by `spatPomp`, facilitates consideration of a range of methods.

# Resources

The `spatPomp` website ([https://kidusasfaw.github.io/spatPomp/](https://kidusasfaw.github.io/spatPomp)) provides links to various resources for users and developers of the package. This includes the following.

1. An extended tutorial [@asfaw24] introduces the mathematical framework behind `spatPomp`, describes the software implementation of this framework, provides pseudocode for various algorithms included in the package, and illustrates some basic usage.
Section 2 explains the elementary methods used to access properties of the `spatPomp` model class.
These elementary methods are the building blocks available to developers for implementing complex algorithms acting on spatPomp models.

2. A tutorial provided as a supplement to [@ning23] focuses specifically on the iterated block particle filter algorithm. This is available at [https://kidusasfaw.github.io/spatPomp/vignettes/ibpf.pdf](https://kidusasfaw.github.io/spatPomp/vignettes/ibpf.pdf).

3. A numerical comparison of spatiotemporal filtering methods by @ionides23, carried out using `spatPomp`,  has source code available at [https://github.com/ionides/bagged_filters](https://github.com/ionides/bagged_filters).

4. A spatiotemporal data analysis of cholera transmission in Haiti [@wheeler24], carried out using `spatPomp`,  has source code available at [https://github.com/jeswheel/haiti_article](https://github.com/jeswheel/haiti_article).

5. A spatiotemporal data analysis of COVID-19 transmission in China [@li24], carried out using `spatPomp`,  has source code available at [https://github.com/jifanli/metapop_article](https://github.com/jifanli/metapop_article).


# Acknowledgments

This work was supported by National Science Foundation grants DMS-1761603 and DMS-1646108, and National Institutes of Health grants 1-U54-GM111274, 1-U01-GM110712 and 1-R01-AI143852.
We recognize those who have participated in the development and testing of `spatPomp`, especially Allister Ho, Zhuoxun Jiang, Jifan Li, Patricia Ning, Eduardo Ochoa, Rahul Subramanian and Jesse Wheeler.
We are grateful to John Lees, Ben Bolker, and the editors at The Journal of Open Source Software for their constructive feedback.

# References
