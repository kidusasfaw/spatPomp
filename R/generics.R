## Adapted Simulation Island Filter (ASIF)
setGeneric("asif",function(object,...)standardGeneric("asif"))

## Adapted Simulation Island Filter (ASIF)
setGeneric("asifir",function(object,...)standardGeneric("asifir"))

## Extract units from spatPomp object
setGeneric("spat_units", function(x,...)standardGeneric("spat_units"))

## Evaluate unit_dmeasure over all units
setGeneric("vec_dmeasure", function(object,...)standardGeneric("vec_dmeasure"))

## Simulate observations over all units given states
setGeneric("vec_rmeasure", function(object,...)standardGeneric("vec_rmeasure"))

## Csnippet utility that frees the user from declaring shorthand pointers
setGeneric("spatPomp_Csnippet", function(object,...)standardGeneric("spatPomp_Csnippet"))

## Adapted simulation for HIPPIE
setGeneric("hippie_pfilter", function(object,...)standardGeneric("hippie_pfilter"))

## HIPPIE
setGeneric("hippie", function(object,...)standardGeneric("hippie"))

## EnKF
setGeneric("genkf",  function (data, ...)standardGeneric("genkf"))

## Block Particle Filter
setGeneric("bpfilter",  function (object, ...)standardGeneric("bpfilter"))


## construct_spatPomp
setGeneric(
  "construct_spatPomp",
  function (data, times, units, ...)
    standardGeneric("construct_spatPomp")
)

