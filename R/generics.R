## Independent Island Filter (IIF)
setGeneric("iif",function(object,...)standardGeneric("iif"))

## Extract units from spatpomp object
setGeneric("unit", function(x,...)standardGeneric("unit"))

## Extract the unit index from spatpomp object
setGeneric("unit_ix", function(x,...)standardGeneric("unit_ix"))

## Evaluate unit_dmeasure over all units
# setGeneric("vec_dmeasure", function(x,...)standardGeneric("vec_dmeasure"))

## basic SMC (particle filter) FOR TEST PURPOSES ONLY
setGeneric("pfilter3",function(object,...)standardGeneric("pfilter3"))

## FOR TEST PURPOSES ONLY
setGeneric("dmeasure3",function(object,...)standardGeneric("dmeasure3"))
