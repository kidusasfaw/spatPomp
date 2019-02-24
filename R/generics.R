## Independent Island Filter (IIF)
setGeneric("iif",function(object,...)standardGeneric("iif"))

## Extract units from spatpomp object
setGeneric("unit", function(x,...)standardGeneric("unit"))

## Extract the unit index from spatpomp object
setGeneric("unit_ix", function(x,...)standardGeneric("unit_ix"))

## Evaluate unit_dmeasure over all units
setGeneric("vec_dmeasure", function(object,...)standardGeneric("vec_dmeasure"))

## Adapted simulation for HIPPIE
setGeneric("hippie_pfilter", function(object,...)standardGeneric("hippie_pfilter"))

## HIPPIE
setGeneric("hippie", function(object,...)standardGeneric("hippie"))

