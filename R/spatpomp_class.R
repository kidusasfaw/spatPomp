## define the pomp class
setClass(
  'spatpomp',
  contains="pomp",
  slots=c(
    units = 'character',
    unit_index = 'character',
    unit_statenames = 'character',
    global_statenames = 'character',
    unit_dmeasure = 'function',
    covar = 'array'
  ),
  prototype=prototype(
    unit_dmeasure = function(y,x,times,params,d,log=FALSE,...)
      stop(sQuote("unit_dmeasure")," not specified")
  )
)
