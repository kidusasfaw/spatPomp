## define the pomp class
setClass(
  'spatpomp',
  contains="pomp",
  slots=c(
    units = 'character',
    unit_statenames = 'character',
    global_statenames = 'character',
    unit_dmeasure = 'function'
  ),
  prototype=prototype(
    unit_dmeasure = function(y,x,times,params,d,log=FALSE,...)
      stop(sQuote("unit_dmeasure")," not specified")
  )
)

setMethod(
  "show",
  signature=signature(object="spatpomp"),
  definition=function (object) {print("donezo")}
)
