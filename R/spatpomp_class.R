## define the pomp class
setClass(
  'spatpomp',
  contains="pomp",
  slots=c(
    units = 'character',
    unit_index = 'character',
    unit_statenames = 'character',
    global_statenames = 'character',
    obstypes = 'character',
    unit_dmeasure = 'pomp_fun'
  ),
  prototype=prototype(
    unit_dmeasure = pomp2:::pomp_fun(slotname="unit_dmeasure")
  )
)
