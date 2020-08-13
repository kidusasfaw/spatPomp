#' An S4 class to represent a spatiotemporal POMP model and data.
#'
#' @slot units A vector containing the spatial units of a spatiotemporal POMP.
#' @slot unit_statenames A vector containing the state names such that appending the unit indices to the
#' unit statenames will result in the each unit's corresponding states.
#' @slot obstypes A vector of observation types for a spatial unit.
#' @slot emeasure A pomp_fun representing the expected measurement for each spatial unit given its states.
#' @slot unit_dmeasure A pomp_fun representing the unit measurement density for each spatial unit.
#' @slot unit_rmeasure A pomp_fun representing the unit observation simulator.
setClass(
  'spatPomp',
  contains="pomp",
  slots=c(
    units = 'character',
    unit_statenames = 'character',
    obstypes = 'character',
    unitname = 'character',
    unit_covarnames = 'character',
    shared_covarnames = 'character',
    unit_emeasure = 'pomp_fun',
    unit_mmeasure = 'pomp_fun',
    unit_vmeasure = 'pomp_fun',
    unit_dmeasure = 'pomp_fun',
    unit_rmeasure = 'pomp_fun'
  ),
  prototype=prototype(
    unitname = as.character(NA),
    unit_covarnames = as.character(NA),
    shared_covarnames = as.character(NA),
    unit_emeasure = pomp:::pomp_fun(slotname="emeasure"),
    unit_mmeasure = pomp:::pomp_fun(slotname="mmeasure"),
    unit_vmeasure = pomp:::pomp_fun(slotname="vmeasure"),
    unit_dmeasure = pomp:::pomp_fun(slotname="unit_dmeasure"),
    unit_rmeasure = pomp:::pomp_fun(slotname="unit_rmeasure")
  )
)
