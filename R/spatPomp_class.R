#' An S4 class to represent a spatiotemporal POMP model and data.
#'
#' @slot units A vector containing the spatial units of a spatiotemporal POMP.
#' @slot unit_index A vector with the spatial units as values and the indices for the units as the names.
#' @slot unit_statenames A vector containing the state names such that appending the unit indices to the
#' unit statenames will result in the each unit's corresponding states.
#' @slot global_statenames A vector containing the statenames shared by all spatial units.
#' @slot obstypes A vector of observation types for a spatial unit.
#' @slot emeasure A pomp_fun representing the expected measurement for each spatial unit given its states.
#' @slot unit_dmeasure A pomp_fun representing the unit measurement density for each spatial unit.
#' @slot unit_rmeasure A pomp_fun representing the unit observation simulator.
setClass(
  'spatPomp',
  contains="pomp",
  slots=c(
    units = 'character',
    unit_index = 'character',
    unit_statenames = 'character',
    global_statenames = 'character',
    obstypes = 'character',
    emeasure = 'pomp_fun',
    mmeasure = 'pomp_fun',
    vmeasure = 'pomp_fun',
    unit_dmeasure = 'pomp_fun',
    unit_rmeasure = 'pomp_fun'
  ),
  prototype=prototype(
    emeasure = pomp:::pomp_fun(slotname="emeasure"),
    mmeasure = pomp:::pomp_fun(slotname="mmeasure"),
    vmeasure = pomp:::pomp_fun(slotname="vmeasure"),
    unit_dmeasure = pomp:::pomp_fun(slotname="unit_dmeasure"),
    unit_rmeasure = pomp:::pomp_fun(slotname="unit_rmeasure")
  )
)
