#' An S4 class to represent a spatiotemporal POMP model and data.
#'
#' @slot unit_names A vector containing the spatial units of a spatiotemporal POMP.
#' @slot unit_statenames A vector containing the state names such that appending the unit indices to the
#' unit statenames will result in the each unit's corresponding states.
#' @slot unit_obsnames A vector of observation types for a spatial unit.
#' @slot eunit_measure A pomp_fun representing the expected measurement for each spatial unit given its states.
#' @slot dunit_measure A pomp_fun representing the unit measurement density for each spatial unit.
#' @slot runit_measure A pomp_fun representing the unit observation simulator.
setClass(
  'spatPomp',
  contains="pomp",
  slots=c(
    unit_names = 'character',
    unit_statenames = 'character',
    unit_obsnames = 'character',
    unitname = 'character',
    unit_covarnames = 'character',
    shared_covarnames = 'character',
    unit_accumvars = 'character',
    eunit_measure = 'pomp_fun',
    munit_measure = 'pomp_fun',
    vunit_measure = 'pomp_fun',
    dunit_measure = 'pomp_fun',
    runit_measure = 'pomp_fun'
  ),
  prototype=prototype(
    unit_names = as.character(NA),
    unit_statenames = as.character(NA),
    unit_obsnames = as.character(NA),
    unitname = as.character(NA),
    unit_covarnames = as.character(NA),
    shared_covarnames = as.character(NA),
    unit_accumvars = as.character(NA),
    eunit_measure = pomp:::pomp_fun(slotname="eunit_measure"),
    munit_measure = pomp:::pomp_fun(slotname="munit_measure"),
    vunit_measure = pomp:::pomp_fun(slotname="vunit_measure"),
    dunit_measure = pomp:::pomp_fun(slotname="dunit_measure"),
    runit_measure = pomp:::pomp_fun(slotname="runit_measure")
  )
)
