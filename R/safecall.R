safecall <- function (...) {
  new("safecall",call=match.call(),envir=parent.frame())
}
