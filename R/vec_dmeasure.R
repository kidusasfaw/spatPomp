## evaluate the unit measurement model density function for each unit

vec_dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  nunits <- length(object@units)
  statenames <- c(object@unit_statenames, object@global_statenames)
  storage.mode(params) <- "double"
  rv <- vector(length = nunits)
  for(i in 1:nunits){
    unit_y <- y[stringr::str_detect(rownames(y), paste0(i,"$")),,drop = FALSE]
    unit_x <- x[c(stringr::str_detect(rownames(x), paste0(i,"$")),object@global_statenames),,drop = FALSE]
    storage.mode(unit_y) <- "double"
    storage.mode(unit_x) <- "double"
    rv[i] <- .Call(do_unit_dmeasure,object,unit_y,unit_x,times,params,log,statenames,.getnativesymbolinfo)
  }
  pompUnload(object)
  rv
}

# setMethod(
#   "vec_dmeasure",
#   signature=signature(object="spatpomp"),
#   definition=function (object, y, x, times, params, log = FALSE, ...)
#     vec_dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
# )
