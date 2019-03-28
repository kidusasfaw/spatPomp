##' Plotting
##'
##' Diagnostic plots.
##'
##' @name plot
##' @rdname plot
##' @include spatpomp_class.R
##' @include igirf.R
##' @aliases plot
##'
NULL

setGeneric(
  "plot",
  function (x, y, ...)
    standardGeneric("plot")
)


##' @export
setMethod(
  "plot",
  signature=signature(x="igirfd_spatpomp"),
  definition=function (x, ...) {
    plot.df <- data.frame(x@traces)
    cn <- colnames(plot.df)
    plot.df <- cbind(c(seq_len(dim(plot.df)[1])), plot.df)
    names(plot.df) <- c("iteration", cn)
    to.gather <- colnames(plot.df)[2:length(colnames(plot.df))]
    to.plot <- plot.df %>% tidyr::gather_(key = "param", val = "value", to.gather) %>% .[-1,]
    ggplot2::ggplot(data = to.plot) +
      geom_line(mapping = aes(x = iteration, y = value)) +
      facet_wrap(~param, ncol = 3, scales = "free")
  }
)
