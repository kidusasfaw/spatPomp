##' Plotting
##'
##' Diagnostic plots.
##'
##' @name plot
##' @rdname plot
##' @include spatPomp_class.R
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
  signature=signature(x="igirfd_spatPomp"),
  definition=function (x, ...) {
    plot.df <- data.frame(x@traces)
    cn <- colnames(plot.df)
    plot.df <- cbind(c(seq_len(dim(plot.df)[1])), plot.df)
    names(plot.df) <- c("iteration", cn)
    to.gather <- colnames(plot.df)[2:length(colnames(plot.df))]
    to.plot <- plot.df %>% tidyr::gather_(key = "param", val = "value", to.gather) %>% .[-1,]
    ggplot2::ggplot(data = to.plot) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = iteration, y = value)) +
      ggplot2::facet_wrap(~param, ncol = 3, scales = "free")
  }
)
