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

##' @export
setMethod(
  "plot",
  signature=signature(x="spatPomp"),
  definition=function (x, log=F, ...) {
    df <- as.data.frame(x)
    if(log) df[x@unit_obsnames] <- log(df[x@unit_obsnames]+1)
    ggplot(data = df,
           mapping = aes(x = !!rlang::sym(x@timename),
                         y = factor(!!rlang::sym(x@unitname), level = unit_names(x)))) +
      geom_tile(mapping = aes(fill = !!rlang::sym(x@unit_obsnames))) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      theme(axis.text.y = element_text(size = 11-(2*floor(length(unit_names(x))/10)),
                                       vjust = 0.5,
                                       hjust=1)) +
      scale_fill_gradientn(colours = terrain.colors(10)) +
      labs(x = "time", y = "unit", fill = ifelse(log,
                                                 paste("log(",x@unit_obsnames,"+1)",sep=""),
                                                 x@unit_obsnames))
  }
)

