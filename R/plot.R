##' Plotting \code{spatPomp} data
##'
##' Visualize \code{spatPomp} data
##' @name plot
##' @rdname plot
##' @include spatPomp_class.R
##' @param x a \code{spatPomp} object
##' @param log should the data be log-transformed before plotting?
##' This helps in contexts where there are spikes that could take away
##' attention from the dynamics illustrated by the rest of the data.
##' @param ncol the number of columns in the grid plot
NULL

##' Plotting output of \code{igirf()}
##'
##' Diagnostic plot for \code{igirf()}
##' @param params the names of the parameters for which the user would like to see a trace plot
##' @name plot-igirfd_spatPomp
##' @rdname plot
##' @aliases plot,igirfd_spatPomp-method
setMethod(
  "plot",
  signature=signature(x="igirfd_spatPomp"),
  definition=function (x, params = names(coef(x)), ncol = 3) {
    if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
    plot.df <- data.frame(x@traces[,c("loglik", params)])
    cn <- colnames(plot.df)
    plot.df <- cbind(c(seq_len(dim(plot.df)[1])), plot.df)
    names(plot.df) <- c("iteration", cn)
    to.gather <- colnames(plot.df)[2:length(colnames(plot.df))]
    to.plot <- plot.df %>% tidyr::gather_(key = "param", val = "value", to.gather) %>% tail(-1)
    ggplot2::ggplot(data = to.plot) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = .data$iteration, y = .data$value)) +
      ggplot2::facet_wrap(~param, ncol = ncol, scales = "free")
  }
)

##' Plotting \code{spatPomp} data
##'
##' Visualize \code{spatPomp} data
##' @name plot-spatPomp
##' @rdname plot
##' @aliases plot,spatPomp-method
##' @export
setMethod(
  "plot",
  signature=signature(x="spatPomp"),
  definition=function (x, log=F) {
    df <- as.data.frame(x)
    if(log) df[x@unit_obsnames] <- log(df[x@unit_obsnames]+1)
    ggplot2::ggplot(data = df,
           mapping = ggplot2::aes(x = !!rlang::sym(x@timename),
                                  y = factor(!!rlang::sym(x@unitname),
                                             levels = unit_names(x)))) +
      ggplot2::geom_tile(mapping = ggplot2::aes(fill = !!rlang::sym(x@unit_obsnames))) +
      ggplot2::scale_x_continuous(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(
        size = 11-(2*floor(length(unit_names(x))/10)),
        vjust = 0.5,
        hjust=1)) +
      ggplot2::scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
      ggplot2::labs(x = "time", y = "unit", fill = ifelse(log,
                                                 paste("log(",x@unit_obsnames,"+1)",sep=""),
                                                 x@unit_obsnames))
  }
)

