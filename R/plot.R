##' Plot methods for \code{spatPomp} objects
##'
##' Visualize data in a \code{spatPomp} object or a derived class.
##' This gives a quick view; the data can be extracted from 
##' the object to make a customized plot.
##'
##' @name plot
##' @rdname plot
##' @include spatPomp_class.R
##' @param x a \code{spatPomp} object
##' @param log should the data be transformed to \code{log10(x+1)} before plotting?
##' This helps in contexts where there are spikes that could take away
##' attention from the dynamics illustrated by the rest of the data. 
##' @param ncol the number of columns in the grid plot
##' @param type for visualizing an object of class \code{spatPomp}, the user
##' can obtain a grid of line plots by default (\code{l}) or a heat map by
##' supplying argument \code{h}.
##' @param plot_unit_names allows suppression of unit names when making
##' a heat map for a large number of units
##' @param params allows selection of a subset of parameters when making
##' a diagnostic plot for a model with many parameters
##' @param ...  for visualizing an object of class \code{spatPomp}, the user
##' can add arguments like \code{nrow} to specify the number of rows in the
##' grid.
##' @return a \code{ggplot} plot of class \sQuote{gg} and \sQuote{ggplot}.
##' @importFrom ggplot2 ggplot geom_line aes facet_wrap
##' @importFrom graphics par
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

##' @name plot-spatPomp
##' @rdname plot
##' @aliases plot,spatPomp-method
##' @export
setMethod(
  "plot",
  signature=signature(x="spatPomp"),
  definition=function (x, type = c('l','h'), log=FALSE, plot_unit_names=TRUE,...) {
    df <- as.data.frame(x)
    type <- match.arg(type)
    if(type == 'l'){
      if(log) df[x@unit_obsnames] <- df[x@unit_obsnames]+1
      unit_nm <- rlang::sym(x@unitname)
      df[[unit_nm]] <- factor(df[[unit_nm]], levels = x@unit_names)
      g <- ggplot2::ggplot(data = df,
        mapping = ggplot2::aes(
	  x = !!rlang::sym(x@timename),
          y = !!rlang::sym(x@unit_obsnames)
	)
      ) +
      ggplot2::labs(y=paste0(x@unit_obsnames,"+1")) +
      ggplot2::scale_y_continuous(transform=if(log) "log10" else "identity") +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(unit_nm, ...)
      return(g)
    }
    if(type == 'h'){
      if(log) df[x@unit_obsnames] <- log10(df[x@unit_obsnames]+1)
      g <- ggplot2::ggplot(data = df,
        mapping = ggplot2::aes(x = !!rlang::sym(x@timename),
          y = factor(!!rlang::sym(x@unitname),
          levels = unit_names(x)))) +
      ggplot2::geom_tile(mapping = ggplot2::aes(fill = !!rlang::sym(x@unit_obsnames))) +
      ggplot2::scale_x_continuous(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::scale_fill_gradient(low = "#000000", high = "#FFFFFF") +
      ggplot2::labs(x = "time",
        y = "unit",
        fill = ifelse(log,
        paste("log10\n(",x@unit_obsnames,"+1)",sep=""),x@unit_obsnames)
      )
      if(plot_unit_names) {
        g <- g + ggplot2::theme(
	  axis.text.y = ggplot2::element_text(
            size = 11-(2*floor(length(unit_names(x))/10)),
            vjust = 0.5, hjust=1))
      } else {
        g <- g + ggplot2::theme(
	  axis.text.y=ggplot2::element_blank(),
	  axis.ticks.y=ggplot2::element_blank())
      }
      return(g)
    }
  }
)

##' @name plot-igirfd_spatPomp
##' @rdname plot
##' @aliases plot,igirfd_spatPomp-method
setMethod(
  "plot",
  signature=signature(x="igirfd_spatPomp"),
  definition=function (x, params = names(coef(x)), ncol = 3) {
    plot.df <- data.frame(x@traces[,c("loglik", params)])
    cn <- colnames(plot.df)
    plot.df <- cbind(c(seq_len(dim(plot.df)[1])), plot.df)
    names(plot.df) <- c("iteration", cn)
    to.gather <- colnames(plot.df)[2:length(colnames(plot.df))]
    to.plot <- plot.df |> tidyr::gather(key = "param", val = "value", tidyr::all_of(to.gather)) |> tail(-1)
    ggplot2::ggplot(data = to.plot) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = .data$iteration, y = .data$value)) +
      ggplot2::facet_wrap(~param, ncol = ncol, scales = "free")
  }
)
